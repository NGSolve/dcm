#include <comp.hpp>
#include "l2cells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"

namespace ngcomp
{

  
  void L2CellTrig ::
  CalcShape (const IntegrationPoint & ip, 
             BareSliceVector<> shape) const 
  {
    double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
    int maxlam = PosMax(lam);
    
    shape.Range(GetNDof()) = 0;
    
    int minvi = (maxlam+1)%3;
    int maxvi = minvi;
    for (int i = 0; i < 3; i++)
      if (i != maxlam)
        {
          if (vnums[i] < vnums[minvi]) minvi = i;
          if (vnums[i] > vnums[maxvi]) maxvi = i;
        }

    int vdir[3];
    vdir[maxlam] = -1;
    vdir[minvi] = 0;
    vdir[maxvi] = 1;
    
    Vec<2> x(lam[minvi], lam[maxvi]);
    Vec<2> xi = MapTrig2Quad (x);

    Mat<2,2> F = DMapQuad2Trig(xi);
    Mat<2,2> F2;    // trafo from vertex permutation
    Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };
    F2.Col(0) = verts[minvi]-verts[maxlam];
    F2.Col(1) = verts[maxvi]-verts[maxlam];

    double trafo = 1./fabs(Det(F2*F));
    ArrayMem<double, 20> polxi(order+1), poleta(order+1);
    LagrangePolynomials(xi(0), GaussIR, polxi);
    LagrangePolynomials(xi(1), GaussIR, poleta);

    auto assign =  [&](int nr, IVec<2> ind)
    {
      shape(nr) = trafo*polxi[ind[0]]*poleta[ind[1]];
    };


    int ii = 0;

    // internal faces into direction of vertices v1, v2
    for (int v1 = 0; v1 < 3; v1++)
      for (int v2 = 0; v2 < v1; v2++)
        {
          if (v1 != maxlam && v2 != maxlam)
            {
              IVec<2> ind = { 0, 0 };
              int dirv1 = vdir[v1];
              int dirv2 = vdir[v2];
              for (int k = 0, kk=ii; k < order+1; k++)
                for (int l = 0; l < order+1; l++, kk++)
                  {
                    ind[dirv1] = k;
                    ind[dirv2] = l;
                    assign(kk, ind);
                  }
            }
          ii += sqr(order+1);
        }
  }


  void L2CellTrig ::
  Interpolate (const ElementTransformation & trafo, 
        const class CoefficientFunction & func, SliceMatrix<> coefs,
        LocalHeap & lh) const 
    {
      throw Exception ("L2CellTrig :: Interpolate not implemented");
      /*
      // auto [pnts, inds, nump] = GetMicroCellPoints<ET_TRIG> (order, true);

      IntegrationRule irtrig;
      for (int v = 0; v < 3; v++)
        for (auto ipx : GaussRadauIR)
          for (auto ipy : GaussRadauIR)
          {
            if (ipx(0) > 1-1e-13) continue;
            Vec<2> xi(ipx(0), ipy(0));
            Vec<2> x = MapQuad2Trig (xi, v);
            Mat<2,2> jac = DMapQuad2Trig(xi, v);
            double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
            irtrig += IntegrationPoint(x(0), x(1), 0, w);
          }
      irtrig += IntegrationPoint(1.0/3, 1.0/3, 0, 0);

      Matrix mat(this->ndof, this->ndof);
      Vector rhs(this->ndof);
      Vec<1> fi;
      for (int i = 0; i < this->ndof; i++)
      {
        // Vec<2> x = pnts[i];
        // IntegrationPoint ip(x(0), x(1), 0, 0);
        MappedIntegrationPoint<2,2> mip(irtrig[i], trafo);
        func.Evaluate (mip, fi);
        rhs(i) = fi(0);
        this->CalcShape(irtrig[i], mat.Row(i));
      }
      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
      */
    }          





  L2Cells::
  L2Cells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {
    switch (ma->GetDimension())
      {
      case 2:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdDual<2,2>>>();
          break;
        }
      }
    Array<double> xi, wi;
    ComputeGaussRule (order+1, xi, wi);
    for (auto i : Range(xi))
      GaussIR.Append (IntegrationPoint (1-xi[i]-1e-14, 0, 0, wi[i]));
  }

  void L2Cells :: Update()
  {
      FESpace::Update();
      size_t ndof = 0;
      
      first_element_dofs.SetSize(ma->GetNE(VOL)+1);
      first_element_dofs[0] = ndof;
      switch (ma->GetDimension())
        {
        case 2:
          {
            for (auto i : Range(ma->GetNE(VOL)))
              first_element_dofs[i+1] = ndof += 3*sqr(order+1);
            break;
          }
        default:
          throw Exception("invalid dimension");
        }
    SetNDof(ndof);
  }
  
  void L2Cells ::
    GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
      dnums.SetSize0();
      if (ei.VB() != VOL) return;
      dnums += Range(first_element_dofs[ei.Nr()], first_element_dofs[ei.Nr()+1]);
    }

  
  FiniteElement & L2Cells ::
    GetFE (ElementId ei, Allocator & alloc) const
  {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
      {
        case ET_TRIG:
          {
            auto trig = new (alloc) L2CellTrig(order, GaussIR);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }
        default:
          throw Exception("element not implemented");
      }
    }


  /*
  shared_ptr<BaseMatrix> H1PrimalCells ::
    GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
    {
      static Timer t("GetMassOperator"); RegionTimer reg(t);
      static Timer tint("integrate");
      static Timer tsp("MakeSparse");

      HeapReset hr(lh);
      int dim = ma->GetDimension();
      switch (dim)
      {
        case 2:
          {
            auto irs = this->GetIntegrationRules();
            IntegrationRule ir = std::move(irs[ET_TRIG]);

            auto & felref = dynamic_cast<const H1PrimalCellTrig&> (GetFE(ElementId(VOL,0), lh));
            Matrix shapes(felref.GetNDof(), ir.Size());
            Matrix shapes_trans(ir.Size(), felref.GetNDof());

            Array<int> rowind, colind;
            Array<double> values;

            Array<short> classnr(ma->GetNE(VOL));
            ma->IterateElements
              (VOL, lh, [&] (auto el, LocalHeap & llh)
               {
               classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
               });

            TableCreator<size_t> creator;
            for ( ; !creator.Done(); creator++)
              for (auto i : Range(classnr))
                creator.Add (classnr[i], i);
            Table<size_t> table = creator.MoveTable();


            for (auto elclass_inds : table)
            {
              if (elclass_inds.Size() == 0) continue;

              HeapReset hr(lh);


              size_t nr = elclass_inds[0];
              auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
              auto & felref = dynamic_cast<const H1PrimalCellTrig&> (GetFE(ElementId(VOL,nr), lh));

              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Col(i));
              for (int i = 0; i < shapes.Height(); i++)
                for (int j = 0; j < shapes.Width(); j++)
                  if (fabs(shapes(i,j)) < 1e-14)
                    shapes(i,j) = 0;
              SmallSparseMatrix spshapes(shapes, 1e-14);
              shapes_trans = Trans(shapes);

              tint.Start();

              MyMutex COOmutex;
              ParallelForRange (elclass_inds.Size(), [&](auto myrange)
                  {
                  auto mylh = lh.Split();
                  Array<int> myrowind, mycolind;
                  Array<double> myvalues;
                  Array<DofId> dofs(felref.GetNDof());
                  Matrix rhoi_shapes_trans(ir.Size(), felref.GetNDof());
                  Matrix<> elmat(felref.GetNDof());

                  for (auto i : myrange)
                  {
                  HeapReset hr(mylh);
                  auto nr = elclass_inds[i];

                  auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                  MappedIntegrationRule<2,2> mir(ir, trafo, mylh);
                  GetDofNrs(ElementId(VOL,nr), dofs);

                  for (size_t i = 0; i < mir.Size(); i++)
                  {
                  double rhoi = 1.;  //coefficient is always 1 right now
                  rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                  rhoi_shapes_trans.Row(i) = rhoi * shapes_trans.Row(i);
                  }

                  spshapes.Mult (rhoi_shapes_trans, elmat);

                  for (int i = 0; i < dofs.Size(); i++)
                    for (int j = 0; j < dofs.Size(); j++)
                      if (fabs(elmat(i,j)) > 1e-14)
                      {
                        myrowind.Append(dofs[i]);
                        mycolind.Append(dofs[j]);
                        myvalues.Append(elmat(i,j));
                      }
                  }

                  COOmutex.lock();
                  rowind += myrowind;
                  colind += mycolind;
                  values += myvalues;
                  COOmutex.unlock();
                  });

              tint.Stop();
            }    


            tsp.Start();
            auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
            tsp.Stop();

            return make_shared<MySuperSparseMatrix> (std::move(*spmat));
          }
        case 3:
          {
            auto irs = this->GetIntegrationRules();
            IntegrationRule ir = std::move(irs[ET_TET]);

            auto & felref = dynamic_cast<const H1PrimalCellTet&> (GetFE(ElementId(VOL,0), lh));
            Matrix shapes(felref.GetNDof(), ir.Size());
            Matrix shapes_trans(ir.Size(), felref.GetNDof());

            Array<int> rowind, colind;
            Array<double> values;

            Array<short> classnr(ma->GetNE(VOL));
            ma->IterateElements
              (VOL, lh, [&] (auto el, LocalHeap & llh)
               {
               classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
               });

            TableCreator<size_t> creator;
            for ( ; !creator.Done(); creator++)
              for (auto i : Range(classnr))
                creator.Add (classnr[i], i);
            Table<size_t> table = creator.MoveTable();


            for (auto elclass_inds : table)
            {
              if (elclass_inds.Size() == 0) continue;

              HeapReset hr(lh);


              size_t nr = elclass_inds[0];
              auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
              auto & felref = dynamic_cast<const H1PrimalCellTet&> (GetFE(ElementId(VOL,nr), lh));

              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Col(i));
              for (int i = 0; i < shapes.Height(); i++)
                for (int j = 0; j < shapes.Width(); j++)
                  if (fabs(shapes(i,j)) < 1e-8)
                    shapes(i,j) = 0;
              SmallSparseMatrix spshapes(shapes, 1e-8);
              shapes_trans = Trans(shapes);

              tint.Start();

              MyMutex COOmutex;
              ParallelForRange (elclass_inds.Size(), [&](auto myrange)
                  {
                  auto mylh = lh.Split();
                  Array<int> myrowind, mycolind;
                  Array<double> myvalues;
                  Array<DofId> dofs(felref.GetNDof());
                  Matrix rhoi_shapes_trans(ir.Size(), felref.GetNDof());
                  Matrix<> elmat(felref.GetNDof());

                  for (auto i : myrange)
                  {
                  HeapReset hr(mylh);
                  auto nr = elclass_inds[i];

                  auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                  MappedIntegrationRule<3,3> mir(ir, trafo, mylh);
                  GetDofNrs(ElementId(VOL,nr), dofs);

                  for (size_t i = 0; i < mir.Size(); i++)
                  {
                  double rhoi = 1.;  //coefficient is always 1 right now
                  rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                  rhoi_shapes_trans.Row(i) = rhoi * shapes_trans.Row(i);
                  }

                  spshapes.Mult (rhoi_shapes_trans, elmat);

                  for (int i = 0; i < dofs.Size(); i++)
                    for (int j = 0; j < dofs.Size(); j++)
                      if (fabs(elmat(i,j)) > 1e-14)
                      {
                        myrowind.Append(dofs[i]);
                        mycolind.Append(dofs[j]);
                        myvalues.Append(elmat(i,j));
                      }
                  }

                  COOmutex.lock();
                  rowind += myrowind;
                  colind += mycolind;
                  values += myvalues;
                  COOmutex.unlock();
                  });

              tint.Stop();
            }    


            tsp.Start();
            auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
            tsp.Stop();

            return make_shared<MySuperSparseMatrix> (std::move(*spmat));
          }
        default:
          throw Exception("MassOperator not implemented for dim!=2,3");
      }
    }
  */


#ifdef COMMENTED
  std::map<ELEMENT_TYPE, IntegrationRule> H1PrimalCells::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;

    IntegrationRule irtrig;
    for (int v = 0; v < 3; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
        {
          Vec<2> xi(ipx(0), ipy(0));
          Vec<2> x = MapQuad2Trig (xi, v);
          Mat<2,2> jac = DMapQuad2Trig(xi, v);
          double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
          irtrig += IntegrationPoint(x(0), x(1), 0, w);
        }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtrig)
        ip.SetWeight(1.0/6);   // fix consistency for lowest order


    if (GaussRadauIR.Size()==-2)
    { // problem: get negative weights

      // A .. corner, B .. on edge, C .. in cell
      Array<int> classify = { 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2 };
      auto f0 = [](double x, double y) { return 1 ; };
      auto f1 = [](double x, double y)
      {
        return x*x ;
        /*
           double lams[3] = { x, y, 1-x-y };
           if (lams[0] > lams[1] && lams[0] > lams[2])
           {
           return x*y;
           }
           return 0.0;
           */
      };
      auto f2 = [](double x, double y)
      {
        // return x*y*(1-x-y) ;
        double lams[3] = { x, y, 1-x-y };
        if (lams[0] > lams[1] && lams[0] > lams[2])
        {
          return x;
        }
        return 0.0;
      };

      Mat<3,3> mat = 0.0;
      Vec<3> rhs = 0.0;
      // sum_j w_j fi(x_j) = \int fi
      for (auto nr : Range(irtrig))
      {
        auto ip = irtrig[nr];
        mat(0, classify[nr]) += f0(ip(0), ip(1));
        mat(1, classify[nr]) += f1(ip(0), ip(1));
        mat(2, classify[nr]) += f2(ip(0), ip(1));
      }
      cout << "mat = " << endl << mat << endl;
      cout << "Inv(mat) = " << endl << Inv(mat) << endl;
      auto irs = ngcomp::GetIntegrationRules(5);
      IntegrationRule irref = std::move(irs[ET_TRIG]);
      // IntegrationRule irref(ET_TRIG, 5);
      for (auto ip : irref)
      {
        rhs(0) += ip.Weight() * f0(ip(0), ip(1));
        rhs(1) += ip.Weight() * f1(ip(0), ip(1));
        rhs(2) += ip.Weight() * f2(ip(0), ip(1));
      }

      cout << "exact integrals =  " << rhs << endl;
      Vec<3> rhs2 = 0.0;
      for (auto ip : irtrig)
      {
        rhs2(0) += ip.Weight() * f0(ip(0), ip(1));
        rhs2(1) += ip.Weight() * f1(ip(0), ip(1));
        rhs2(2) += ip.Weight() * f2(ip(0), ip(1));
      }
      cout << "GR - integrals =  " << rhs2 << endl;        

      cout.precision(16);
      Vec<3> sol = Inv(mat) * rhs;
      cout << "sol = " << sol << endl;
      cout << "ir = " << irtrig << endl;

      Vec<3> rsol { irtrig[0].Weight(), irtrig[1].Weight(), irtrig[3].Weight() };
      cout << "residual = " << mat*rsol << "-" << rhs << " = " << mat*rsol-rhs << endl;

      for (auto i : Range(irtrig))
        irtrig[i].SetWeight(sol(classify[i]));
    }



    rules[ET_TRIG] = std::move(irtrig);

    IntegrationRule irtet;
    for (int v = 0; v < 4; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          for (auto ipz : GaussRadauIR)
          {
            Vec<3> xi(ipx(0), ipy(0), ipz(0));
            Vec<3> x = MapHex2Tet (xi, v);
            Mat<3,3> jac = DMapHex2Tet(xi, v);
            double w = ipx.Weight()*ipy.Weight()*ipz.Weight()*fabs(Det(jac));
            irtet += IntegrationPoint(x(0), x(1), x(2), w);
          }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtet)
        ip.SetWeight(1.0/24);   // fix consistency for lowest order





    rules[ET_TET] = std::move(irtet);

    IntegrationRule irseg;
    for (auto & ip : GaussRadauIR)
    {
      double x = ip(0), w = ip.Weight();
      irseg.Append ( IntegrationPoint(x/2, 0, 0, 0.5*w));
      irseg.Append ( IntegrationPoint((x+1)/2, 0, 0, 0.5*w));
    }
    rules[ET_SEGM] = std::move(irseg);

    return rules;
  }
#endif
}
