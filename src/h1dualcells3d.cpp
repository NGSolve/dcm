#include <comp.hpp>
#include "h1dualcells3d.hpp"
#include "dcs_common.hpp"
#include "hdivprimalcells.hpp"
#include "supersparse.hpp"

namespace ngcomp
{
  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;


  class H1DualCellTet : public ScalarFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1DualCellTet (int order, const IntegrationRule & _GaussRadauIR)
      : ScalarFiniteElement<3> (4*(order+1)*(order+1)*(order+1), order),
      GaussRadauIR(_GaussRadauIR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);
      
      shape.Range(GetNDof()) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      ArrayMem<double, 20> polxi(order+1), poleta(order+1), polzeta(order+1);   
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);
      
      auto assign =  [&](int nr, IVec<3> ind)
        {
          shape(nr) = polxi[ind[0]]*poleta[ind[1]]*polzeta[ind[2]];
        };

      
      assign(maxlam, { 0, 0, 0 });
      
      int ii = 4;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < order; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += order;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];
                  // if (vnums[v1] > vnums[v2]) Swap (v1,v2); // optional ? 

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  for (int k = 0, kk=ii; k < order; k++)
                    for (int l = 0; l < order; l++, kk++)
                      {
                        ind[dirv1] = k+1;
                        ind[dirv2] = l+1;
                        assign(kk, ind);
                      }
                }
              ii += sqr(order);
            }
        }
      
      ii += maxlam*order*order*order;
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          for (int k = 0; k < order; k++)
            assign(ii++, { i+1, j+1, k+1 });
    }

    
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      auto dshape = baredshape.AddSize(GetNDof(), 3);
      dshape = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = Trans(Inv(F2*F));

      
      ArrayMem<AutoDiff<1>, 20> polxi(order+1), poleta(order+1), polzeta(order+1);   
      LagrangePolynomials(AutoDiff<1>(xi(0),0), GaussRadauIR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), GaussRadauIR, poleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), GaussRadauIR, polzeta);

      auto assign =  [&](int nr, IVec<3> i)
        {
          dshape.Row(nr) =
          trafo*Vec<3>(polxi[i[0]].DValue(0)*poleta[i[1]].Value()*polzeta[i[2]].Value(),
                       polxi[i[0]].Value()*poleta[i[1]].DValue(0)*polzeta[i[2]].Value(),
                       polxi[i[0]].Value()*poleta[i[1]].Value()*polzeta[i[2]].DValue(0));
        };
      
      
      assign(maxlam, { 0, 0, 0 });
      
      int ii = 4;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < order; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += order;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];
                  // if (vnums[v1] > vnums[v2]) Swap (v1,v2); // optional ? 

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  for (int k = 0, kk=ii; k < order; k++)
                    for (int l = 0; l < order; l++, kk++)
                      {
                        ind[dirv1] = k+1;
                        ind[dirv2] = l+1;
                        assign(kk, ind);
                      }
                }
              ii += sqr(order);
            }
        }
      
      ii += maxlam*order*order*order;
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          for (int k = 0; k < order; k++)
            assign(ii++, { i+1, j+1, k+1 });
    } 
  };
  

  H1DualCells3D::
  H1DualCells3D (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();

    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      GaussRadauIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
  }


  void H1DualCells3D :: Update()
  {
    FESpace::Update();
    size_t ndof = ma->GetNV();
      
    first_edge_dofs.SetSize(ma->GetNEdges()+1);
    first_edge_dofs[0] = ndof;
    for (auto i : Range(ma->GetNEdges()))
      first_edge_dofs[i+1] = ndof += 2*order;

    first_face_dofs.SetSize(ma->GetNFaces()+1);
    first_face_dofs[0] = ndof;
    for (auto i : Range(ma->GetNFaces()))
      first_face_dofs[i+1] = ndof += 3*sqr(order);

    first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
    first_cell_dofs[0] = ndof;
    for (auto i : Range(ma->GetNE(VOL)))
      first_cell_dofs[i+1] = ndof += 4*order*order*order;
      
    SetNDof(ndof);
  }
    
  void H1DualCells3D ::
  GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    if (ei.VB() != VOL) return;
        
    auto el = ma->GetElement(ei);
    for (auto v : el.Vertices())
      dnums.Append(v);
    for (auto ed : el.Edges())
      dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
    for (auto fa : el.Faces())
      dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
    dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
  }


  FiniteElement & H1DualCells3D ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) H1DualCellTet(order, GaussRadauIR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      default:
        throw Exception("element not implemented");
      }
  }

  shared_ptr<BaseMatrix> H1DualCells3D::GetGradientOperator3D (bool dual) const
  {
    //throw Exception(" does not work yet...");
    cout << "called GetGradientOperator" << endl;
    Flags flags;
    flags.SetFlag ("order", order-1);
    flags.SetFlag ("dirichlet", this->flags.GetNumListFlag("dirichlet"));
    flags.SetFlag ("uniform_order", true);
    auto feshdiv = make_shared<HDivPrimalCells>(ma, flags);
    feshdiv -> Update();
    feshdiv -> FinalizeUpdate();
    LocalHeap lh(10*1000*1000);
    shared_ptr<BaseMatrix> sum;
    
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
    
    IntegrationRule irsegm(ET_SEGM, 2*order+10);

    auto altshapes = true;
    
    auto [pnts,ind,num] = GetNanoPoints<ET_TET> (order+1);
    auto [lines,line_ends] = GetNanoLines<ET_TET> (order+1);

    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        
        auto & felh1 = dynamic_cast<const ScalarFiniteElement<3>&> (GetFE (ei, lh));
        auto & felhdiv = dynamic_cast<const HDivPrimalCellTet&> (feshdiv->GetFE (ei, lh));
        auto & trafo = ma->GetTrafo(ei, lh);

        //cout << felhdiv.GetNDof() << endl; 
        Matrix<> gradmat(felhdiv.GetNDof(), felh1.GetNDof());
        Matrix<> mixedmass(felhdiv.GetNDof(), felhdiv.GetNDof());
        Matrix<> lumpedmass(felhdiv.GetNDof(), felhdiv.GetNDof());
        
        gradmat = 0.;
        mixedmass = 0.;
        lumpedmass = 0.;

        Matrix shapehdiv(felhdiv.GetNDof(), 3);        
        Vector shapeh1(felh1.GetNDof());

        //cout << felhdiv.GetNDof() << "?=" << lines.Size() << endl;
        
        bool altshapeshdiv = true;


          // gradmat
          for (auto i : Range(line_ends))
          {
            IntegrationPoint p(pnts[line_ends[i][0]][0], pnts[line_ends[i][0]][1],pnts[line_ends[i][0]][2]);
            felh1.CalcShape(p, shapeh1);    
            gradmat.Row(i) += -shapeh1;
            IntegrationPoint q(pnts[line_ends[i][1]][0], pnts[line_ends[i][1]][1],pnts[line_ends[i][1]][2]);
            felh1.CalcShape(q, shapeh1);   
            gradmat.Row(i) += shapeh1;
          }
        //cout << "gradmat = " << gradmat << endl;
        //throw Exception("test");
      
          // mixedmass
          for (auto i : Range(lines))
          {
            //cout << i << endl;
            for (auto e : lines[i])
              {
                //cout << e[0] << endl;
                Vec<3> micropts[2] =
                  { pnts[e[0]], pnts[e[1]]};
                //cout << micropts[1] << micropts[0]<< endl;
                Vec<3> tau = micropts[1]-micropts[0];
            
                  for (auto ip : irsegm)
                    {
                      Vec<3> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), p(2));
                      if (altshapeshdiv)
                        felhdiv.CalcAltShape(ipseg, shapehdiv);    
                      else
                        felhdiv.CalcShape(ipseg, shapehdiv);    
                    mixedmass.Row(i) += ip.Weight()*(shapehdiv*tau);
                    }
                  }
          }
              
        //cout << "mixedmass = " << mixedmass.Cols(0,10) << endl;
        CalcInverse (mixedmass);
        
        //cout << "mixedmassinv = " << mixedmass.Cols(0,10) << endl;
        //throw Exception("test");
            // integrate
            // lumpedmass ...
            //lumped integration is the same (except order=0 fixup)
            //auto irs = feshdiv->GetIntegrationRules();
            //auto irs = GetIntegrationRules(2*order+6);
            //IntegrationRule ir = std::move(irs[ET_TRIG]);
            
            //auto irs = ngcomp::GetIntegrationRules(2*order+4);
            auto irs = feshdiv->GetIntegrationRules();
            IntegrationRule ir = std::move(irs[ET_TET]);            
            
            Matrix shapes(felhdiv.GetNDof(), ir.Size()*3);
            Matrix Piolashapes(felhdiv.GetNDof(), ir.Size()*3);
        
            for (int i = 0; i < ir.Size(); i++)
              {
                felhdiv.CalcAltShape (ir[i], shapes.Cols(3*i, 3*i+3));
                felhdiv.CalcAltShape (ir[i], Piolashapes.Cols(3*i, 3*i+3));
                //felhdiv.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
                //felhdiv.CalcPiolaShape (ir[i], Piolashapes.Cols(3*i, 3*i+3));
                Piolashapes.Cols(3*i, 3*i+3) *= ir[i].Weight();
              }

            lumpedmass = shapes * Trans (Piolashapes);
        
            //cout << "lumpedmass = " << lumpedmass << endl;

        Matrix<> finalmat(felhdiv.GetNDof(), felh1.GetNDof());
        if (dual)
          finalmat = lumpedmass * mixedmass * gradmat;
        else
          finalmat = mixedmass * gradmat;
        //cout << "finalmat = " << finalmat << endl;
        
        Table<DofId> xdofs(elclass_inds.Size(), felh1.GetNDof()),
          ydofs(elclass_inds.Size(), felhdiv.GetNDof());
        
        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            feshdiv->GetDofNrs(ei, dnumsy);
            GetDofNrs(ei, dnumsx);
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }
        auto mat = make_shared<T_ConstEBE>
          ( feshdiv->GetNDof(),GetNDof(),
           finalmat, std::move(ydofs),std::move(xdofs));
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }
    return sum;
    
  }

  shared_ptr<BaseMatrix> H1DualCells3D ::
  GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    auto irs = GetIntegrationRules();
    IntegrationRule ir = std::move(irs[ET_TET]);

    auto & felref = dynamic_cast<const H1DualCellTet&> (GetFE(ElementId(VOL,0), lh));
    Matrix shapes(felref.GetNDof(), ir.Size());
    Matrix rhoi_shapes(felref.GetNDof(), ir.Size());
    
    for (int i = 0; i < ir.Size(); i++)
      felref.CalcShape (ir[i], shapes.Col(i));

    for (int i = 0; i < shapes.Height(); i++)
      for (int j = 0; j < shapes.Width(); j++)
        if (fabs(shapes(i,j)) < 1e-8)
          shapes(i,j) = 0;
    rhoi_shapes = shapes;
    
    Array<int> rowind, colind;
    Array<double> values;
    Matrix<> elmat(felref.GetNDof());
    Array<DofId> dofs(felref.GetNDof());
    
    for (size_t nr : Range(ma->GetNE()))
      {
        if (defon && !defon->Mask()[ma->GetElIndex(ElementId(VOL,nr))]) continue;
        
        HeapReset hr(lh);
        MappedIntegrationRule<3,3> mir(ir, ma->GetTrafo(ElementId(VOL,nr), lh), lh);
        GetDofNrs(ElementId(VOL,nr), dofs);

        FlatVector<> rhovals(mir.Size(), lh);
        rhovals = 1;
        rho->Evaluate (mir, rhovals.AsMatrix(mir.Size(), 1));
        
        for (size_t i = 0; i < mir.Size(); i++)
          {
            double rhoi = rhovals(i);
              
            //Mat<2,2> F = mir[i].GetJacobian();
            rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
            rhoi_shapes.Col(i) = rhoi * shapes.Col(i);
          }
        elmat = rhoi_shapes * Trans(shapes);
        
        for (int i = 0; i < dofs.Size(); i++)
          for (int j = 0; j < dofs.Size(); j++)
            if (fabs(elmat(i,j)) > 1e-8)
              {
                rowind.Append(dofs[i]);
                colind.Append(dofs[j]);
                values.Append(elmat(i,j));
              }
      }

    auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
    
    return make_shared<MySuperSparseMatrix> (std::move(*spmat));
  }



  std::map<ELEMENT_TYPE, IntegrationRule> H1DualCells3D::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    
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
    
    rules[ET_TET] = std::move(irtet);
    return rules;
  }    
}
