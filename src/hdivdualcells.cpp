#include <comp.hpp>
#include "hdivdualcells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"


namespace ngcomp
{
  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;
  
  class HDivDualCellTrig : public HDivCellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIR;
    const IntegrationRule & tangentialIR;
  public:
    HDivDualCellTrig (const IntegrationRule & _GaussRadauIR,
                       const IntegrationRule & _tangentialIR)
      : HCurlCellFiniteElement<2> (3*2*_GaussRadauIR.Size()*_tangentialIR.Size(), _GaussRadauIR.Size()-1),
      GaussRadauIR(_GaussRadauIR), tangentialIR(_tangentialIR)
    { ; }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const
    {
      ;
    }
    

    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                SliceMatrix<> divshape) const
    { ; } 

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const
    { ; }

    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    { ; }
  };


  class HDivDualCellTet : public HDivCellFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIR;
    const IntegrationRule & tangentialIR;
  public:
    HDivDualCellTet (const IntegrationRule & _GaussRadauIR,
                      const IntegrationRule & _tangentialIR)
      : HDivCellFiniteElement<3> (4*3*sqr(_GaussRadauIR.Size())*_tangentialIR.Size(), _GaussRadauIR.Size()-1),
      GaussRadauIR(_GaussRadauIR), tangentialIR(_tangentialIR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const
    { ; }

    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const
    { ; } 




    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               SliceMatrix<> shape) const
    { ; }




    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    { ; }
  };



  

  template <int D>
  class DiffOpAltShape : public DiffOp<DiffOpAltShape<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcAltShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<D> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }
  };








  

  HDivDualCells::
  HDivDualCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {
    // cout << "HCurlDual ctor, collocated = " << collocated << endl;

    if (ma->GetDimension()==2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivEdge<2>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShape<2>>> ());
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShape<3>>> ());
      }
      
    
    
    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      {
        GaussRadauIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        tangentialIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      }
  }


  void HCurlDualCells :: Update()
  {
    FESpace::Update();
    size_t ndof = 0;

    int ndt = tangentialIR.Size();
    int ndGR = GaussRadauIR.Size();

    first_edge_dofs.SetSize(ma->GetNEdges()+1);
    first_edge_dofs[0] = ndof;
    for (auto i : Range(ma->GetNEdges()))
      first_edge_dofs[i+1] = ndof += 2*ndt;

    first_face_dofs.SetSize(ma->GetNFaces()+1);
    first_face_dofs[0] = ndof;
    for (auto i : Range(ma->GetNFaces()))
      first_face_dofs[i+1] = ndof += 3*2*(ndGR-1)*ndt;

    if (ma->GetDimension() == 3)
      {
        first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
        first_cell_dofs[0] = ndof;
        for (auto i : Range(ma->GetNE(VOL)))
          first_cell_dofs[i+1] = ndof += 4*3*(ndGR-1)*(ndGR-1)*ndt;
      }
    SetNDof(ndof);
  }
    
  void HCurlDualCells ::
  GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    auto el = ma->GetElement(ei);
    for (auto ed : el.Edges())
      dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
    for (auto fa : el.Faces())
      dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
    if (ma->GetDimension()==3 && ei.VB() == VOL)
      dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
    //cout << dnums << endl;
  }


  FiniteElement & HCurlDualCells ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) HCurlDualCellTet(GaussRadauIR, GaussRadauIR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      case ET_TRIG:
        {
          auto trig = new (alloc) HCurlDualCellTrig(GaussRadauIR, GaussRadauIR);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      default:
        throw Exception("element not implemented");
      }
  }



  shared_ptr<BaseMatrix> HCurlDualCells ::
  GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
  {
    if (!collocated)
      throw Exception("only for collocated spaces");

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
          IntegrationRule ir = move(irs[ET_TRIG]);
          auto & felref = dynamic_cast<const HCurlDualCellTrig&> (GetFE(ElementId(VOL,0), lh));
          Matrix shapes(felref.GetNDof(), 2*ir.Size());
          Matrix shapes_trans(2*ir.Size(), felref.GetNDof());
      
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
              auto & felref = dynamic_cast<const HCurlDualCellTrig&> (GetFE(ElementId(VOL,nr), lh));
          
              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Cols(2*i, 2*i+2));
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
                                  Matrix rhoi_shapes_trans(2*ir.Size(), felref.GetNDof());
                                  Matrix<> elmat(felref.GetNDof());
                                  Matrix rhomat(ir.Size(), 4);
                                  Matrix rhoscal(ir.Size(), 1);
              
                                  for (auto i : myrange)
                                    {
                                      HeapReset hr(mylh);
                                      auto nr = elclass_inds[i];

                                      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                                      MappedIntegrationRule<2,2> mir(ir, trafo, mylh);
                                      GetDofNrs(ElementId(VOL,nr), dofs);

                  
                                      if (rho)
                                        {
                                          if (rho->Dimension() == 1)
                                            rho->Evaluate(mir, rhoscal);
                                          else
                                            rho->Evaluate(mir, rhomat);
                                        }
                  
                                      for (size_t i = 0; i < mir.Size(); i++)
                                        {
                                          Mat<2,2> rhoi = Id<2>();
                                          if (rho)
                                            {
                                              if (rho->Dimension() == 1)
                                                rhoi *= rhoscal(i, 0);
                                              else
                                                rhoi = rhomat.Row(i).AsMatrix(2,2);
                                            }
                      
                                          Mat<2,2> Finv = mir[i].GetJacobianInverse();
                                          rhoi = Finv * rhoi * Trans(Finv);
                                          rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                                          rhoi_shapes_trans.Rows(2*i, 2*i+2) = rhoi * shapes_trans.Rows(2*i, 2*i+2);
                                        }
                  
                                      spshapes.Mult (rhoi_shapes_trans, elmat);
                  
                                      for (int i = 0; i < dofs.Size(); i++)
                                        for (int j = 0; j < dofs.Size(); j++)
                                          if (fabs(elmat(i,j)) > 1e-10)
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
      
          return make_shared<MySuperSparseMatrix> (move(*spmat));
        }
      case 3:
        {
          auto irs = this->GetIntegrationRules();
          IntegrationRule ir = move(irs[ET_TET]);
          auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,0), lh));
          Matrix shapes(felref.GetNDof(), 3*ir.Size());
          Matrix shapes_trans(3*ir.Size(), felref.GetNDof());
      
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
              auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,nr), lh));
          
              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
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
                                  Matrix rhoi_shapes_trans(3*ir.Size(), felref.GetNDof());
                                  Matrix<> elmat(felref.GetNDof());
                                  Matrix rhomat(ir.Size(), 9);
                                  Matrix rhoscal(ir.Size(), 1);
              
                                  for (auto i : myrange)
                                    {
                                      HeapReset hr(mylh);
                                      auto nr = elclass_inds[i];

                                      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                                      MappedIntegrationRule<3,3> mir(ir, trafo, mylh);
                                      GetDofNrs(ElementId(VOL,nr), dofs);

                  
                                      if (rho)
                                        {
                                          if (rho->Dimension() == 1)
                                            rho->Evaluate(mir, rhoscal);
                                          else
                                            rho->Evaluate(mir, rhomat);
                                        }
                  
                                      for (size_t i = 0; i < mir.Size(); i++)
                                        {
                                          Mat<3,3> rhoi = Id<3>();
                                          if (rho)
                                            {
                                              if (rho->Dimension() == 1)
                                                rhoi *= rhoscal(i, 0);
                                              else
                                                rhoi = rhomat.Row(i).AsMatrix(3,3);
                                            }
                      
                                          Mat<3,3> Finv = mir[i].GetJacobianInverse();
                                          rhoi = Finv * rhoi * Trans(Finv);
                                          rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                                          rhoi_shapes_trans.Rows(3*i, 3*i+3) = rhoi * shapes_trans.Rows(3*i, 3*i+3);
                                        }
                  
                                      spshapes.Mult (rhoi_shapes_trans, elmat);
                  
                                      for (int i = 0; i < dofs.Size(); i++)
                                        for (int j = 0; j < dofs.Size(); j++)
                                          if (fabs(elmat(i,j)) > 1e-10)
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
      
          return make_shared<MySuperSparseMatrix> (move(*spmat));
        }
      default:
        throw Exception("MassOperator not implemented for dim!=2,3");
      }

    /*
      for (size_t nr : Range(ma->GetNE()))
      {
      HeapReset hr(lh);

      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
      auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,nr), lh));

      for (int i = 0; i < ir.Size(); i++)
      felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
      for (int i = 0; i < shapes.Height(); i++)
      for (int j = 0; j < shapes.Width(); j++)
      if (fabs(shapes(i,j)) < 1e-8)
      shapes(i,j) = 0;
        
        
      MappedIntegrationRule<3,3> mir(ir, trafo, lh);
      GetDofNrs(ElementId(VOL,nr), dofs);
        
      for (size_t i = 0; i < mir.Size(); i++)
      {
      Mat<3,3> rhoi = Id<3>();
      Mat<3,3> Finv = mir[i].GetJacobianInverse();
      rhoi = Finv * rhoi * Trans(Finv);
      rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
      rhoi_shapes.Cols(3*i, 3*i+3) = shapes.Cols(3*i, 3*i+3) * rhoi;
      }
        
      elmat = rhoi_shapes * Trans(shapes);
        
      for (int i = 0; i < dofs.Size(); i++)
      for (int j = 0; j < dofs.Size(); j++)
      if (fabs(elmat(i,j)) > 1e-10)
      {
      rowind.Append(dofs[i]);
      colind.Append(dofs[j]);
      values.Append(elmat(i,j));
      }
      }
    */

  }

  

  
  std::map<ELEMENT_TYPE, IntegrationRule> HCurlDualCells::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    
    IntegrationRule irtet;
    irtet.SetDim(3);
    int i = 0;
    for (int v = 0; v < 4; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          for (auto ipz : GaussRadauIR)
            {
              Vec<3> xi(ipx(0), ipy(0), ipz(0));
              Vec<3> x = MapHex2Tet (xi, v);
              Mat<3,3> jac = DMapHex2Tet(xi, v);
              double w = ipx.Weight()*ipy.Weight()*ipz.Weight()*fabs(Det(jac));
              IntegrationPoint ip(x(0), x(1), x(2), w);
              ip.SetNr(i);
              irtet += ip;
              i++;
            }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtet)
        ip.SetWeight(1.0/24);   // fix consistency for lowest order
    
    rules[ET_TET] = move(irtet);


    IntegrationRule irtrig;
    for (int v = 0; v < 3; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          {
            Vec<2> xi(ipx(0), ipy(0));
            Vec<2> x = MapQuad2Trig (xi, v);
            Mat<2,2> jac = DMapQuad2Trig (xi, v);
            double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
            irtrig += IntegrationPoint(x(0), x(1), 0, w);
          }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtrig)
        ip.SetWeight(1.0/6);   // fix consistency for lowest order

    if (GaussRadauIR.Size()==-2)
      {
        // A .. corner, B .. on edge, C .. in cell
        // Array<int> classify = { 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1 };
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

        auto irs = ngcomp::GetIntegrationRules(5);
        IntegrationRule irref = move(irs[ET_TRIG]);
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
        
        Vec<3> sol = Inv(mat) * rhs;
        cout << "sol = " << sol << endl;
        // cout << "ir = " << irtrig << endl;

        for (auto i : Range(irtrig))
          irtrig[i].SetWeight(sol(classify[i]));
      }



    if (GaussRadauIR.Size()==2)
      {
        // integrate gradients precisely 
        Array<int> classify = { 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2 };

        Array<double> xi, wi;
        ComputeGaussRadauRule (order+1, xi, wi);
        IntegrationRule GaussRadauIRplus;
        for (auto i : Range(xi))
          GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));

        HCurlDualCellPotentialTrig trigpot(GaussRadauIRplus);
        trigpot.SetVertexNumbers (Array<int> { 1, 2, 3 } );
        
        Vector shape(trigpot.GetNDof());
        Matrix dshape(trigpot.GetNDof(), 2);
        Vector sumshape(trigpot.GetNDof());
        Matrix sumdshape(trigpot.GetNDof(), 2);

        sumshape = 0.0;
        sumdshape = 0.0;
        
        // sum_j w_j fi(x_j) = \int fi
        for (auto nr : Range(irtrig))
          {
            auto ip = irtrig[nr];
            // trigpot.CalcShape(ip, shape);
            trigpot.CalcDShape(ip, dshape);
            /*
            sumshape += 
            mat(0, classify[nr]) += f0(ip(0), ip(1));
            mat(1, classify[nr]) += f1(ip(0), ip(1));
            mat(2, classify[nr]) += f2(ip(0), ip(1));
            */
          }

        auto irs = ngcomp::GetIntegrationRules(5);
        IntegrationRule irref = move(irs[ET_TRIG]);
        // IntegrationRule irref(ET_TRIG, 5);
        for (auto ip : irref)
          {
            trigpot.CalcDShape(ip, dshape);
            sumdshape += ip.Weight() * dshape;
          }
        
        // cout << "exact integrals =  " << endl << sumdshape << endl;
        /*
        Vec<3> rhs2 = 0.0;
        for (auto ip : irtrig)
          {
            rhs2(0) += ip.Weight() * f0(ip(0), ip(1));
            rhs2(1) += ip.Weight() * f1(ip(0), ip(1));
            rhs2(2) += ip.Weight() * f2(ip(0), ip(1));
          }
        cout << "GR - integrals =  " << rhs2 << endl;        
        
        Vec<3> sol = Inv(mat) * rhs;
        cout << "sol = " << sol << endl;
        // cout << "ir = " << irtrig << endl;

        for (auto i : Range(irtrig))
          irtrig[i].SetWeight(sol(classify[i]));
        */
      }


    
    rules[ET_TRIG] = move(irtrig);

    

    
    return rules;
  }





  /* ******************** Potential space **************************** */


  

  
  
  
  class HCurlDualCellsPotential3D : public FESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIRplus;
  public:
    HCurlDualCellsPotential3D (shared_ptr<MeshAccess> ama, const Flags & flags)
      : FESpace (ama, flags)            
    {
      if (ma->GetDimension() == 2)
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        }
      else
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        }
      
      Array<double> xi, wi;
      ComputeGaussRadauRule (order, xi, wi);
      for (auto i : Range(xi))
        GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));
      /*
        cout << "HCurlDual - Potential ctor" << endl;
        cout << "order = " << order << endl;
        cout << "IR = " << GaussRadauIRplus << endl;
      */
    }

    string GetClassName () const override { return "hcurldualcellspotential3d"; }

    void Update() override
    {
      FESpace::Update();
      size_t ndof = 0;
      
      ndof += ma->GetNV();

      first_edge_dofs.SetSize(ma->GetNEdges()+1);
      first_edge_dofs[0] = ndof;
      for (auto i : Range(ma->GetNEdges()))
        first_edge_dofs[i+1] = ndof += 2*order-1;

      first_face_dofs.SetSize(ma->GetNFaces()+1);
      first_face_dofs[0] = ndof;
      for (auto i : Range(ma->GetNFaces()))
        first_face_dofs[i+1] = ndof += 3*sqr(order-1) + 3*(order-1); 

      if (ma->GetDimension() == 3)
        {
          first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
          first_cell_dofs[0] = ndof;
          for (auto i : Range(ma->GetNE(VOL)))
            first_cell_dofs[i+1] = ndof += 4*sqr(order-1)*(order-1) + 6*sqr(order-1);
        }
      
      SetNDof(ndof);
    }
    
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override
    {
      dnums.SetSize0();
      auto el = ma->GetElement(ei);
      dnums += el.Vertices();
      for (auto ed : el.Edges())
        dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
      for (auto fa : el.Faces())
        dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
      if (ma->GetDimension()==3)
        if (ei.VB() == VOL)
          dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
    }

    
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override
    {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
        {
        case ET_TET:
          {
            auto tet = new (alloc) HCurlDualCellPotentialTet(GaussRadauIRplus);
            tet->SetVertexNumbers (ngel.vertices);
            return *tet;
          }
        case ET_TRIG:
          {
            auto trig = new (alloc) HCurlDualCellPotentialTrig(GaussRadauIRplus);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }

        default:
          throw Exception("element not implemented");
        }
    }
  };



  shared_ptr<FESpace> HCurlDualCells::GetPotentialSpace() const
  {
    Flags flags;
    flags.SetFlag ("order", order+1);
    flags.SetFlag ("dirichlet", this->flags.GetNumListFlag("dirichlet"));

    auto pot = make_shared<HCurlDualCellsPotential3D>(ma, flags);
    pot -> Update();
    pot -> FinalizeUpdate();
    return pot;
  }

  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator(bool dual) const
  {
    // copy from bilinearform.cpp, line 875 (assemble geometry free)
    // void BilinearForm :: AssembleGF (LocalHeap & lh)

    if (ma -> GetDimension() == 2)
      {
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;
        Flags h1flags;
        h1flags.SetFlag ("order", order);
        auto fescurl = make_shared<H1PrimalCells> (ma, h1flags);
        fescurl->Update();
        fescurl->FinalizeUpdate();
        
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

        IntegrationRule irquad(ET_QUAD, 2*order+4);
        IntegrationRule irsegm(ET_SEGM, 2*order+4);
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(VOL,elclass_inds[0]);
        
            auto & fel = dynamic_cast<const HCurlFiniteElement<2>&> (GetFE (ei, lh));
            auto & felc = dynamic_cast<const ScalarFiniteElement<2>&> (fescurl->GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);

            // cout << "fel.ndof = " << fel.GetNDof() << endl;
            // cout << "felc.ndof = " << felc.GetNDof() << endl;
        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TRIG> (order+1);
            auto faces = GetMicroCellFaces<ET_TRIG> (order+1);
            auto facebnds = GetMicroCellFaceBoundaries<ET_TRIG> (order+1);

            Matrix<> curlmat(faces.Size(), fel.GetNDof());
            Matrix<> mixedmass(faces.Size(), felc.GetNDof());
            Matrix<> lumpedmass(felc.GetNDof(), felc.GetNDof());

            curlmat = 0.;
            mixedmass = 0.;
            lumpedmass = 0.;
            Matrix shape(fel.GetNDof(), 2);
            Vector shapecurl(felc.GetNDof());

            // curlmat ...
            for (auto i : Range(facebnds))
              for (auto fb : facebnds[i])
                {
                  Vec<2> micropts[2] = {
                    pnts[fb[0]],
                    pnts[fb[1]]
                  };
                  Vec<2> tau = micropts[1]-micropts[0];
                  for (auto ip : irsegm)
                    {
                      Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), 0, 0);
                      fel.CalcShape(ipseg, shape);
                      curlmat.Row(i) += ip.Weight() * shape*tau;
                    }
                }
            // cout << "curlmat = " << endl << curlmat << endl;

        
            // mixedmass
            for (auto i : Range(faces))
              for (auto f : faces[i])
                {
                  Vec<2> micropts[4] =
                    { pnts[f[0]], pnts[f[1]], pnts[f[2]], pnts[f[3]] };
              
                  for (auto ip : irquad)
                    {
                      AutoDiff<2> adx(ip(0), 0);
                      AutoDiff<2> ady(ip(1), 1);
                      AutoDiff<2> adlams[] =
                        { (1-adx)*(1-ady), adx*(1-ady), adx*ady, (1-adx)*ady };
                  
                      Vec<2, AutoDiff<2>> adp = AutoDiff<2>(0.0);
                      for (int j = 0; j < 4; j++)
                        adp += adlams[j] * micropts[j];

                      Mat<2,2> jac;                  
                      for (int k = 0; k < 2; k++)
                        for (int l = 0; l < 2; l++)
                          jac(k,l) = adp(k).DValue(l);
                  
                      IntegrationPoint iptrig(adp(0).Value(), adp(1).Value(), 0, 0);
                      felc.CalcShape(iptrig, shapecurl);
                  
                      mixedmass.Row(i) += ip.Weight()*Det(jac) * shapecurl;
                    }
                }

            // cout << "mixedmass = " << mixedmass << endl;
        
            // integrate
            // lumpedmass ...
            auto irs = fescurl->GetIntegrationRules();
            IntegrationRule ir = move(irs[ET_TRIG]);
            Matrix shapes(felc.GetNDof(), ir.Size());
            Matrix rhoi_shapes_trans(ir.Size(), felc.GetNDof());
        
            for (int i = 0; i < ir.Size(); i++)
              felc.CalcShape (ir[i], shapes.Col(i));

            // cout << "shapes = " << endl << shapes << endl;
        
            MappedIntegrationRule<2,2> mir(ir, trafo, lh);
            for (size_t i = 0; i < mir.Size(); i++)
              rhoi_shapes_trans.Row(i) = mir[i].IP().Weight() * shapes.Col(i);
                  
            lumpedmass = shapes * rhoi_shapes_trans;
            // cout << "lumpedmass = " << endl << lumpedmass << endl;
        
            // cout << "mixedmass = " << endl << mixedmass << endl;

            CalcInverse (mixedmass);
            // cout << "mixedmass^-1 = " << mixedmass << endl;
            Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
            tmp = mixedmass * curlmat;

            if (dual)
              curlmat = lumpedmass * tmp;
            else
              curlmat = tmp;              

            // cout << "final mat = " << curlmat << endl;
        
            Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
              ydofs(elclass_inds.Size(), felc.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsx);
                fescurl->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            /*
              cout << "elmat = " << curlmat << endl;
              cout << "xdofs = " << xdofs << endl;
              cout << "ydofs = " << ydofs << endl;
            */
        
            auto mat = make_shared<T_ConstEBE>
              (fescurl->GetNDof(), GetNDof(),
               curlmat, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        return sum;
      }

    else // dim == 3
      
      {
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;
        Flags h1flags;
        h1flags.SetFlag ("order", order);
        auto fescurl = make_shared<HCurlPrimalCells3D> (ma, h1flags);
        fescurl->Update();
        fescurl->FinalizeUpdate();
        
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

        IntegrationRule irquad(ET_QUAD, 2*order+4);
        IntegrationRule irsegm(ET_SEGM, 2*order+4);
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(VOL,elclass_inds[0]);
        
            auto & fel = dynamic_cast<const HCurlFiniteElement<3>&> (GetFE (ei, lh));
            auto & felc = dynamic_cast<const HCurlCellFiniteElement<3>&> (fescurl->GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);

        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TET> (order+1);
            auto faces = GetMicroCellFaces<ET_TET> (order+1);
            auto facebnds = GetMicroCellFaceBoundaries<ET_TET> (order+1);

            Matrix<> curlmat(faces.Size(), fel.GetNDof());
            Matrix<> mixedmass(faces.Size(), felc.GetNDof());
            Matrix<> lumpedmass(felc.GetNDof(), felc.GetNDof());

            cout << "num faces = " << faces.Size() << " =?= " << felc.GetNDof() << " felc.ndof" << endl;
            
            curlmat = 0.;
            mixedmass = 0.;
            lumpedmass = 0.;
            Matrix shape(fel.GetNDof(), 3);
            Matrix shapecurl(felc.GetNDof(), 3);

            // curlmat ...
            for (auto i : Range(facebnds))
              for (auto fb : facebnds[i])
                {
                  Vec<3> micropts[2] = {
                    pnts[fb[0]],
                    pnts[fb[1]]
                  };
                  Vec<3> tau = micropts[1]-micropts[0];
                  for (auto ip : irsegm)
                    {
                      Vec<3> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), p(2), 0);
                      fel.CalcShape(ipseg, shape);
                      curlmat.Row(i) += ip.Weight() * shape*tau;
                    }
                }
            // cout << "curlmat = " << endl << curlmat << endl;

        
            // mixedmass
            for (auto i : Range(faces))
              for (auto f : faces[i])
                {
                  Vec<3> micropts[4] =
                    { pnts[f[0]], pnts[f[1]], pnts[f[2]], pnts[f[3]] };
              
                  for (auto ip : irquad)
                    {
                      AutoDiff<2> adx(ip(0), 0);
                      AutoDiff<2> ady(ip(1), 1);
                      AutoDiff<2> adlams[] =
                        { (1-adx)*(1-ady), adx*(1-ady), adx*ady, (1-adx)*ady };
                  
                      Vec<3, AutoDiff<2>> adp = AutoDiff<2>(0.0);
                      for (int j = 0; j < 4; j++)
                        adp += adlams[j] * micropts[j];

                      Mat<3,2> jac;                  
                      for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 2; l++)
                          jac(k,l) = adp(k).DValue(l);
                  
                      IntegrationPoint iptrig(adp(0).Value(), adp(1).Value(), adp(2).Value(), 0);
                      // felc.CalcPiolaShape(iptrig, shapecurl);
                      felc.CalcAltShape(iptrig, shapecurl);
                  
                      // mixedmass.Row(i) += ip.Weight()*L2Norm(Cross(jac.Col(0), jac.Col(1))) * shapecurl;
                      mixedmass.Row(i) += ip.Weight()* shapecurl * Cross(jac.Col(0), jac.Col(1));

                      /*
                      *testout << "mixed, i = " << i << endl;
                      *testout << "jac = " << jac << ", n = " << Cross(jac.Col(0), jac.Col(1)) << endl;
                      *testout << "shape = " << shapecurl << endl;
                      *testout << "shapen = " << shapecurl*Cross(jac.Col(0), jac.Col(1)) << endl;                      
                      */
                    }
                }

            // cout << "mixedmass = " << mixedmass << endl;
        
            // integrate
            // lumpedmass ...
            /*
            auto irs = fescurl->GetIntegrationRules();
            IntegrationRule ir = move(irs[ET_TET]);
            */
            auto irs = ngcomp::GetIntegrationRules(2*order+4);
            IntegrationRule ir = move(irs[ET_TET]);            
            
            Matrix shapes(felc.GetNDof(), ir.Size()*3);
            Matrix Piolashapes(felc.GetNDof(), ir.Size()*3);
            Matrix rhoi_shapes_trans(ir.Size()*3, felc.GetNDof());
        
            for (int i = 0; i < ir.Size(); i++)
              {
                felc.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
                // felc.CalcPiolaShape (ir[i], Piolashapes.Cols(3*i, 3*i+3));
                felc.CalcAltShape (ir[i], Piolashapes.Cols(3*i, 3*i+3));
                Piolashapes.Cols(3*i, 3*i+3) *= ir[i].Weight();
              }

            /*
            MappedIntegrationRule<3,3> mir(ir, trafo, lh);
            for (size_t i = 0; i < mir.Size(); i++)
              {
                Mat<3,3> rhoi = Id<3>();
                // Mat<3,3> Finv = mir[i].GetJacobianInverse();
                // rhoi = Finv * rhoi * Trans(Finv);
                rhoi *= ir[i].Weight(); //  * mir[i].GetJacobiDet();
                rhoi_shapes_trans.Rows(3*i,3*i+3) = rhoi * Trans(shapes.Cols(3*i,3*i+3));
              } 
            lumpedmass = shapes * rhoi_shapes_trans;
            */

            lumpedmass = shapes * Trans (Piolashapes);
            
            *testout << "curlmat = " << curlmat << endl;
            *testout << "mixedmass = " << mixedmass << endl;
            *testout << "lumpedmass = " << lumpedmass << endl;
            
            CalcInverse (mixedmass);
            *testout << "inv mixedmass = " << mixedmass << endl;            
            Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
            tmp = mixedmass * curlmat;

            if (dual)
              curlmat = lumpedmass * tmp;
            else
              curlmat = tmp;              
            

            *testout << "final mat = " << curlmat << endl;
        
            Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
              ydofs(elclass_inds.Size(), felc.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsx);
                fescurl->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            /*
              cout << "elmat = " << curlmat << endl;
              cout << "xdofs = " << xdofs << endl;
              cout << "ydofs = " << ydofs << endl;
            */
        
            auto mat = make_shared<T_ConstEBE>
              (fescurl->GetNDof(), GetNDof(),
               curlmat, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        return sum;
      }

  }
}
