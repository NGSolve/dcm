#include <comp.hpp>
#include "hdivdualcells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"
#include "intrules.hpp"


namespace ngcomp
{
  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;
  
  class HDivDualCellTrig : public HDivCellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & TransversalIR; //transversal to vector
    const IntegrationRule & NormalIR;  // in direction of vector
  public:
    HDivDualCellTrig (const IntegrationRule & _GaussRadauIR,
                       const IntegrationRule & _tangentialIR)
      : HDivFiniteElement<2> (3*2*_TransversalIR.Size()*_NormalIR.Size(), _NormalIR.Size()-1),      
      TransversalIR(_TransversalRadauIR), NormalIR(_NormalIR)
    { ; }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    { ; }
    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> divshape) const override
    { ; } 

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
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
      : HDivCellFiniteElement<3> (4*3*_TransversalIR.Size()*sqr(_NormalIR.Size()), _NormalIR.Size()-1),
      TransversalIR(_TransversalIR), NormalIR(_normalIR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    { ; }

    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> curlshape) const override
    { ; } 



    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const override
    { ; }




    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    { ; }
  };




  HDivDualCells::
  HDivDualCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {

    if (ma->GetDimension()==2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivEdge<2>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<2>>> ());
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<3>>> ());
      }
      
    
    
    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      {
        TransversalIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        NormalIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      }
  }


  void HDivDualCells :: Update()
  {
    FESpace::Update();
    size_t ndof = 0;

    int ndT = TransversalIR.Size();
    int ndN = NormalIR.Size();

    if (ma->GetDimension() == 2)
      {
        first_edge_dofs.SetSize(ma->GetNEdges()+1);
        first_edge_dofs[0] = ndof;
        for (auto i : Range(ma->GetNEdges()))
          first_edge_dofs[i+1] = ndof += 2*ndN;


        first_face_dofs.SetSize(ma->GetNFaces()+1);
        first_face_dofs[0] = ndof;
        for (auto i : Range(ma->GetNFaces()))
          first_face_dofs[i+1] = ndof += 3*2*(ndT-1)*ndN;
      }
    if (ma->GetDimension() == 3)
      {
        first_face_dofs.SetSize(ma->GetNFaces()+1);
        first_face_dofs[0] = ndof;
        for (auto i : Range(ma->GetNFaces()))
          first_face_dofs[i+1] = ndof += 3*ndN*ndN;

        first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
        first_cell_dofs[0] = ndof;
        for (auto i : Range(ma->GetNE(VOL)))
          first_cell_dofs[i+1] = ndof += 4*3*ndN*ndN*(ndT-1);
      }
    SetNDof(ndof);
  }
    
  void HDivDualCells ::
  GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    auto el = ma->GetElement(ei);
    for (auto fa : el.Faces())
      dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
    if (ma->GetDimension()==3 && ei.VB() == VOL)
      dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
  }


  FiniteElement & HCurlDualCells ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) HDivDualCellTet(GaussRadauIR, GaussRadauIR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      case ET_TRIG:
        {
          auto trig = new (alloc) HDivDualCellTrig(GaussRadauIR, GaussRadauIR);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      default:
        throw Exception("element not implemented");
      }
  }



  shared_ptr<BaseMatrix> HDivDualCells ::
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
          auto & felref = dynamic_cast<const HDivDualCellTrig&> (GetFE(ElementId(VOL,0), lh));
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
              auto & felref = dynamic_cast<const HDivDualCellTrig&> (GetFE(ElementId(VOL,nr), lh));
          
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
                      
                                          Mat<2,2> F = mir[i].GetJacobian();
                                          rhoi = Trans(F) * rhoi * F;
                                          rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
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
          auto & felref = dynamic_cast<const HDivDualCellTet&> (GetFE(ElementId(VOL,0), lh));
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
              auto & felref = dynamic_cast<const HDivDualCellTet&> (GetFE(ElementId(VOL,nr), lh));
          
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
                      
                                          Mat<3,3> F = mir[i].GetJacobian();
                                          rhoi = TransF * rhoi * F;
                                          rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
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

  }

  

  
  std::map<ELEMENT_TYPE, IntegrationRule> HDivDualCells::GetIntegrationRules() const
  {
      std::map<ELEMENT_TYPE, IntegrationRule> rules;
      rules[ET_TRIG] = PrimalCellIR(NormalIR,true);

      rules[ET_TET] = PrimalVolIR(NormalIR,true);

      auto irseg = SelectIntegrationRule(ET_SEGM, 2*order+2);
      rules[ET_SEGM] = PrimalSegmIR(irseg);



      return rules;
  }

};
