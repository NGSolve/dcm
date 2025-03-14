#ifndef H1PRIMALCELLS_HPP
#define H1PRIMALCELLS_HPP

#include "h1cells.hpp"

namespace ngcomp
{

  class H1PrimalCellSegm : public H1CellFiniteElement<1>, public VertexOrientedFE<ET_SEGM>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1PrimalCellSegm (int order, const IntegrationRule & _GaussRadauIR);
    using VertexOrientedFE<ET_SEGM>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override;
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const override;
  };



  class H1PrimalCellTrig : public H1CellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1PrimalCellTrig (int order, const IntegrationRule & _GaussRadauIR)
      : H1CellFiniteElement<2> (1+3*order+3*order*order, order),
      GaussRadauIR(_GaussRadauIR)
    {
      VertexOrientedFE<ET_TRIG>::SetVertexNumbers( Array<int>({ 0, 1, 2 }) );
    }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override;
    virtual void CalcL2Shape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override;
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const override;
    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override;
  };

  
  class H1PrimalCells : public FESpace
  {
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
  public:
    H1PrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "h1primalcells"; }

    // documentation
    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "H1PrimalCells.";
      docu.long_docu =
        R"raw_string(Continuous on primal cells.
)raw_string";      
      return docu;
    }

    void Update() override;

    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    shared_ptr<BaseMatrix> GetRotOperator(bool dual = true) const
    {
    if (ma -> GetDimension() == 2)
      return GetRotOperator2DNano(dual);
    else
      throw Exception("not implemented");
    }

    virtual shared_ptr<BaseMatrix> GetMassOperator(
        shared_ptr<CoefficientFunction> rho, 
        shared_ptr<Region> defon, 
        LocalHeap & lh) const override;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(bool fix_lo = true) const;
  protected:
    shared_ptr<BaseMatrix> GetRotOperator2DNano(bool dual) const;
  };

}

#endif
