namespace ngcomp
{
  
  class L2CellTrig : public ScalarFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussIR;
  public:
    L2CellTrig (int order, const IntegrationRule & _GaussIR)
      : ScalarFiniteElement<2> (3*sqr(order+1), order),
        GaussIR(_GaussIR)
    {
      VertexOrientedFE<ET_TRIG>::SetVertexNumbers( Array<int>({ 0, 1, 2 }) );
    }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
                                
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override;

    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const override
    {
      throw Exception("L2CellTrig::CalcDShape not overloaded");
    }

    
    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override;
  };

  class L2Cells : public FESpace
  {
    Array<DofId> first_element_dofs;
    IntegrationRule GaussIR;
  public:

    L2Cells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "l2cells"; }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "L2 on cells";
      docu.long_docu =
        R"raw_string(1/J*Qk
)raw_string";      
      return docu;
    }

    void Update() override;
    
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    // virtual shared_ptr<BaseMatrix>
    // GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;
  };
}
