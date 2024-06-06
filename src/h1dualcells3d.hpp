namespace ngcomp
{

  class H1DualCells3D : public FESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
  public:
    H1DualCells3D (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "h1dualcells"; }

    // documentation
    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "H1DualCells3D.";
      docu.long_docu =
        R"raw_string(Continuous on dual cells.
)raw_string";      
      return docu;
    }

    void Update() override;
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;


    //std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(optional<int> intorder) const;

    virtual shared_ptr<BaseMatrix> GetGradientOperator3D(bool dual = true) const;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;

    virtual shared_ptr<BaseMatrix>
    GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;

  };
}
