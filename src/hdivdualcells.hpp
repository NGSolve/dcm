namespace ngcomp
{

  class HDivDualCells : public FESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
    IntegrationRule tangentialIR;    // Gauss or Gauss-Radau
    bool collocated;   // collocated or staggered grid
  public:
    HDivDualCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "hcurldualcells3d"; }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "HDivDualCells3D.";
      docu.long_docu =
        R"raw_string(normal continuous on dual cells.
)raw_string";
      
      return docu;
    }

    void Update() override;
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual shared_ptr<BaseMatrix> GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
    
    //shared_ptr<FESpace> GetPotentialSpace() const;
    //shared_ptr<BaseMatrix> GetDivOperator() const;
    
  };
}
