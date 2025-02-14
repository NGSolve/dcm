#ifndef HDIVDUALCELLS_HPP
#define HDIVDUALCELLS_HPP

#include "hdivcells.hpp"


namespace ngcomp
{

  class HDivDualCells : public FESpace
  {
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule TransversalIR;
    IntegrationRule NormalIR;    // Gauss or Gauss-Radau
  public:
    HDivDualCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "hdivdualcells"; }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "HDivDualCells.";
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

#endif
