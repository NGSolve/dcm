#ifndef H1DUALCELLS_HPP
#define H1DUALCELLS_HPP

#include "h1cells.hpp"

namespace ngcomp
{

  class H1DualCells : public FESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
  public:

    H1DualCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "h1dualcells"; }

    // documentation
    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "H1 on Dual Cells";
      docu.long_docu =
        R"raw_string(Continuous on dual cells
)raw_string";      
      return docu;
    }

    // organzize the FESpace, called after every mesh update
    void Update() override;
    
    // dof-numbers for element-id ei
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    
    // generate FiniteElement for element-id ei
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual shared_ptr<BaseMatrix>
    GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;

    virtual shared_ptr<BaseMatrix> GetGradientOperator2D(bool dual = true) const;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(bool fix_lo = true) const;
  };
}


#endif
