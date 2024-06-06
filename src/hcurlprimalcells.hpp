#include "hcurlcells.hpp"

namespace ngcomp
{
 
  
  class HCurlPrimalCells : public FESpace
  {
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
    IntegrationRule GaussRadauIRinv;
    IntegrationRule tangentialIR;    // Gauss or Gauss-Radau
    bool collocated;   // collocated or staggered grid
  public:
    HCurlPrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "hcurlprimalcells3d"; }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "HCurlPrimalCells3D.";
      docu.long_docu =
        R"raw_string(tangentially continuous on primal cells.
for collocated: polynomial order is order 
for staggered: polynomial order is between order and order+1
)raw_string";
      docu.Arg("collocated") = "bool = True\n"
        "  Use same orders in all directios";
      
      return docu;
    }

    void Update() override;
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual shared_ptr<BaseMatrix> GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };
}
