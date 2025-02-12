#ifndef HCURLDUALCELLS_HPP
#define HCURLDUALCELLS_HPP


#include "hcurlcells.hpp"

namespace ngcomp
{


  class HCurlDualCells : public FESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIR;
    IntegrationRule tangentialIR;    // Gauss or Gauss-Radau
    bool collocated;   // collocated or staggered grid
  public:
    HCurlDualCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "hcurldualcells3d"; }

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "HCurlDualCells3D.";
      docu.long_docu =
        R"raw_string(tangentially continuous on dual cells.
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
    
    shared_ptr<FESpace> GetPotentialSpace(bool include_central, bool distributional=false) const;
    shared_ptr<BaseMatrix> GetGradientOperator(bool nanocell = false) const;
    shared_ptr<BaseMatrix> GetCurlOperator(bool dual = true, bool altshapes  = true, bool lumping = true, bool Kronecker=false) const;
    shared_ptr<BaseMatrix> GetRotOperator(bool altshapes = true) const;
    shared_ptr<BaseMatrix> GetNanoOperator(bool altshapes = true) const;
    shared_ptr<BaseMatrix> GetNanoOperatorRect(shared_ptr<FESpace> feshcurl, bool altshapes = true) const;
    shared_ptr<BaseMatrix> GetMassOperator2DNano(bool altshapes = true) const;

  protected:
    shared_ptr<BaseMatrix> GetCurlOperator3DNano(bool dual, bool altshapes, bool lumping) const;
    shared_ptr<BaseMatrix> GetCurlOperator3D(bool dual) const;
    shared_ptr<BaseMatrix> GetCurlOperator2DNano(bool dual, bool altshapes, bool lumping) const;
    shared_ptr<BaseMatrix> GetCurlOperator2D(bool dual) const;
    shared_ptr<BaseMatrix> GetCurlOperator2DKronecker(bool dual) const;
    shared_ptr<BaseMatrix> GetGradientOperator3D() const;
    
    shared_ptr<BaseMatrix> GetRotOperator2DNano(bool altshapes) const;
  };




  
  class HCurlDualCellsPotential3D : public PotentialFESpace
  {
    Array<DofId> first_edge_dofs;
    Array<DofId> first_face_dofs;
    Array<DofId> first_cell_dofs;
    IntegrationRule GaussRadauIRplus;
    bool include_central = false;
    bool distributional = false;
  public:
    HCurlDualCellsPotential3D (shared_ptr<MeshAccess> ama, const Flags & flags);
    string GetClassName () const override;
    void Update() override;
    void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    BitArray GetInnerDofs() const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
  };
  
}

#endif
