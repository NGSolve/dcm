#include "hdivcells.hpp"

/*
#include "smallsparse.hpp"
#include "intrules.hpp"
#include "supersparse.hpp"
*/



namespace ngcomp
{
  // Necessary in header file for some experimental H1 gradients?
  class HDivPrimalCellTrig : public HDivCellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & IR;
  public:
    HDivPrimalCellTrig (const IntegrationRule & _IR)
      : HDivCellFiniteElement<2> (3*_IR.Size()*(2*_IR.Size()-1), _IR.Size()-1),      
      IR(_IR)
    { ; }
    
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcAltCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               SliceVector<> divshape) const;

    Mat<2,2> GetPiola(const IntegrationPoint & ip) const;
    double GetDet(const IntegrationPoint & ip) const;
  };

  class HDivPrimalCellTet : public HDivCellFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & IR;
  public:
    HDivPrimalCellTet (const IntegrationRule & _IR)
      : HDivCellFiniteElement<3> (6*_IR.Size()*_IR.Size()*(2*_IR.Size()-1), _IR.Size()-1),      
      IR(_IR)
    { ; }
    
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const;

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                               SliceVector<> divshape) const;

    virtual void CalcAltDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> divshape) const;
    //Mat<3,3> GetPiola(const IntegrationPoint & ip) const;
    //double GetDet(const IntegrationPoint & ip) const;
  };


  class HDivPrimalCells : public FESpace
  {
    Array<DofId> first_element_dofs;
    IntegrationRule IR;
    //bool uniform_order;
    
  public:

    HDivPrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags);

    string GetClassName () const override { return "HDivPrimalCells"; }

    // documentation
    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.short_docu = "HDiv conforming space on primal cells which is non-smooth but normal continuous on dual edges.";
      docu.long_docu =
        R"raw_string(HDiv conforming space on primal cells which is non-smooth but normal continuous on dual edges.
)raw_string";      
      return docu;
    }

    void Update() override;
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override
    {
      dnums.SetSize0();
      if (ei.VB() != VOL) return;
      dnums += Range(first_element_dofs[ei.Nr()], first_element_dofs[ei.Nr()+1]);
    }
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    
    virtual shared_ptr<BaseMatrix> GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const override;
    //shared_ptr<BaseMatrix> ConvertGR2GOperator() const;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };
}    
