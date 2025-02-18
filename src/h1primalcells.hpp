namespace ngcomp
{


  enum MAPPINGTYPEH1 { NONE, L2, L2POLYNOMIAL };

  class H1PrimalCellSegm : public ScalarFiniteElement<1>, public VertexOrientedFE<ET_SEGM>
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



  class H1PrimalCellTrig : public ScalarFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1PrimalCellTrig (int order, const IntegrationRule & _GaussRadauIR)
      : ScalarFiniteElement<2> (1+3*order+3*order*order, order),
      GaussRadauIR(_GaussRadauIR)
    {
      VertexOrientedFE<ET_TRIG>::SetVertexNumbers( Array<int>({ 0, 1, 2 }) );
    }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override;
    virtual void CalcL2AltShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const;
    virtual void CalcL2Shape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const;
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape, MAPPINGTYPEH1 type) const 
    {
      switch (type)
        {
        case NONE:
          CalcShape(ip, shape); break;
        case L2:
          CalcL2Shape(ip, shape); break;          
        case L2POLYNOMIAL:
          CalcL2AltShape(ip, shape); break;          
        };
    }
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

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  protected:
    shared_ptr<BaseMatrix> GetRotOperator2DNano(bool dual) const;
  };



  class DiffOpL2Shape2D : public DiffOp<DiffOpL2Shape2D>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const H1PrimalCellTrig&>(fel).CalcL2Shape (mip.IP(), mat.Row(0));
      Mat<2,2> F = mip.GetJacobian();
      double trafo = 1./fabs(Det(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<1> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }

    static int DimRef() { return 2; }
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const H1PrimalCellTrig&>(fel).CalcL2Shape (ip, mat.Row(0));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = 1/fabs(Det(static_cast<const MappedIntegrationPoint<2,2>&>(mip).GetJacobian()));
    }

    
  };

}
