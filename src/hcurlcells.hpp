#ifndef HCURLCELLS_HPP
#define HCURLCELLS_HPP




namespace ngcomp
{
 
  enum MAPPINGTYPE { POLYNOMIAL, COVARIANT, PIOLA, PIOLAPOLYNOMIAL };
  
  template <int DIM>
  class HCurlCellFiniteElement : public HCurlFiniteElement<DIM>
  {
  public:
    using  HCurlFiniteElement<DIM> :: HCurlFiniteElement;
    using  HCurlFiniteElement<DIM> :: CalcShape;
    
    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    { throw Exception ("CalcAltShape not overloaded"); }
    virtual void CalcAltCurlShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    { throw Exception ("CalcAltShape not overloaded"); }
    virtual void CalcPiolaShape (const IntegrationPoint & ip, 
                                 BareSliceMatrix<> shape) const
    { throw Exception ("CalcPiolaShape not overloaded"); }

    virtual void CalcPiolaAltShape (const IntegrationPoint & ip, 
                                 BareSliceMatrix<> shape) const
    { throw Exception ("CalcPiolaAltShape not overloaded"); }

    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceMatrix<> shape, bool alt) const
    {
      if (alt)
        CalcAltShape(ip, shape);
      else
        CalcShape(ip, shape);
    }

    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceMatrix<> shape, MAPPINGTYPE type) const
    {
      switch (type)
        {
        case POLYNOMIAL:
          CalcAltShape(ip, shape); break;
        case COVARIANT:
          CalcShape(ip, shape); break;
        case PIOLA:
          CalcPiolaShape(ip, shape); break;          
        case PIOLAPOLYNOMIAL:
          CalcPiolaAltShape(ip, shape); break;          
        };
    }

    // du1/dx, du1/dy, du1/dz, du2/dx, ...
    virtual void CalcGradShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> shape) const
    { throw Exception ("CalcGradShape not overloaded"); }
  };





  template <int D>
  class DiffOpAltShape : public DiffOp<DiffOpAltShape<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<D> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }

    static int DimRef() { return D; }
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = Trans(static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobianInverse());
    }

    
  };

  template <int D>
  class DiffOpAltCurl : public DiffOp<DiffOpAltCurl<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltCurlShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = 1/Det(F)*F;
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<D> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }

    static int DimRef() { return D; }
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltCurlShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      mat = 1/Det(F)*F;
    }

    
  };

  template <int D>
  class DiffOpPiolaShape : public DiffOp<DiffOpPiolaShape<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcPiolaShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<D> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }


    static int DimRef() { return D; }     

    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcPiolaShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = Trans(static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobianInverse());
    }

    
  };



  template <int D>
  class DiffOpHodge : public DiffOp<DiffOpHodge<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      // static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcPiolaShape (mip.IP(), Trans(mat));
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltShape (mip.IP(), Trans(mat));
      // static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcShape (mip.IP(), Trans(mat));      
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = 1/Det(F) * F; // Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<D> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }

    static int DimRef() { return D; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcAltShape (ip, Trans(mat));
      // static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcPiolaShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      mat = 1/Det(F) * F;   
    }
    
  };




  template <int D>
  class DiffOpGradHCurl : public DiffOp<DiffOpGradHCurl<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 1 };

    static bool SupportsVB (VorB checkvb) { return true; }

    static IVec<2> GetDimensions() { return { D,D }; }

    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcGradShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          auto matshape = mat.Col(i).AsMatrix(D,D);
          Mat<D,D> transshape = trafo * matshape * Trans(trafo);
          mat.Col(i).AsMatrix(D,D) = transshape;
        }
    }

    static int DimRef() { return D*D; }
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const HCurlCellFiniteElement<D>&>(fel).CalcGradShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      throw Exception("DiffOpGradHCurl :: CalcTransformationMatrix needs fixing");
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      auto invFT = Trans(Inv(F));
      for (int i = 0; i < D; i++)
        mat.Rows(i*D,i*D+D).Cols(i+D,i*D+D) = invFT;
    }
  };

  
    


  

  class PotentialFESpace : public FESpace
  {
  public:
    using FESpace::FESpace;
    virtual BitArray GetInnerDofs() const = 0;
  };

  
};



#endif


