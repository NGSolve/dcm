#ifndef HDIVCELLS_HPP
#define HDIVCELLS_HPP

#include "dcs_common.hpp"



namespace ngcomp
{
 
  
  template <int DIM>
  class HDivCellFiniteElement : public HDivFiniteElement<DIM>
  {
  public:
    using  HDivFiniteElement<DIM> :: HDivFiniteElement;
    using  HDivFiniteElement<DIM> :: CalcShape;
    
    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    { throw Exception ("CalcAltShape not overloaded"); }

    virtual void CalcPiolaShape (const IntegrationPoint & ip, 
                                 BareSliceMatrix<> shape) const
    { throw Exception ("CalcPiolaShape not overloaded"); }

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
        case PIOLA:
          CalcShape(ip, shape); break;          
        //case COVARIANT:
        //  CalcCovShape(ip, shape); break;
        };
    }

    // du1/dx, du1/dy, du1/dz, du2/dx, ...
    virtual void CalcGradShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> shape) const
    { throw Exception ("CalcGradShape not overloaded"); }
    
  };





  template <int D>
  class DiffOpAltShapeHDiv : public DiffOp<DiffOpAltShapeHDiv<D>>
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcAltShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      Mat<D,D> trafo = 1./Det(F)*F;
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcAltShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      mat = 1./Det(F)*F;
    }

    
  };

  /*

  template <int D>
  class DiffOpCovShape : public DiffOp<DiffOpCovShape<D>>
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcCovShape (mip.IP(), Trans(mat));
      Mat<D,D> Finv = mip.GetJacobianInverse();
      Mat<D,D> trafo = Trans(Finv);
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcCovShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      //Mat<3,3> F = mip.GetJacobian();
      Mat<D,D> Finv = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobianInverse();
      mat = Trans(Finv);
    }

  };
*/

  /*
  template <int D>
  class DiffOpAltDiv : public DiffOp<DiffOpAltDiv<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HDivPrimalCellTet&>(fel).CalcAltDivShape (mip.IP(), Trans(mat));
      Mat<D,D> F = mip.GetJacobian();
      double trafo = 1./Det(F);
      mat *= trafo;
    }

    static int DimRef() { return D; }
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
     
      static_cast<const HDivPrimalCellTet&>(fel).CalcAltDivShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      //Mat<3,3> F = mip.GetJacobian();
      Mat<3,3> F = static_cast<const MappedIntegrationPoint<3,3>&>(mip).GetJacobian();
      mat = 1./Det(F);
    }

    
  };
  */

/*
  template <int D>
  class DiffOpGradHDiv : public DiffOp<DiffOpGradHDiv<D>>
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcGradShape (mip.IP(), Trans(mat));
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
      static_cast<const HDivCellFiniteElement<D>&>(fel).CalcGradShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      throw Exception("DiffOpGradHDiv :: CalcTransformationMatrix needs fixing");
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      auto invFT = Trans(Inv(F));
      for (int i = 0; i < D; i++)
        mat.Rows(i*D,i*D+D).Cols(i+D,i*D+D) = invFT;
    }
  };

 */ 

  
};



#endif


