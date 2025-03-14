#ifndef H1CELLS_HPP
#define H1CELLS_HPP
#include "dcs_common.hpp"



namespace ngcomp
{
 
  
  template <int DIM>
  class H1CellFiniteElement : public ScalarFiniteElement<DIM>
  {
  public:
    using  ScalarFiniteElement<DIM> :: ScalarFiniteElement;
    using  ScalarFiniteElement<DIM> :: CalcShape;
    
    virtual void CalcL2Shape (const IntegrationPoint & ip, 
                               BareSliceVector<> shape) const
    { throw Exception ("CalcL2Shape not overloaded"); }

    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceVector<> shape, MAPPINGTYPE type) const
    {
      switch (type)
        {
        case POLYNOMIAL:
          CalcShape(ip, shape); break;
        case L2:
          CalcL2Shape(ip, shape); break;
        };
    }
  };




  template <int D>
  class DiffOpL2ShapeH1 : public DiffOp<DiffOpL2ShapeH1<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      // static_cast<const H1CellFiniteElement<D>&>(fel).CalcPiolaShape (mip.IP(), Trans(mat));
      static_cast<const H1CellFiniteElement<D>&>(fel).CalcL2Shape (mip.IP(), mat.Row(0));
      Mat<D,D> F = mip.GetJacobian();
      double trafo = 1./fabs(Det(F)); // Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<1> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
    }

    static int DimRef() { return D; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      static_cast<const H1CellFiniteElement<D>&>(fel).CalcL2Shape (ip, mat.Row(0));
      // static_cast<const H1CellFiniteElement<D>&>(fel).CalcPiolaShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = 1/fabs(Det(static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian()));
    }
    
  };

};



#endif


