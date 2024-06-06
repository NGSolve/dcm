#ifndef FILE_SUPERSPARSE
#define FILE_SUPERSPARSE


namespace ngla {

  class MySuperSparseMatrix : public SparseMatrix<double>
  {
    Array<double> blocks1;
    Array<int> blocks1ind;

    Array<Mat<2,2>> blocks2;
    Array<array<int,2>> blocks2ind;

    Array<Mat<3,3>> blocks3;
    Array<array<int,3>> blocks3ind;

    Array<Matrix<double>> gen_blocks;
    Array<Array<int>> gen_ind;

  public:
    MySuperSparseMatrix (SparseMatrix<double> && m);
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;    
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
  };
  

  
  
  class MySuperSparseBlockJacobi : public BaseMatrix
  {
    shared_ptr<const MySuperSparseMatrix> mat;

    Array<double> blocks1;
    Array<int> blocks1ind;

    Array<Mat<2,2>> blocks2;
    Array<array<int,2>> blocks2ind;

    Array<Mat<3,3>> blocks3;
    Array<array<int,3>> blocks3ind;

    Array<Matrix<double>> gen_blocks;
    Array<Array<int>> gen_ind;
    
  public:
    MySuperSparseBlockJacobi (shared_ptr<const MySuperSparseMatrix> amat,
                              shared_ptr<Table<int>> blocks);
    
    bool IsComplex() const override { return false; } 
    int VHeight() const override;
    int VWidth() const override;

    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
  };
}

#endif
