#include <comp.hpp>
#include "smallsparse.hpp"
#include "intrules.hpp"
#include "supersparse.hpp"


namespace ngla {


  MySuperSparseMatrix :: MySuperSparseMatrix (SparseMatrix<double> && m)
    : SparseMatrix<double> (m)
  {

    size_t n = Height();

    // find dofs in the same block
    Array<int> cluster(n);
    for (auto i : Range(n)) cluster[i] = i;
    bool done = false;
    while (!done)
      {
        done = true;
        for (auto i : Range(n))
          for (auto j : GetRowIndices(i))
            {
              if (cluster[i] != cluster[j])
                {
                  done = false;
                  auto minv = min(cluster[i], cluster[j]);
                  cluster[i] = minv;
                  cluster[j] = minv;
                }
            }
      }

    TableCreator<int> creator(n);
    for ( ; !creator.Done(); creator++)
      {
        for (auto i : Range(n))
          creator.Add (cluster[i], i);
      }

    auto blocks = make_shared<Table<int>>(creator.MoveTable());

    auto & constmat = (const SparseMatrix<double>&) *this;
    
    for (auto i : Range(*blocks))
      {
        auto blocki = (*blocks)[i];

        switch (blocki.Size())
          {
          case 0:
            break;

          case 1:
            {
              blocks1.Append (constmat(blocki[0], blocki[0]));
              blocks1ind.Append (blocki[0]);
              break;
            }

          case 2:
            {
              constexpr int BS = 2;
              array<int,BS> blockind { blocki[0], blocki[1] };
              Mat<BS,BS> block;
              for (int j = 0; j < BS; j++)
                for (int k = 0; k < BS; k++)
                  block(j,k) = constmat(blocki[j], blocki[k]);
              blocks2.Append (block);
              blocks2ind.Append (blockind);
              break;
            }

          case 3:
            {
              constexpr int BS = 3;
              array<int,BS> blockind { blocki[0], blocki[1], blocki[2] };
              Mat<BS,BS> block;
              for (int j = 0; j < BS; j++)
                for (int k = 0; k < BS; k++)
                  block(j,k) = constmat(blocki[j], blocki[k]);
              blocks3.Append (block);
              blocks3ind.Append (blockind);
              break;
            }

          default:
            {
              int bs = blocki.Size();
              // cout << "bs = " << bs << endl;
              Matrix block(bs);
              for (int j = 0; j < bs; j++)
                for (int k = 0; k < bs; k++)
                  block(j,k) = constmat(blocki[j], blocki[k]);

              gen_blocks.Append (block);
              gen_ind.Append (Array<int> (blocki));
              // throw Exception ("block size "+ToString(blocki.Size())+" not implemented");
            }
          }
      }
  }


  void MySuperSparseMatrix ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    // SparseMatrix<double>::MultAdd (s, x, y);
    // return;
    
    static Timer t("supersparse::MultAdd"); RegionTimer reg(t);
    /*
    static Timer t1("supersparse::MultAdd1");
    static Timer t2("supersparse::MultAdd2");    
    static Timer t3("supersparse::MultAdd3");
    */
    
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    FlatArray<double> b1(blocks1);
    FlatArray<int> i1(blocks1ind);
    ParallelFor (Range(blocks1), [b1,i1,fx,fy,s] (size_t i)
                 {
                   fy(i1[i]) += s*b1[i] * fx(i1[i]);
                 }, TasksPerThread(4));
    FlatArray<Mat<2,2>> b2(blocks2);
    FlatArray<array<int,2>> i2(blocks2ind);
    // t2.Start();
    ParallelFor (Range(blocks2), [b2,i2,fx,fy,s] (size_t i)
                 {
                   Vec<2> vx(fx(i2[i][0]), fx(i2[i][1]));
                   Vec<2> vy = s * b2[i] * vx;
                   fy(i2[i][0]) += vy(0);
                   fy(i2[i][1]) += vy(1);
                 }, TasksPerThread(4));
    // t2.Stop();
    FlatArray<Mat<3,3>> b3(blocks3);
    FlatArray<array<int,3>> i3(blocks3ind);
    // t3.Start();
    ParallelFor (Range(blocks3), [b3,i3,fx,fy,s] (size_t i)
                 {
                   Vec<3> vx(fx(i3[i][0]), fx(i3[i][1]), fx(i3[i][2]));
                   Vec<3> vy = s * b3[i] * vx;
                   fy(i3[i][0]) += vy(0);
                   fy(i3[i][1]) += vy(1);
                   fy(i3[i][2]) += vy(2);
                 }, TasksPerThread(4));
    // t3.Stop();

    ParallelFor (Range(gen_blocks), [&] (size_t i)
                 {
                   FlatArray<int> ind = gen_ind[i];
                   VectorMem<20> vx(ind.Size());
                   VectorMem<20> vy(ind.Size());

                   vx = fx(ind);
                   vy = s * gen_blocks[i] * vx;
                   fy(ind) += vy;
                 }, TasksPerThread(4));

    
  }








  

  
  shared_ptr<BaseMatrix> MySuperSparseMatrix ::
  InverseMatrix (shared_ptr<BitArray> subset) const
  {
    static Timer t("supersparse::inverse"); RegionTimer reg(t);
    cout << IM(3) << "MySuperSpasrseMatrix::Inverse called" << endl;
    size_t n = Height();

    // find dofs in the same block
    Array<int> cluster(n);
    for (auto i : Range(n)) cluster[i] = i;
    bool done = false;
    while (!done)
      {
        done = true;
        for (auto i : Range(n))
          for (auto j : GetRowIndices(i))
            {
              if (cluster[i] != cluster[j])
                {
                  done = false;
                  auto minv = min(cluster[i], cluster[j]);
                  cluster[i] = minv;
                  cluster[j] = minv;
                }
            }
      }

    TableCreator<int> creator(n);
    for ( ; !creator.Done(); creator++)
      {
        if (!subset)
          for (auto i : Range(n))
            creator.Add (cluster[i], i);
        else
          for (auto i : Range(n))
            if (subset->Test(i))
              creator.Add (cluster[i], i);
      }

    auto blocks = make_shared<Table<int>>(creator.MoveTable());

    return make_shared<MySuperSparseBlockJacobi>
      (dynamic_pointer_cast<const MySuperSparseMatrix>
       (this->shared_from_this()), blocks);
      
    // return CreateBlockJacobiPrecond(blocks);
  }



  
  MySuperSparseBlockJacobi ::
  MySuperSparseBlockJacobi (shared_ptr<const MySuperSparseMatrix> amat,
                            shared_ptr<Table<int>> blocks)
    : mat(amat)
  {
    for (auto i : Range(*blocks))
      {
        auto blocki = (*blocks)[i];

        switch (blocki.Size())
          {
          case 0:
            break;

          case 1:
            {
              blocks1.Append (1.0/(*mat)(blocki[0], blocki[0]));
              blocks1ind.Append (blocki[0]);
              break;
            }

          case 2:
            {
              constexpr int BS = 2;
              array<int,BS> blockind { blocki[0], blocki[1] };
              Mat<BS,BS> block;
              for (int j = 0; j < BS; j++)
                for (int k = 0; k < BS; k++)
                  block(j,k) = (*mat)(blocki[j], blocki[k]);
              blocks2.Append (Inv(block));
              blocks2ind.Append (blockind);
              break;
            }

          case 3:
            {
              constexpr int BS = 3;
              array<int,BS> blockind { blocki[0], blocki[1], blocki[2] };
              Mat<BS,BS> block;
              for (int j = 0; j < BS; j++)
                for (int k = 0; k < BS; k++)
                  block(j,k) = (*mat)(blocki[j], blocki[k]);
              blocks3.Append (Inv(block));
              blocks3ind.Append (blockind);
              break;
            }

          default:
            {
              int bs = blocki.Size();
              // cout << "bs = " << bs << endl;
              Matrix block(bs);
              for (int j = 0; j < bs; j++)
                for (int k = 0; k < bs; k++)
                  block(j,k) = (*mat)(blocki[j], blocki[k]);
              CalcInverse (block);

              gen_blocks.Append (block);
              gen_ind.Append (Array<int> (blocki));
              // throw Exception ("block size "+ToString(blocki.Size())+" not implemented");
            }
          }
      }
  }

  void MySuperSparseBlockJacobi ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    static Timer t("supersparseBJ::MultAdd"); RegionTimer reg(t);
    static Timer t1("supersparseBJ::MultAdd1");
    static Timer t2("supersparseBJ::MultAdd2");    
    static Timer t3("supersparseBJ::MultAdd3");
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    FlatArray<double> b1(blocks1);
    FlatArray<int> i1(blocks1ind);
    t1.Start();
    ParallelFor (Range(blocks1), [b1,i1,fx,fy,s] (size_t i)
                 {
                   fy(i1[i]) += s*b1[i] * fx(i1[i]);
                 }, TasksPerThread(4));
    t1.Stop();
    FlatArray<Mat<2,2>> b2(blocks2);
    FlatArray<array<int,2>> i2(blocks2ind);
    t2.Start();
    ParallelFor (Range(blocks2), [b2,i2,fx,fy,s] (size_t i)
                 {
                   Vec<2> vx(fx(i2[i][0]), fx(i2[i][1]));
                   Vec<2> vy = s * b2[i] * vx;
                   fy(i2[i][0]) += vy(0);
                   fy(i2[i][1]) += vy(1);
                 }, TasksPerThread(4));
    t2.Stop();
    FlatArray<Mat<3,3>> b3(blocks3);
    FlatArray<array<int,3>> i3(blocks3ind);
    t3.Start();
    ParallelFor (Range(blocks3), [b3,i3,fx,fy,s] (size_t i)
                 {
                   Vec<3> vx(fx(i3[i][0]), fx(i3[i][1]), fx(i3[i][2]));
                   Vec<3> vy = s * b3[i] * vx;
                   fy(i3[i][0]) += vy(0);
                   fy(i3[i][1]) += vy(1);
                   fy(i3[i][2]) += vy(2);
                 }, TasksPerThread(4));
    t3.Stop();

    ParallelFor (Range(gen_blocks), [&] (size_t i)
                 {
                   FlatArray<int> ind = gen_ind[i];
                   VectorMem<20> vx(ind.Size());
                   VectorMem<20> vy(ind.Size());

                   vx = fx(ind);
                   vy = s * gen_blocks[i] * vx;
                   fy(ind) += vy;
                 }, TasksPerThread(4));

    
  }

  
  
  int MySuperSparseBlockJacobi :: VHeight() const { return mat -> VWidth(); }
  int MySuperSparseBlockJacobi ::VWidth() const { return mat -> VHeight(); }
  
  AutoVector MySuperSparseBlockJacobi ::CreateRowVector () const { return mat->CreateColVector(); }
  AutoVector MySuperSparseBlockJacobi ::CreateColVector () const { return mat->CreateRowVector(); }


}
