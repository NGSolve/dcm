// ngscxx -c smallsparse.cpp ; ngsld smallsparse.o -lngstd -lngbla


#include <bla.hpp>
using namespace ngbla;


namespace ngbla
{

  void MyTranspose (SliceMatrix<> a, SliceMatrix<> b);

  
  /*
  void MyTranspose (SliceMatrix<> a, SliceMatrix<> b)
  {
    size_t j = 0;
    size_t ha = a.Height();
    size_t wa = a.Width();
    size_t da = a.Dist();
    size_t db = b.Dist();
    for ( ; j+4 <= wa; j+=4)
      {
        size_t i = 0;
        for ( ; i+4 <= ha; i+=4)
          {
            double * pa = &a(i,j);
            double * pb = &b(j,i);
            SIMD<double,4> a0(pa);
            SIMD<double,4> a1(pa+1*da);
            SIMD<double,4> a2(pa+2*da);
            SIMD<double,4> a3(pa+3*da);
            SIMD<double,4> b0, b1, b2, b3;
            SIMDTranspose(a0,a1,a2,a3, b0,b1,b2,b3);
            b0.Store(pb);
            b1.Store(pb+1*db);
            b2.Store(pb+2*db);
            b3.Store(pb+3*db);
          }
        for ( ; i < ha; i++)
          {
            double * pa = &a(i,j);
            double * pb = &b(j,i);
            pb[0] = pa[0];
            pb[db] = pa[1];
            pb[2*db] = pa[2];
            pb[3*db] = pa[3];
          }
      }
    for ( ; j < wa; j++)
      b.Row(j) = a.Col(j);
  }
  */
  
}

void Add2 (double v1, SliceVector<> x1,
           double v2, SliceVector<> x2,
           SliceVector<> y)
{
  // y = v1 * x1 + v2 * x2;
  auto delta = x2.Addr(0)-x1.Addr(0);
  double * py = y.Addr(0);
  double * px = x1.Addr(0);
  for (size_t i = 0; i < y.Size(); i++)
    {
      py[0] = v1 * px[0] + v2 * px[delta];
      px += x1.Dist();
      py += y.Dist();
    }
}

void Add2 (double v1, FlatVector<> x1,
           double v2, FlatVector<> x2,
           FlatVector<> y)
{
  // y = v1 * x1 + v2 * x2;
  auto delta = x2.Addr(0)-x1.Addr(0);
  double * py = y.Addr(0);
  double * px = x1.Addr(0);
  for (size_t i = 0; i < y.Size(); i++)
    {
      py[0] = v1 * px[0] + v2 * px[delta];
      px ++;
      py ++;
    }
}


class SmallSparseMatrix
{
  Array<int> cols;
  Array<int> rows;
  Array<int> vals;
  Array<int> index;

public:
  SmallSparseMatrix (SliceMatrix<> a, double tol = 1e-12)
  {
    index.Append(0);
    for (int i = 0; i < a.Height(); i++)
      {
        for (int j = 0; j< a.Width(); j++)
          if (fabs (a(i,j)) > tol)
            {
              vals.Append(a(i,j));
              cols.Append(j);
              rows.Append(i);   
            }
        index.Append(vals.Size());
      }
    ;
  }

  int NZE() const { return index.Last(); }

  void Mult(FlatVector<> x, FlatVector<> y)
  {
    for (int i = 0; i+1 < index.Size(); i++)
      {
        double sum = 0;
        for (auto j : Range(index[i], index[i+1]))
          sum += vals[j] * x[cols[j]];
        y[i] = sum;
      }
  }


  void Mult(SliceMatrix<> x, SliceMatrix<> y)
  {
    /*
    for (int i = 0; i < x.Height(); i++)
      Mult(x.Row(i), y.Row(i));
    */
    /*
    y = 0.0;
    for (size_t i = 0; i < vals.Size(); i++)
      y.Row(rows[i]) += vals[i] * x.Row(cols[i]);
    */
    
    for (int i = 0; i+1 < index.Size(); i++)
      {
        auto r = Range(index[i], index[i+1]);
        switch (r.Size())
          {
          case 0:
            y.Row(i) = 0.0; break;
          case 1:
            y.Row(i) = vals[r.First()] * x.Row(cols[r.First()]); break;
          case 2:
            // y.Row(i) = vals[r.First()] * x.Row(cols[r.First()]) + vals[r.First()+1] * x.Row(cols[r.First()+1]);
            Add2 (vals[r.First()], x.Row(cols[r.First()]), vals[r.First()+1], x.Row(cols[r.First()+1]), y.Row(i));            
            break;
          default:
            y.Row(i) = 0;
            for (auto j : r)
              y.Row(i) += vals[j] * x.Row(cols[j]);
          }
      }
  }
  

  // y^T = this * x^T
  void MultTrans(SliceMatrix<> x, SliceMatrix<> y)
  {
    /*
    for (int i = 0; i < x.Height(); i++)
      Mult(x.Row(i), y.Row(i));
    */

    /*
    y = 0.0;
    for (size_t i = 0; i < vals.Size(); i++)
      y.Col(rows[i]) += vals[i] * x.Col(cols[i]);
    */

    for (int i = 0; i+1 < index.Size(); i++)
      {
        auto r = Range(index[i], index[i+1]);
        switch (r.Size())
          {
          case 0:
            y.Col(i) = 0.0; break;
          case 1:
            y.Col(i) = vals[r.First()] * x.Col(cols[r.First()]);
            break;
          case 2:
            {
              // y.Col(i) = vals[r.First()] * x.Col(cols[r.First()]) + vals[r.First()+1] * x.Col(cols[r.First()+1]);
              auto yc = y.Col(i);
              double val1 = vals[r.First()];
              double val2 = vals[r.First()+1];
              auto xc1 = x.Col(cols[r.First()]);
              auto xc2 = x.Col(cols[r.First()+1]);
              // extern int myvar(); myvar();
              Add2 (val1, xc1, val2, xc2, yc);
              // yc = val1 * xc1 + val2 * xc2;
              // extern int myvar2(); myvar2();
              break;
            }
          default:
            y.Col(i) = 0;
            for (auto j : r)
              y.Col(i) += vals[j] * x.Col(cols[j]);
            ;
          }
      }

    
  }


  
  friend ostream & operator<< (ostream & ost, SmallSparseMatrix & mat);
};

ostream & operator<< (ostream & ost, SmallSparseMatrix & mat)
{
  ost << "cols = " << endl << mat.cols << endl;
  ost << "vals = " << endl << mat.vals << endl;
  ost << "index = " << endl << mat.index << endl;
  return ost;
}


// for looking into the assembly code
void Func1 (SmallSparseMatrix a, FlatVector<> x, FlatVector<> y)
{
  a.Mult (x, y);
}

void Func2 (SmallSparseMatrix a, SliceMatrix<> x, SliceMatrix<> y)
{
  a.Mult (x, y);
}

void Func3 (SmallSparseMatrix a, SliceMatrix<> x, SliceMatrix<> y)
{
  a.MultTrans (x, y);
}




int main()
{
  size_t N = 1000;
  
  Matrix a = Identity(30);
  a.Diag(1) = 2;
  Matrix b(a.Width(), N);
  Matrix c(a.Width(), N);  

  for (auto & val : b.AsVector())
    val = rand();
  // cout << "val = " << endl << b << endl;

  SmallSparseMatrix spa(a);
  // cout << "spa = " << endl << spa << endl;


  // ********* Timing mult sparse  *********
  
  {
    static Timer t("mult");
    size_t flops = spa.NZE() * N;
    size_t runs = 1e8 / flops + 1;
    
    t.Start();
    for (size_t i = 0; i < runs; i++)
      spa.Mult (b, c);
    t.Stop();
    t.AddFlops (runs * flops);
    cout << "c = " << endl << L2Norm(c) << endl;
  }

  // ********* Timing multtrans sparse  *********

  {
    Matrix bt = Trans(b);
    Matrix ct = Trans(c);
    
    static Timer t("multtrans");
    size_t flops = spa.NZE() * N;
    size_t runs = 1e8 / flops + 1;
    
    t.Start();
    for (size_t i = 0; i < runs; i++)
      spa.MultTrans (bt, ct);
    t.Stop();
    t.AddFlops (runs * flops);
    cout << "ct = " << endl << L2Norm(ct) << endl;
  }


  // ********* Timing dense *********

  {
    static Timer t("dense");
    size_t flops = a.Height()*a.Width() * N;
    size_t runs = 1e8 / flops + 1;
    
    t.Start();
    for (size_t i = 0; i < runs; i++)
      c = a*b;
    t.Stop();
    t.AddFlops (runs * flops);
    cout << "c = " << endl << L2Norm(c) << endl;
  }

  // ********* Timing Trans *********

  {
    Matrix bt = Trans(b);    
    static Timer t("transpose");
    size_t flops = b.Height() * b.Width();
    size_t runs = 1e8 / flops + 1;
    
    t.Start();
    for (size_t i = 0; i < runs; i++)
      // bt = Trans(b);
      MyTranspose (bt, b);
    t.Stop();
    t.AddFlops (runs * flops);
    cout << "b = " << endl << L2Norm(bt) << endl;
  }


  
  NgProfiler::Print(stdout);
}



/*
  JS results on Apple M1:
  
  size_t N = 1000;
  Matrix a = Identity(30);
  a.Diag(1) = 2;
  Matrix b(a.Width(), N);
  Matrix c(a.Width(), N);  

  job 8186 calls        1, time 0.0009 sec, MFlops = 12680.65 dense
  job 8187 calls        1, time 0.0087 sec, MFlops = 1155.01 multtrans
  job 8188 calls        1, time 0.0015 sec, MFlops = 6524.23 mult


  // Jul 19, 15:08   (fix multtrans, more runs)
  job 8185 calls        1, time 0.0257 sec, MFlops = 3898.16 transpose
  job 8186 calls        1, time 0.0051 sec, MFlops = 19876.87 dense
  job 8187 calls        1, time 0.0671 sec, MFlops = 1490.15 multtrans
  job 8188 calls        1, time 0.0134 sec, MFlops = 7482.05 mult
  
 */

