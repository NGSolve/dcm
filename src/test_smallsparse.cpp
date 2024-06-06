// ngscxx -c test_smallsparse.cpp ; ngsld test_smallsparse.o -lngcore -lngbla


#include <bla.hpp>
#include "smallsparse.hpp"

using namespace ngbla;


namespace ngbla
{
  // from ngsolve - ngbla
  void MyTranspose (SliceMatrix<> a, SliceMatrix<> b);
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
  a.Mult (Trans(x), Trans(y));
}




int main()
{
  
  for (size_t N = 12; N <= 5000; N *= 2)
    {

      
      Matrix<> a = Identity(30);
      a.Diag(1) = 2;
      Matrix<> b(a.Width(), N);
      Matrix<> c(a.Width(), N);  
      
      for (auto & val : b.AsVector())
        val = rand();
      // cout << "val = " << endl << b << endl;
      
      SmallSparseMatrix spa(a);
      // cout << "spa = " << endl << spa << endl;
      
      Matrix bt = Trans(b);
      Matrix ct = Trans(c);
      
      // testing
      
      spa.Mult (b, c);
      cout << "err = " << L2Norm (c - a*b) << endl;
      
      spa.Mult (Trans(bt), Trans(ct));
      cout << "err = " << L2Norm (Trans(ct) - a*b) << endl;
      
      // ********* Timing mult sparse  *********
      
      {
        Timer t("mult" + ToString(N));
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
        Timer t("multtrans" + ToString(N));
        size_t flops = spa.NZE() * N;
        size_t runs = 1e8 / flops + 1;
        
        t.Start();
        for (size_t i = 0; i < runs; i++)
          // spa.MultTrans (bt, ct);
          spa.Mult (Trans(bt), Trans(ct));
        t.Stop();
        t.AddFlops (runs * flops);
        cout << "ct = " << endl << L2Norm(ct) << endl;
      }
      
      
      // ********* Timing dense *********
      
      {
        Timer t("dense" + ToString(N));
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
        Timer t("transpose" + ToString(N));
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

  // Jul 20, 8:00    (unify MultIntern)
  job 8185 calls        1, time 0.0230 sec, MFlops = 4345.41 transpose
  job 8186 calls        1, time 0.0056 sec, MFlops = 18022.68 dense
  job 8187 calls        1, time 0.0267 sec, MFlops = 3752.93 multtrans
  job 8188 calls        1, time 0.0142 sec, MFlops = 7044.33 mult






ob 8153 calls        1, time 0.0355 sec, MFlops = 2820.27 transpose3072
job 8154 calls        1, time 0.0050 sec, MFlops = 20314.08 dense3072
job 8155 calls        1, time 0.0823 sec, MFlops = 1215.73 multtrans3072
job 8156 calls        1, time 0.0089 sec, MFlops = 11286.94 mult3072
job 8157 calls        1, time 0.0326 sec, MFlops = 3071.84 transpose1536
job 8158 calls        1, time 0.0050 sec, MFlops = 20263.62 dense1536
job 8159 calls        1, time 0.0795 sec, MFlops = 1258.31 multtrans1536
job 8160 calls        1, time 0.0087 sec, MFlops = 11485.84 mult1536
job 8161 calls        1, time 0.0235 sec, MFlops = 4258.17 transpose768
job 8162 calls        1, time 0.0049 sec, MFlops = 20376.01 dense768
job 8163 calls        1, time 0.0306 sec, MFlops = 3263.56 multtrans768
job 8164 calls        1, time 0.0090 sec, MFlops = 11162.07 mult768
job 8165 calls        1, time 0.0214 sec, MFlops = 4665.74 transpose384
job 8166 calls        1, time 0.0049 sec, MFlops = 20528.73 dense384
job 8167 calls        1, time 0.0197 sec, MFlops = 5081.28 multtrans384
job 8168 calls        1, time 0.0081 sec, MFlops = 12392.13 mult384
job 8169 calls        1, time 0.0157 sec, MFlops = 6356.14 transpose192
job 8170 calls        1, time 0.0049 sec, MFlops = 20526.10 dense192
job 8171 calls        1, time 0.0209 sec, MFlops = 4792.58 multtrans192
job 8172 calls        1, time 0.0087 sec, MFlops = 11484.59 mult192
job 8173 calls        1, time 0.0156 sec, MFlops = 6417.10 transpose96
job 8174 calls        1, time 0.0049 sec, MFlops = 20548.24 dense96
job 8175 calls        1, time 0.0190 sec, MFlops = 5254.12 multtrans96
job 8176 calls        1, time 0.0100 sec, MFlops = 9986.80 mult96
job 8177 calls        1, time 0.0154 sec, MFlops = 6509.99 transpose48
job 8178 calls        1, time 0.0050 sec, MFlops = 20032.51 dense48
job 8179 calls        1, time 0.0197 sec, MFlops = 5072.87 multtrans48
job 8180 calls        1, time 0.0127 sec, MFlops = 7851.63 mult48
job 8181 calls        1, time 0.0158 sec, MFlops = 6332.58 transpose24
job 8182 calls        1, time 0.0050 sec, MFlops = 20166.69 dense24
job 8183 calls        1, time 0.0221 sec, MFlops = 4524.74 multtrans24
job 8184 calls        1, time 0.0166 sec, MFlops = 6023.49 mult24
job 8185 calls        1, time 0.0165 sec, MFlops = 6048.49 transpose12
job 8186 calls        1, time 0.0069 sec, MFlops = 14463.14 dense12
job 8187 calls        1, time 0.0300 sec, MFlops = 3335.06 multtrans12
job 8188 calls        1, time 0.0287 sec, MFlops = 3482.50 mult12


 */

