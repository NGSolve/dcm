
namespace ngbla
{


  class SmallSparseMatrix
  {
    Array<int> cols;
    Array<int> rows;
    Array<double> vals;
    Array<int> index;




    template <bool ADD>
    static INLINE void Add1 (double v1, SliceVector<> x1,
                      SliceVector<> y) 
    {
      double * py = y.Addr(0);
      double * px = x1.Addr(0);
      size_t n = y.Size();
      for (size_t i = 0; i < n; i++, px += x1.Dist(), py += y.Dist())
        {
          if constexpr (ADD)
                         py[0] += v1 * px[0];
          else
            py[0] = v1 * px[0];
        }
    }

    template <bool ADD>
    static INLINE void Add2 (double v1, SliceVector<> x1,
                             double v2, SliceVector<> x2,
                             SliceVector<> y)
    {
      // y = v1 * x1 + v2 * x2;
      auto delta = x2.Addr(0)-x1.Addr(0);
      double * py = y.Addr(0);
      double * px = x1.Addr(0);
      size_t n = y.Size();
      for (size_t i = 0; i < n; i++, px += x1.Dist(), py += y.Dist())
        {
          if constexpr (ADD)
                         py[0] += v1 * px[0] + v2 * px[delta];
          else
            py[0] = v1 * px[0] + v2 * px[delta];
        }
    }

    template <bool ADD>
    static INLINE void Add3 (double v1, SliceVector<> x1,
                             double v2, SliceVector<> x2,
                             double v3, SliceVector<> x3,
                             SliceVector<> y)
    {
      auto delta2 = x2.Addr(0)-x1.Addr(0);
      auto delta3 = x3.Addr(0)-x1.Addr(0);
      double * py = y.Addr(0);
      double * px = x1.Addr(0);
      size_t n = y.Size();
      for (size_t i = 0; i < n; i++, px += x1.Dist(), py += y.Dist())
        {
          if constexpr (ADD)
                         py[0] += v1 * px[0] + v2 * px[delta2] + v3 * px[delta3];
          else
            py[0] = v1 * px[0] + v2 * px[delta2] + v3*px[delta3];
        }
    }

    template <bool ADD>
    static INLINE void Add4 (double v1, SliceVector<> x1,
                             double v2, SliceVector<> x2,
                             double v3, SliceVector<> x3,
                             double v4, SliceVector<> x4,
                             SliceVector<> y)
    {
      auto delta2 = x2.Addr(0)-x1.Addr(0);
      auto delta3 = x3.Addr(0)-x1.Addr(0);
      auto delta4 = x4.Addr(0)-x1.Addr(0);
      double * py = y.Addr(0);
      double * px = x1.Addr(0);
      size_t n = y.Size();
      for (size_t i = 0; i < n; i++, px += x1.Dist(), py += y.Dist())
        {
          if constexpr (ADD)
                         py[0] += v1 * px[0] + v2 * px[delta2] + v3 * px[delta3] + v4 * px[delta4];
          else
            py[0] = v1 * px[0] + v2 * px[delta2] + v3*px[delta3] + v4*px[delta4];
        }
    }


    template <bool ADD>
    static INLINE void Add1 (double v1, FlatVector<> x1,
                             FlatVector<> y)
    {
      // y = v1 * x1 
      auto delta1 = x1.Addr(0)-y.Addr(0);
      double * py = y.Addr(0);
      size_t n = y.Size();
      size_t i = 0;
      constexpr int SW = 2*SIMD<double>::Size();
      for ( ; i+SW <= n; i+=SW, py+=SW)
        {
          SIMD<double,SW> x1(py+delta1);
          SIMD<double,SW> y;
          if constexpr (ADD)
                         y = SIMD<double,SW>(py);
          else
            y = SIMD<double,SW>(0.0);
          y += v1 * x1;
          y.Store(py);
        }

      if constexpr (ADD)
                     for (size_t i = 0; i < (n%SW); i++, py++)
                       py[0] += v1 * py[delta1];
      else
        for (size_t i = 0; i < (n%SW); i++, py++)
          py[0] = v1 * py[delta1];    
    }

    template <bool ADD>
    static INLINE void Add2 (double v1, FlatVector<> x1,
                             double v2, FlatVector<> x2,
                             FlatVector<> y)
    {
      // y = v1 * x1 + v2 * x2;
      auto delta1 = x1.Addr(0)-y.Addr(0);
      auto delta2 = x2.Addr(0)-y.Addr(0);
      double * py = y.Addr(0);
      size_t n = y.Size();
      size_t i = 0;
      constexpr int SW = 2*SIMD<double>::Size();
      for ( ; i+SW <= n; i+=SW, py+=SW)
        {
          SIMD<double,SW> x1(py+delta1);
          SIMD<double,SW> x2(py+delta2);
          SIMD<double,SW> y;
          if constexpr (ADD)
                         y = SIMD<double,SW>(py) + v1 * x1 + v2 * x2;
          else
            y = v1 * x1 + v2 * x2; 
          y.Store(py);
        }

      if constexpr (ADD)  
                     for (size_t i = 0; i < (n%SW); i++, py++)
                       py[0] += v1 * py[delta1] + v2 * py[delta2];
      else
        for (size_t i = 0; i < (n%SW); i++, py++)
          py[0] = v1 * py[delta1] + v2 * py[delta2];
    }

    template <bool ADD>
    static INLINE void Add3 (double v1, FlatVector<> x1,
                             double v2, FlatVector<> x2,
                             double v3, FlatVector<> x3,
                             FlatVector<> y)
    {
      // y = v1 * x1 + v2 * x2;
      auto delta1 = x1.Addr(0)-y.Addr(0);
      auto delta2 = x2.Addr(0)-y.Addr(0);
      auto delta3 = x3.Addr(0)-y.Addr(0);
      double * py = y.Addr(0);
      size_t n = y.Size();
      size_t i = 0;
      constexpr int SW = 2*SIMD<double>::Size();
      for ( ; i+SW <= n; i+=SW, py+=SW)
        {
          SIMD<double,SW> x1(py+delta1);
          SIMD<double,SW> x2(py+delta2);
          SIMD<double,SW> x3(py+delta3);
          SIMD<double,SW> y;
          if constexpr (ADD)
                         y = SIMD<double,SW>(py) + v1 * x1 + v2 * x2 + v3 * x3;
          else
            y = v1 * x1 + v2 * x2 + v3 * x3; 
          y.Store(py);
        }

      if constexpr (ADD)  
                     for (size_t i = 0; i < (n%SW); i++, py++)
                       py[0] += v1 * py[delta1] + v2 * py[delta2] + v3 * py[delta3];
      else
        for (size_t i = 0; i < (n%SW); i++, py++)
          py[0] = v1 * py[delta1] + v2 * py[delta2] + v3 * py[delta3];
    }


    template <bool ADD>
    static INLINE void Add4 (double v1, FlatVector<> x1,
                             double v2, FlatVector<> x2,
                             double v3, FlatVector<> x3,
                             double v4, FlatVector<> x4,
                             FlatVector<> y)
    {
      // y = v1 * x1 + v2 * x2;
      auto delta1 = x1.Addr(0)-y.Addr(0);
      auto delta2 = x2.Addr(0)-y.Addr(0);
      auto delta3 = x3.Addr(0)-y.Addr(0);
      auto delta4 = x4.Addr(0)-y.Addr(0);
      double * py = y.Addr(0);
      size_t n = y.Size();
      size_t i = 0;
      constexpr int SW = 2*SIMD<double>::Size();
      for ( ; i+SW <= n; i+=SW, py+=SW)
        {
          SIMD<double,SW> x1(py+delta1);
          SIMD<double,SW> x2(py+delta2);
          SIMD<double,SW> x3(py+delta3);
          SIMD<double,SW> x4(py+delta4);
          SIMD<double,SW> y;
          if constexpr (ADD)
                         y = SIMD<double,SW>(py) + v1 * x1 + v2 * x2 + v3 * x3 + v4 * x4;
          else
            y = v1 * x1 + v2 * x2 + v3 * x3 + v4 * x4; 
          y.Store(py);
        }

      if constexpr (ADD)  
                     for (size_t i = 0; i < (n%SW); i++, py++)
                       py[0] += v1 * py[delta1] + v2 * py[delta2] + v3 * py[delta3] + v4 * py[delta4];
      else
        for (size_t i = 0; i < (n%SW); i++, py++)
          py[0] = v1 * py[delta1] + v2 * py[delta2] + v3 * py[delta3] + v4 * py[delta4];
    }




    
  public:
    SmallSparseMatrix() = default;
    
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

    template <ORDERING ORDER>
    void MultIntern(SliceMatrix<double, ORDER> x, SliceMatrix<double, ORDER> y) const
    {
      /*
        y = 0.0;
        for (size_t i = 0; i < vals.Size(); i++)
        y.Row(rows[i]) += vals[i] * x.Row(cols[i]);
      */
    
      for (size_t i = 0; i+1 < index.Size(); i++)
        {
          auto r = Range(index[i], index[i+1]);
          switch (r.Size())
            {
            case 0:
              y.Row(i) = 0.0; break;
              
            case 1:
              Add1<false> (vals[r.First()], x.Row(cols[r.First()]),
                           y.Row(i));
              break;
              
            case 2:
              Add2<false> (vals[r.First()], x.Row(cols[r.First()]),
                           vals[r.First()+1], x.Row(cols[r.First()+1]),
                           y.Row(i));
              break;

            case 3:
              Add3<false> (vals[r.First()], x.Row(cols[r.First()]),
                           vals[r.First()+1], x.Row(cols[r.First()+1]),
                           vals[r.First()+2], x.Row(cols[r.First()+2]),                           
                           y.Row(i));
              break;
              
            default:
              Add4<false> (vals[r.First()], x.Row(cols[r.First()]),
                           vals[r.First()+1], x.Row(cols[r.First()+1]),
                           vals[r.First()+2], x.Row(cols[r.First()+2]),
                           vals[r.First()+3], x.Row(cols[r.First()+3]),
                           y.Row(i));

              for (auto j : r.Modify(4,0))
                Add1<true> (vals[j], x.Row(cols[j]), y.Row(i));
            }
        }
    }
  
    void Mult(SliceMatrix<double> x, SliceMatrix<double> y) const
    {
      MultIntern (x,y);
    }
    void Mult(SliceMatrix<double,ColMajor> x, SliceMatrix<double,ColMajor> y) const
    {
      MultIntern (x,y);
    }




    template <ORDERING ORDER>
    void MultAddIntern(double val, SliceMatrix<double, ORDER> x, SliceMatrix<double, ORDER> y) const
    {
      for (size_t i = 0; i+1 < index.Size(); i++)
        {
          auto r = Range(index[i], index[i+1]);
          switch (r.Size())
            {
            case 0:
              break;
              
            case 1:
              Add1<true> (val*vals[r.First()], x.Row(cols[r.First()]),
                           y.Row(i));
              break;
              
            case 2:
              Add2<true> (val*vals[r.First()], x.Row(cols[r.First()]),
                          val*vals[r.First()+1], x.Row(cols[r.First()+1]),
                          y.Row(i));
              break;
              
            case 3:
              Add3<true> (val*vals[r.First()], x.Row(cols[r.First()]),
                          val*vals[r.First()+1], x.Row(cols[r.First()+1]),
                          val*vals[r.First()+2], x.Row(cols[r.First()+2]),
                          y.Row(i));
              break;

            default:
              Add4<true> (val*vals[r.First()], x.Row(cols[r.First()]),
                          val*vals[r.First()+1], x.Row(cols[r.First()+1]),
                          val*vals[r.First()+2], x.Row(cols[r.First()+2]),
                          val*vals[r.First()+3], x.Row(cols[r.First()+3]),                          
                          y.Row(i));

              for (auto j : r.Modify(4,0))
                Add1<true> (val*vals[j], x.Row(cols[j]), y.Row(i));
            }
        }
    }
  
    void MultAdd (double val, SliceMatrix<double> x, SliceMatrix<double> y) const
    {
      MultAddIntern (val, x, y);
    }
    void MultAdd (double val, SliceMatrix<double,ColMajor> x, SliceMatrix<double,ColMajor> y) const
    {
      MultAddIntern (val, x, y);
    }



    
#ifdef HISTORIC
    // y^T = this * x^T
    // use  spmat.Mult (Trans(b), Trans(c)) instead
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
                Add2<false> (vals[r.First()], x.Col(cols[r.First()]),
                             vals[r.First()+1], x.Col(cols[r.First()+1]),
                             y.Col(i));              
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
#endif
    

  
    friend ostream & operator<< (ostream & ost, SmallSparseMatrix & mat);
  };

  inline ostream & operator<< (ostream & ost, SmallSparseMatrix & mat)
  {
    ost << "cols = " << endl << mat.cols << endl;
    ost << "vals = " << endl << mat.vals << endl;
    ost << "index = " << endl << mat.index << endl;
    return ost;
  }

}
