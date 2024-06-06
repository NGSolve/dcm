#ifndef FILE_DCSCOMMON
#define FILE_DCSCOMMON

#include "intrules.hpp"

namespace ngcomp {


  template <typename T>
  size_t PosMax (const T & ar)
  {
    size_t pos = 0;
    size_t i = 0;
    for (auto & v : ar)
      {
        if (v > ar[pos]) pos = i;
        i++;
      }
    return pos;
  }

  
  template <typename T>
  inline Vec<2,T> MapQuad2Trig (Vec<2,T> xi, int vnum = 2)
  {
    // map to the first micro-hex:

    Vec<2,T> x = 0.5 * xi;
    // now (1,1) is mapping to 1/2, correct to 1/3:
    x -= 1.0/6 * xi(0)*xi(1) * Vec<2>(1,1);

    Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };

    return x(0) * verts[(vnum+1)%3] + x(1) * verts[(vnum+2)%3] + (1-x(0)-x(1)) * verts[vnum];
  }

  inline Mat<2,2> DMapQuad2Trig (Vec<2> xi, int vnum = 2)
  {
    Vec<2,AutoDiff<2>> adxi(AutoDiff<2>(xi(0),0),
                            AutoDiff<2>(xi(1),1));
    auto x = MapQuad2Trig (adxi, vnum);
    Mat<2,2> jacobi;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        jacobi(i,j) = x(i).DValue(j);
    return jacobi;
  }



  template <typename T>
  inline Vec<2,T> MapQuad2Trig (Vec<2,T> xi, int vnum, int vy)
  {
    int vx = 3-vnum-vy;
    // map to the first micro-hex:

    Vec<2,T> x = 0.5 * xi;
    // now (1,1) is mapping to 1/2, correct to 1/3:
    x -= 1.0/6 * xi(0)*xi(1) * Vec<2>(1,1);

    Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };

    return x(0) * verts[vx] + x(1) * verts[vy] + (1-x(0)-x(1)) * verts[vnum];
  }

  inline Mat<2,2> DMapQuad2Trig (Vec<2> xi, int vnum, int vy)
  {
    Vec<2,AutoDiff<2>> adxi(AutoDiff<2>(xi(0),0),
                            AutoDiff<2>(xi(1),1));
    auto x = MapQuad2Trig (adxi, vnum, vy);
    Mat<2,2> jacobi;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        jacobi(i,j) = x(i).DValue(j);
    return jacobi;
  }



  
  inline Vec<2> MapTrig2Quad (Vec<2> x)
  {
      double xi  = (3+2*x(0)-2*x(1))/2 - sqrt( sqr( (2*x(1)-2*x(0)-3)/2 ) - 6*x(0) );      
      double eta = (3+2*x(1)-2*x(0))/2 - sqrt( sqr( (2*x(0)-2*x(1)-3)/2 ) - 6*x(1) );
      return Vec<2> (xi,eta);
  }
  
  template <typename T>
  inline Vec<3,T> MapHex2Tet (Vec<3,T> xi, int vnum)
  {
    // map to the first micro-hex:

    Vec<3,T> x = 0.5 * xi;
    // fix face-midepoints:
    x(0) -= 1.0/6 * xi(0)*(xi(1)+xi(2));
    x(1) -= 1.0/6 * xi(1)*(xi(0)+xi(2));
    x(2) -= 1.0/6 * xi(2)*(xi(0)+xi(1));
    // now (1,1,1) is mapping to 1/6, correct to 1/4:

    x += 1.0/12 * xi(0)*xi(1)*xi(2) * Vec<3>(1,1,1);

    Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };

    return x(0) * verts[(vnum+1)%4] + x(1) * verts[(vnum+2)%4] + x(2) * verts[(vnum+3)%4] + (1-x(0)-x(1)-x(2)) * verts[vnum];
  }

  inline Mat<3,3> DMapHex2Tet (Vec<3> xi, int vnum)
  {
    Vec<3,AutoDiff<3>> adxi(AutoDiff<3>(xi(0),0),
                            AutoDiff<3>(xi(1),1),
                            AutoDiff<3>(xi(2),2));
    auto x = MapHex2Tet (adxi, vnum);
    Mat<3,3> jacobi;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobi(i,j) = x(i).DValue(j);
    return jacobi;
  }
  
  template <typename T>
  inline Vec<3,T> MapHex2Tet (Vec<3,T> xi)
  {
    // map to the first micro-hex:

    Vec<3,T> x = 0.5 * xi;
    // fix face-midepoints:
    x(0) -= 1.0/6 * xi(0)*(xi(1)+xi(2));
    x(1) -= 1.0/6 * xi(1)*(xi(0)+xi(2));
    x(2) -= 1.0/6 * xi(2)*(xi(0)+xi(1));
    // now (1,1,1) is mapping to 1/6, correct to 1/4:

    x += 1.0/12 * xi(0)*xi(1)*xi(2) * Vec<3>(1,1,1);
    return x;
  }

  inline Mat<3,3> DMapHex2Tet (Vec<3> xi)
  {
    Vec<3,AutoDiff<3>> adxi(AutoDiff<3>(xi(0),0),
                            AutoDiff<3>(xi(1),1),
                            AutoDiff<3>(xi(2),2));
    auto x = MapHex2Tet (adxi);
    Mat<3,3> jacobi;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobi(i,j) = x(i).DValue(j);
    return jacobi;
  }


  // new

  template <typename T>
  inline Vec<3,T> MapHex2Tet (Vec<3,T> xi, int vnum, int vz)
  {
    // map to the first micro-hex:

    Vec<3,T> x = 0.5 * xi;
    // fix face-midepoints:
    x(0) -= 1.0/6 * xi(0)*(xi(1)+xi(2));
    x(1) -= 1.0/6 * xi(1)*(xi(0)+xi(2));
    x(2) -= 1.0/6 * xi(2)*(xi(0)+xi(1));
    // now (1,1,1) is mapping to 1/6, correct to 1/4:

    x += 1.0/12 * xi(0)*xi(1)*xi(2) * Vec<3>(1,1,1);

    Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };

    // return x(0) * verts[(vnum+1)%4] + x(1) * verts[(vnum+2)%4] + x(2) * verts[(vnum+3)%4] + (1-x(0)-x(1)-x(2)) * verts[vnum];

    int vx = (vz+1)%4;
    if (vx == vnum) vx = (vx+1)%4;
    int vy = 6-vnum-vx-vz;

    int perm[] = { vx, vy, vz, vnum };
    int sign = 1;
    for (int i = 1; i < 4; i++)
      for (int j = 0; j < i; j++)
        if (perm[j] > perm[i]) sign = -sign;

    if (sign == -1) Swap (vx, vy);
    
    return x(0) * verts[vx] + x(1) * verts[vy] + x(2) * verts[vz] + (1-x(0)-x(1)-x(2)) * verts[vnum];
  }

  inline Mat<3,3> DMapHex2Tet (Vec<3> xi, int vnum, int vz)
  {
    Vec<3,AutoDiff<3>> adxi(AutoDiff<3>(xi(0),0),
                            AutoDiff<3>(xi(1),1),
                            AutoDiff<3>(xi(2),2));
    auto x = MapHex2Tet (adxi, vnum, vz);
    Mat<3,3> jacobi;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobi(i,j) = x(i).DValue(j);
    return jacobi;
  }
  
  // end new



  
  inline Vec<3> MapTet2Hex (Vec<3> x)
  {
    // Newton solve,
    Vec<3> xi = 2*x;
    for (int i = 0; i < 30; i++)
      {
        Vec<3> f = MapHex2Tet(xi, 3);
        Mat<3,3> df = DMapHex2Tet(xi, 3);
        Vec<3> w = Inv(df) * (f-x);
        xi -= w;

        // if (i > 8 && L2Norm(w) > 1e-10)
        // cerr << "newton diverges, i = " << i << ", w = " << w << endl;

        if (L2Norm(w) < 1e-12 && i < 20)
          i = 28;     // one more step
      }
    return xi;
    /*
    // Newton with homotopy
    // f = lam*MapHex2Tet(xi) + (1-lam)*0.5*xi

    Vec<3> xi = 2*x;
    int nsteps = 5;
    for (int step = 1; step <= nsteps; step++)
      {
        double lam = double(step)/nsteps;
        for (int i = 0; i < 30; i++)
          {
            Vec<3> f = lam*MapHex2Tet(xi, 3) + (1-lam)*0.5*xi;
            Mat<3,3> df = lam*DMapHex2Tet(xi, 3) + (1-lam)*0.5*Id<3>();
            Vec<3> w = Inv(df) * (f-x);
            xi -= w;
            
            if (i > 3 && L2Norm(w) > 1e-10)
              {
                cerr << "newton diverges, lam = " << lam << ", i = " << i << ", w = " << w << endl;
                cerr << "inv = " << Inv(df) << endl;
              }
          }
      }
    return xi;
    */
  }



  template <int DIM>
  class MicroJacobiDet : public CoefficientFunctionNoDerivative
  {
  public:
    MicroJacobiDet () : CoefficientFunctionNoDerivative(1,false)
    { 
      // SetDimensions(Array<int>({DIMR,DIMS}));
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & mip) const override 
    {
      IntegrationPoint ip = mip.IP();
      auto & trafo = mip.GetTransformation();
      Array<int> vnums = { 1, 2, 3 };
      trafo.GetSort(vnums);
      

      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      int minvi = (maxlam+1)%3;
      int maxvi = minvi;
      for (int i = 0; i < 3; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      
      Vec<2> x(lam[minvi], lam[maxvi]);
      Vec<2> xi = MapTrig2Quad (x);
      Mat<2,2> F = DMapQuad2Trig(xi);
      return Det(F);
    }

    /*
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override 
    {
      if (ip.DimSpace() != DIMR)
        throw Exception("illegal dim!");
      res = static_cast<const MappedIntegrationPoint<DIMS,DIMR>&>(ip).GetJacobian().AsVector();
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> res) const override 
    {
      if (ir[0].DimSpace() != DIMR)
      	throw Exception("illegal dim!");
      for (int i = 0; i < ir.Size(); i++)
      	res.Row(i).AddSize(DIMS*DIMR) = static_cast<const MappedIntegrationPoint<DIMS,DIMR>&>(ir[i]).GetJacobian().AsVector();
    }
    */
    
    /*virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override 
    {
      values.AddSize(D*D, ir.Size()) = Trans(ir.GetJacobian());
      }*/
  };
  
}
#endif
