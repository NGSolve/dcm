#include <comp.hpp>
#include "hcurldualcells.hpp"
#include "h1primalcells.hpp"
#include "hcurlprimalcells.hpp"
#include "dcs_common.hpp" 
#include "supersparse.hpp"
#include "smallsparse.hpp"


namespace ngcomp
{
  
  class HCurlDualCellTrig : public HCurlCellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIR;
    const IntegrationRule & tangentialIR;
  public:
    HCurlDualCellTrig (const IntegrationRule & _GaussRadauIR,
                       const IntegrationRule & _tangentialIR)
      : HCurlCellFiniteElement<2> (3*2*_GaussRadauIR.Size()*_tangentialIR.Size(), _GaussRadauIR.Size()-1),
      GaussRadauIR(_GaussRadauIR), tangentialIR(_tangentialIR)
    { ; }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    {
      //cout << "HCurlDualCellTrig.HCurlDualCellTrig called" <<endl;

      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.AddSize(ndof, 2) = 0;
      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;

      if (vnums[minv]>vnums[maxv])
        {
          minv = (maxlam+2)%3;
          maxv = (maxlam+1)%3;
        }

      //cout << "maxlam: local: " << maxlam <<  ", global " << vnums[maxlam] << endl;
      //cout << "minv: local: " << minv <<  ", global " << vnums[minv] << endl;
      //cout << "maxv: local: " << maxv <<  ", global " << vnums[maxv] << endl;
      Vec<2> x(lam[minv],lam[maxv]);
      Vec<2> xi = MapTrig2Quad(x);

      Mat<2,2> F = DMapQuad2Trig(xi);
      Mat<2,2> F2;

      Vec<2> verts[] = { Vec<2>(1, 0 ), Vec<2>( 0, 1 ), Vec<2>( 0, 0)  };      
      F2.Col(0) = verts[minv]-verts[maxlam];
      F2.Col(1) = verts[maxv]-verts[maxlam];

      Mat<2> trafo = Trans(Inv(F2*F));


      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);

      auto assign_shape = [&](int nr, int ix, int iy,int dir)
        {
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafo*Vec<2>(polxitang[ix]*poleta[iy],
                                           0);
              //shape.Row(nr) = Vec<2>(polxitang[ix]*poleta[iy],0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafo*Vec<2>(0,
                                           polxi[ix]*poletatang[iy]);
              //shape.Row(nr) = Vec<2>(0,polxi[ix]*poletatang[iy]);
              break;
            }
            break;
            }
        };

      
      int ii = 0;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          //cout << "edge " << i << endl;
          for (int j = 0; j < 2; j++)
            {
              //cout << "vertex " << e[j] <<endl;
              if (e[j] == maxlam)
                {
                  if (e[1-j] == minv)
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, k, 0, 0);
                    }
                  else
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, 0, k, 1);
                    }
                }
              ii+=ndt;
            }
        }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
        {

          if (f[i] == maxlam)
            {
              int v1 = f[(i+1)%3];
              int kk = ii;
              for (int l = 1; l < ndGR; l++)
                {
                  for (int k = 0; k < ndt; k++)
                    {
                      if (v1 == minv)
                        {
                          assign_shape(kk++, k, l,0);
                          assign_shape(kk++, l, k,1);
                        }
                      else 
                        {
                          assign_shape(kk++, l, k,1);
                          assign_shape(kk++, k, l,0);
                        }

                    }
                }
            }
          ii += 2*(ndGR-1)*ndt;
        }
      //cout << "CalcShape done, shape = " << shape << endl;
     
    }
    

    
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1)};
      int maxlam = PosMax(lam);

      curlshape = 0;

      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;
      if (vnums[minv] > vnums[maxv]) {
        minv = (maxlam+2)%3;
        maxv = (maxlam+1)%3;
      }

      
      Vec<2> x(lam[minv],lam[maxv]);
      Vec<2> xi = MapTrig2Quad(x);

      Mat<2,2> F = DMapQuad2Trig(xi);
      Mat<2,2> F2;

      Vec<2> verts[] = { Vec<2>(1, 0 ), Vec<2>( 0, 1 ), Vec<2>( 0, 0)  };      
      F2.Col(0) = verts[minv]-verts[maxlam];
      F2.Col(1) = verts[maxv]-verts[maxlam];

      //Mat<2> trafo = Trans(Inv(F2*F));
      double trafo = 1.0/Det(F2*F); 

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<AutoDiff<1>, 20> polxi(ndGR), poleta(ndGR);
      ArrayMem<AutoDiff<1>, 20> polxitang(ndt), poletatang(ndt);

      LagrangePolynomials(AutoDiff<1>(xi(0), 0), GaussRadauIR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1), 0), GaussRadauIR, poleta);

      LagrangePolynomials(AutoDiff<1>(xi(0), 0), tangentialIR, polxitang);
      LagrangePolynomials(AutoDiff<1>(xi(1), 0), tangentialIR, poletatang);

      auto assign =  [&](int nr, int ix, int iy, int dir)
        {
          switch (dir)
            {
            case 0:
            {
              curlshape.Row(nr) = trafo*(-polxitang[ix].Value()*poleta[iy].DValue(0));
              break;
            }
            case 1:
            {
              curlshape.Row(nr) = trafo*(polxi[ix].DValue(0)*poletatang[iy].Value());
              break;
            }
            break;
            }
        };


      int ii = 0;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          //cout << "edge " << i << endl;
          for (int j = 0; j < 2; j++)
            {
              //cout << "vertex " << e[j] <<endl;
              if (e[j] == maxlam)
                {
                  if (e[1-j] == minv)
                    {
                      for (int k = 0; k < ndt; k++)
                        assign(ii+k, k, 0, 0);
                    }
                  else
                    {
                      for (int k = 0; k < ndt; k++)
                        assign(ii+k, 0, k, 1);
                    }
                }
              ii+=ndt;
            }
        }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
        {

          if (f[i] == maxlam)
            {
              int v1 = f[(i+1)%3];
              int kk = ii;
              for (int l = 1; l < ndGR; l++)
                {
                  for (int k = 0; k < ndt; k++)
                    {
                      if (v1 == minv)
                        {
                          assign(kk++, k, l,0);
                          assign(kk++, l, k,1);
                        }
                      else 
                        {
                          assign(kk++, l, k,1);
                          assign(kk++, k, l,0);
                        }

                    }
                }
            }
          ii += 2*(ndGR-1)*ndt;
        }
        
    }
    virtual void CalcPiolaAltShape (const IntegrationPoint & ip, 
                                    BareSliceMatrix<> shape) const override
    {
      //cout << "HCurlDualCellTrig.HCurlDualCellTrig called" <<endl;

      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      //shape = 0;
      //
      shape.AddSize(ndof, 2) = 0;
      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;

      if (vnums[minv]>vnums[maxv])
        {
          minv = (maxlam+2)%3;
          maxv = (maxlam+1)%3;
        }

      //cout << "maxlam: local: " << maxlam <<  ", global " << vnums[maxlam] << endl;
      //cout << "minv: local: " << minv <<  ", global " << vnums[minv] << endl;
      //cout << "maxv: local: " << maxv <<  ", global " << vnums[maxv] << endl;
      Vec<2> x(lam[minv],lam[maxv]);
      Vec<2> xi = MapTrig2Quad(x);

      //Mat<2,2> F = DMapQuad2Trig(xi);
      Mat<2,2> F2;

      Vec<2> verts[] = { Vec<2>(1, 0 ), Vec<2>( 0, 1 ), Vec<2>( 0, 0)  };      
      F2.Col(0) = verts[minv]-verts[maxlam];
      F2.Col(1) = verts[maxv]-verts[maxlam];

      //Mat<2> trafo = Trans(Inv(F2*F));


      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);

      auto assign_shape = [&](int nr, int ix, int iy,int dir)
        {
          Vec<2> xinode(GaussRadauIR[ix](0), GaussRadauIR[iy](0));

          Mat<2,2> F = DMapQuad2Trig(xinode);

          Mat<2> trafo = 1/Det(F2*F)* (F2*F);
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafo*Vec<2>(polxitang[ix]*poleta[iy],
                                           0);
              //shape.Row(nr) = Vec<2>(polxitang[ix]*poleta[iy],0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafo*Vec<2>(0,
                                           polxi[ix]*poletatang[iy]);
              //shape.Row(nr) = Vec<2>(0,polxi[ix]*poletatang[iy]);
              break;
            }
            break;
            }
        };

      
      int ii = 0;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          //cout << "edge " << i << endl;
          for (int j = 0; j < 2; j++)
            {
              //cout << "vertex " << e[j] <<endl;
              if (e[j] == maxlam)
                {
                  if (e[1-j] == minv)
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, k, 0, 0);
                    }
                  else
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, 0, k, 1);
                    }
                }
              ii+=ndt;
            }
        }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
        {

          if (f[i] == maxlam)
            {
              int v1 = f[(i+1)%3];
              int kk = ii;
              for (int l = 1; l < ndGR; l++)
                {
                  for (int k = 0; k < ndt; k++)
                    {
                      if (v1 == minv)
                        {
                          assign_shape(kk++, k, l,0);
                          assign_shape(kk++, l, k,1);
                        }
                      else 
                        {
                          assign_shape(kk++, l, k,1);
                          assign_shape(kk++, k, l,0);
                        }

                    }
                }
            }
          ii += 2*(ndGR-1)*ndt;
        }
      //cout << "CalcShape done, shape = " << shape << endl;
     
    }

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const override
    {
      //cout << "HCurlDualCellTrig.HCurlDualCellTrig called" <<endl;

      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 2) = 0;
      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;

      if (vnums[minv]>vnums[maxv])
        {
          minv = (maxlam+2)%3;
          maxv = (maxlam+1)%3;
        }

      //cout << "maxlam: local: " << maxlam <<  ", global " << vnums[maxlam] << endl;
      //cout << "minv: local: " << minv <<  ", global " << vnums[minv] << endl;
      //cout << "maxv: local: " << maxv <<  ", global " << vnums[maxv] << endl;
      Vec<2> x(lam[minv],lam[maxv]);
      Vec<2> xi = MapTrig2Quad(x);

      Mat<2,2> F = DMapQuad2Trig(xi);
      Mat<2,2> F2;

      Vec<2> verts[] = { Vec<2>(1, 0 ), Vec<2>( 0, 1 ), Vec<2>( 0, 0)  };      
      F2.Col(0) = verts[minv]-verts[maxlam];
      F2.Col(1) = verts[maxv]-verts[maxlam];

      Mat<2> trafo = Trans(Inv(F2*F));


      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);

      auto assign_shape = [&](int nr, int ix, int iy,int dir)
        {
          Vec<2> xinode(GaussRadauIR[ix](0), GaussRadauIR[iy](0));

          Mat<2,2> F = DMapQuad2Trig(xinode);

          Mat<2> trafo = Trans(Inv(F2*F));
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafo*Vec<2>(polxitang[ix]*poleta[iy],
                                           0);
              //shape.Row(nr) = Vec<2>(polxitang[ix]*poleta[iy],0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafo*Vec<2>(0,
                                           polxi[ix]*poletatang[iy]);
              //shape.Row(nr) = Vec<2>(0,polxi[ix]*poletatang[iy]);
              break;
            }
            break;
            }
        };

      
      int ii = 0;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          //cout << "edge " << i << endl;
          for (int j = 0; j < 2; j++)
            {
              //cout << "vertex " << e[j] <<endl;
              if (e[j] == maxlam)
                {
                  if (e[1-j] == minv)
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, k, 0, 0);
                    }
                  else
                    {
                      for (int k = 0; k < ndt; k++)
                        assign_shape(ii+k, 0, k, 1);
                    }
                }
              ii+=ndt;
            }
        }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
        {

          if (f[i] == maxlam)
            {
              int v1 = f[(i+1)%3];
              int kk = ii;
              for (int l = 1; l < ndGR; l++)
                {
                  for (int k = 0; k < ndt; k++)
                    {
                      if (v1 == minv)
                        {
                          assign_shape(kk++, k, l,0);
                          assign_shape(kk++, l, k,1);
                        }
                      else 
                        {
                          assign_shape(kk++, l, k,1);
                          assign_shape(kk++, k, l,0);
                        }

                    }
                }
            }
          ii += 2*(ndGR-1)*ndt;
        }
      //cout << "CalcShape done, shape = " << shape << endl;
     
    }

    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override
    {
      Matrix mat(ndof, ndof);
      Vector rhs(ndof);

      mat = 0.0;
      rhs = 0.0;
      IntegrationRule ir(ET_SEGM, order+1);
      
      auto [pnts, inds, nump] = GetMicroCellPoints<ET_TRIG> (order+1, true);
      auto [edges, nume] = GetMicroCellEdges<ET_TRIG> (order+1, true);

      Vec<2> fi;
      Matrix shape(ndof, 2);
      
      for (int i = 0; i < nume; i++) // sum over real edges, only 
        {
          Vec<2> p0 = pnts[edges[i][0]];
          Vec<2> p1 = pnts[edges[i][1]];

          for (auto ip : ir)
            {
              double t = ip(0);
              Vec<2> tauref = p1-p0;
              Vec<2> x = p0+t*(p1-p0);

              IntegrationPoint ip2d(x(0), x(1), 0, 0);
              MappedIntegrationPoint<2,2> mip(ip2d, trafo);
              Vec<2> tau = mip.GetJacobian() * tauref;
              
              func.Evaluate (mip, fi);

              CalcShape (ip2d, shape);

              mat.Row(i) += ip.Weight() * shape*tauref;
              rhs(i) += ip.Weight() * InnerProduct(fi,tau);
            }
        }

      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }
  };


  class HCurlDualCellTet : public HCurlCellFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIR;
    const IntegrationRule & tangentialIR;
  public:
    HCurlDualCellTet (const IntegrationRule & _GaussRadauIR,
                      const IntegrationRule & _tangentialIR)
      : HCurlCellFiniteElement<3> (4*3*sqr(_GaussRadauIR.Size())*_tangentialIR.Size(), _GaussRadauIR.Size()-1),
      GaussRadauIR(_GaussRadauIR), tangentialIR(_tangentialIR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = Trans(Inv(F2*F));

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt), polzetatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);
      LagrangePolynomials(xi(2), tangentialIR, polzetatang);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafo*Vec<3>(polxitang[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafo*Vec<3>(0, polxi[i[0]]*poletatang[i[1]]*polzeta[i[2]], 0);
              break;
            }
            case 2:
            {
              shape.Row(nr) = trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzetatang[i[2]]);
              break;
            }
            break;
            }
        };


        
      
      int ii = 0;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndt; k++)
                    {
                      ind[dirv2] = k;
                      assign(ii+k, ind, dirv2);
                    }
                }
              ii += ndt;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR; l++)
                    for (int k = 0; k < ndt; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind, dirv1);
                        ind[dirv1] = l;
                        ind[dirv2] = k;
                        assign(kk++, ind, dirv2);
                      }
                }
              ii += 2*ndt*(ndGR-1);
            }
        }
      
      ii += maxlam*3*sqr(ndGR-1)*ndt;
      for (int i = 1; i < ndGR; i++)
        for (int j = 1; j < ndGR; j++)
          for (int k = 0; k < ndt; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }
    }

    
    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      curlshape = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      // Mat<3> trafo = Trans(Inv(F2*F));
      Mat<3> trafo = 1.0/Det(F2*F) * (F2*F);

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<AutoDiff<1>, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      ArrayMem<AutoDiff<1>, 20> polxitang(ndt), poletatang(ndt), polzetatang(ndt);         
      LagrangePolynomials(AutoDiff<1>(xi(0), 0), GaussRadauIR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1), 0), GaussRadauIR, poleta);
      LagrangePolynomials(AutoDiff<1>(xi(2), 0), GaussRadauIR, polzeta);

      LagrangePolynomials(AutoDiff<1>(xi(0), 0), tangentialIR, polxitang);
      LagrangePolynomials(AutoDiff<1>(xi(1), 0), tangentialIR, poletatang);
      LagrangePolynomials(AutoDiff<1>(xi(2), 0), tangentialIR, polzetatang);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          switch (dir)
            {
            case 0:
            {
              curlshape.Row(nr) = trafo*Vec<3>(0,
                                               polxitang[i[0]].Value()*poleta[i[1]].Value()*polzeta[i[2]].DValue(0),
                                               -polxitang[i[0]].Value()*poleta[i[1]].DValue(0)*polzeta[i[2]].Value());
              break;
            }
            case 1:
            {
              curlshape.Row(nr) = trafo*Vec<3>(-polxi[i[0]].Value()*poletatang[i[1]].Value()*polzeta[i[2]].DValue(0),
                                               0,
                                               polxi[i[0]].DValue(0)*poletatang[i[1]].Value()*polzeta[i[2]].Value());
              break;
            }
            case 2:
            {
              curlshape.Row(nr) = trafo*Vec<3>(polxi[i[0]].Value()*poleta[i[1]].DValue(0)*polzetatang[i[2]].Value(),
                                               -polxi[i[0]].DValue(0)*poleta[i[1]].Value()*polzetatang[i[2]].Value(),
                                               0);
              break;
            }
            break;
            }
        };


        
      
      int ii = 0;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndt; k++)
                    {
                      ind[dirv2] = k;
                      assign(ii+k, ind, dirv2);
                    }
                }
              ii += ndt;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR; l++)
                    for (int k = 0; k < ndt; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind, dirv1);
                        ind[dirv1] = l;
                        ind[dirv2] = k;
                        assign(kk++, ind, dirv2);
                      }
                }
              ii += 2*ndt*(ndGR-1);
            }
        }
      
      ii += maxlam*3*sqr(ndGR-1)*ndt;
      for (int i = 1; i < ndGR; i++)
        for (int j = 1; j < ndGR; j++)
          for (int k = 0; k < ndt; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }
    }




    virtual void CalcAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = Trans(Inv(F2*F));

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt), polzetatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);
      LagrangePolynomials(xi(2), tangentialIR, polzetatang);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          Vec<3> xinode(GaussRadauIR[i[0]](0), GaussRadauIR[i[1]](0), GaussRadauIR[i[2]](0));
          Mat<3,3> F = DMapHex2Tet(xinode);          
          Mat<3> trafo = Trans(Inv(F2*F));
          
          
          switch (dir)
            {
            case 0:
              {
                shape.Row(nr) = trafo*Vec<3>(polxitang[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
                break;
              }
            case 1:
              {
                shape.Row(nr) = trafo*Vec<3>(0, polxi[i[0]]*poletatang[i[1]]*polzeta[i[2]], 0);
                break;
              }
            case 2:
              {
                shape.Row(nr) = trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzetatang[i[2]]);
                break;
              }
              break;
            }
        };


        
      
      int ii = 0;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndt; k++)
                    {
                      ind[dirv2] = k;
                      assign(ii+k, ind, dirv2);
                    }
                }
              ii += ndt;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR; l++)
                    for (int k = 0; k < ndt; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind, dirv1);
                        ind[dirv1] = l;
                        ind[dirv2] = k;
                        assign(kk++, ind, dirv2);
                      }
                }
              ii += 2*ndt*(ndGR-1);
            }
        }
      
      ii += maxlam*3*sqr(ndGR-1)*ndt;
      for (int i = 1; i < ndGR; i++)
        for (int j = 1; j < ndGR; j++)
          for (int k = 0; k < ndt; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }
    }

  virtual void CalcPiolaShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      // Mat<3> trafo = Trans(Inv(F2*F));
      Mat<3> trafox = 1/Det(F2*F) * (F2*F);

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt), polzetatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);
      LagrangePolynomials(xi(2), tangentialIR, polzetatang);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          Vec<3> xinode(GaussRadauIR[i[0]](0), GaussRadauIR[i[1]](0), GaussRadauIR[i[2]](0));
          Mat<3,3> F = DMapHex2Tet(xinode);          
          Mat<3> trafo1 = Det(F2*F) *Inv(F2*F) * Trans(Inv(F2*F));
          Mat<3> trafo = trafox * trafo1;
          
          switch (dir)
            {
            case 0:
              {
                shape.Row(nr) = trafo*Vec<3>(polxitang[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
                break;
              }
            case 1:
              {
                shape.Row(nr) = trafo*Vec<3>(0, polxi[i[0]]*poletatang[i[1]]*polzeta[i[2]], 0);
                break;
              }
            case 2:
              {
                shape.Row(nr) = trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzetatang[i[2]]);
                break;
              }
              break;
            }
        };


        
      
      int ii = 0;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndt; k++)
                    {
                      ind[dirv2] = k;
                      assign(ii+k, ind, dirv2);
                    }
                }
              ii += ndt;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR; l++)
                    for (int k = 0; k < ndt; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind, dirv1);
                        ind[dirv1] = l;
                        ind[dirv2] = k;
                        assign(kk++, ind, dirv2);
                      }
                }
              ii += 2*ndt*(ndGR-1);
            }
        }
      
      ii += maxlam*3*sqr(ndGR-1)*ndt;
      for (int i = 1; i < ndGR; i++)
        for (int j = 1; j < ndGR; j++)
          for (int k = 0; k < ndt; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }
    }

  virtual void CalcPiolaAltShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      // Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      //Mat<3> trafox = 1/Det(F2*F) * (F2*F);

      int ndGR = GaussRadauIR.Size();
      int ndt = tangentialIR.Size();

      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      ArrayMem<double, 20> polxitang(ndt), poletatang(ndt), polzetatang(ndt);         
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);

      LagrangePolynomials(xi(0), tangentialIR, polxitang);
      LagrangePolynomials(xi(1), tangentialIR, poletatang);
      LagrangePolynomials(xi(2), tangentialIR, polzetatang);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          Vec<3> xinode(GaussRadauIR[i[0]](0), GaussRadauIR[i[1]](0), GaussRadauIR[i[2]](0));
          Mat<3,3> F = DMapHex2Tet(xinode);          
          Mat<3> trafo = 1/Det(F2*F) * (F2*F);
          
          switch (dir)
            {
            case 0:
              {
                shape.Row(nr) = trafo*Vec<3>(polxitang[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
                break;
              }
            case 1:
              {
                shape.Row(nr) = trafo*Vec<3>(0, polxi[i[0]]*poletatang[i[1]]*polzeta[i[2]], 0);
                break;
              }
            case 2:
              {
                shape.Row(nr) = trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzetatang[i[2]]);
                break;
              }
              break;
            }
        };


        
      
      int ii = 0;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndt; k++)
                    {
                      ind[dirv2] = k;
                      assign(ii+k, ind, dirv2);
                    }
                }
              ii += ndt;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR; l++)
                    for (int k = 0; k < ndt; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind, dirv1);
                        ind[dirv1] = l;
                        ind[dirv2] = k;
                        assign(kk++, ind, dirv2);
                      }
                }
              ii += 2*ndt*(ndGR-1);
            }
        }
      
      ii += maxlam*3*sqr(ndGR-1)*ndt;
      for (int i = 1; i < ndGR; i++)
        for (int j = 1; j < ndGR; j++)
          for (int k = 0; k < ndt; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }
    }


    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    {
      Matrix mat(ndof, ndof);
      Vector rhs(ndof);

      mat = 0.0;
      rhs = 0.0;
      IntegrationRule ir(ET_SEGM, order+1);
      
      auto [pnts, inds, nump] = GetMicroCellPoints<ET_TET> (order+1, true);
      auto [edges, nume] = GetMicroCellEdges<ET_TET> (order+1, true);

      Vec<3> fi;
      Matrix shape(ndof, 3);
      
      for (int i = 0; i < nume; i++) // sum over real edges, only 
        {
          Vec<3> p0 = pnts[edges[i][0]];
          Vec<3> p1 = pnts[edges[i][1]];

          for (auto ip : ir)
            {
              double t = ip(0);
              Vec<3> tauref = p1-p0;
              Vec<3> x = p0+t*(p1-p0);

              IntegrationPoint ip3d(x(0), x(1), x(2), 0);
              MappedIntegrationPoint<3,3> mip(ip3d, trafo);
              Vec<3> tau = mip.GetJacobian() * tauref;
              
              func.Evaluate (mip, fi);

              CalcShape (ip3d, shape);

              mat.Row(i) += ip.Weight() * shape*tauref;
              rhs(i) += ip.Weight() * InnerProduct(fi,tau);
            }
        }

      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }
  };



  


  class HCurlDualCellPotentialTet : public ScalarFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIRplus;
  public:
    HCurlDualCellPotentialTet (const IntegrationRule & _GaussRadauIRplus)
      : ScalarFiniteElement<3> (4 + 6         // vertex + edge MP
                                + 12*(_GaussRadauIRplus.Size()-2)   // edges
                                + 12*sqr(_GaussRadauIRplus.Size()-2) // face-quad inner
                                + 12*(_GaussRadauIRplus.Size()-2)    // edges inside face
                                + 4*sqr(_GaussRadauIRplus.Size()-2)*(_GaussRadauIRplus.Size()-2) // hex inner
                                + 6*sqr(_GaussRadauIRplus.Size()-2)  // faces inside hex
                                , _GaussRadauIRplus.Size()-1),
      GaussRadauIRplus(_GaussRadauIRplus)
    { ; }
    
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TET; }


    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override
    {
      auto [pnts, inds, nump] = GetMicroCellPoints<ET_TET> (order, true);

      Matrix mat(this->ndof, this->ndof);
      Vector rhs(this->ndof);

      Vec<1> fi;
      for (int i = 0; i < nump; i++)
        {
          Vec<3> x = pnts[i];
          IntegrationPoint ip(x(0), x(1), x(2), 0);
          MappedIntegrationPoint<3,3> mip(ip, trafo);
          func.Evaluate (mip, fi);
          rhs(i) = fi(0);
          this->CalcShape(ip, mat.Row(i));
        }

      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }      

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      shape.Range(ndof) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = Trans(Inv(F2*F));

      int ndGR = GaussRadauIRplus.Size();
      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      LagrangePolynomials(xi(0), GaussRadauIRplus, polxi);
      LagrangePolynomials(xi(1), GaussRadauIRplus, poleta);
      LagrangePolynomials(xi(2), GaussRadauIRplus, polzeta);

      auto assign =  [&](int nr, IVec<3> i)
        {
          shape(nr) = polxi[i[0]]*poleta[i[1]]*polzeta[i[2]];
        };


      for (int i = 0; i < 4; i++)
        {
          IVec<3> ind = { 0, 0, 0 };
          if (i == maxlam)          
            assign(i, ind);
        }

      int ii = 4;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);

          // shape at edge mid-point
          for (int j = 0; j < 2; j++)
            if (e[j] == maxlam)
              {
                IVec<3> ind = { 0, 0, 0 };
                int dirv2 = vdir[e[1-j]];
                ind[dirv2] = ndGR-1;
                assign(ii, ind);
              }
          ii++;

          // inner shapes on half-edges
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndGR-2; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR-1; l++)
                    for (int k = 1; k < ndGR-1; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind);
                      }
                }
              ii += sqr(ndGR-2);
            }

          // edge in face between f[j] and f[j+1]
          
          for (int j = 0; j < 3; j++)
            {
              int v0 = f[j];
              int v1 = f[(j+1)%3];
              int v2 = f[(j+2)%3];

              if (v0 == maxlam || v1 == maxlam)
                {
                  int ve2 = (v0 == maxlam) ? v1 : v0;
                  int vop = v2;
                  
                  IVec<3> ind = { 0, 0, 0 };
                  int dirvop = vdir[vop];
                  int dirve2 = vdir[ve2];
                  for (int k = 1; k < ndGR-1; k++)
                    {
                      ind[dirvop] = k;
                      ind[dirve2] = ndGR-1;
                      assign(ii+k-1, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      int iiT = ii+maxlam*sqr(ndGR-2)*(ndGR-2);
      for (int i = 1; i < ndGR-1; i++)
        for (int j = 1; j < ndGR-1; j++)
          for (int k = 1; k < ndGR-1; k++)
            assign(iiT++, { i, j, k } );
      ii += 4*sqr(ndGR-2)*(ndGR-2);

      // faces inside tet
      for (int v0 = 0; v0 < 3; v0++)
        for (int v1 = v0+1; v1 < 4; v1++)
          {
            if (maxlam == v0 || maxlam == v1)
              {
                int vop = (maxlam==v0) ? v1 : v0;
                int v2 = 0;
                if (v2 == v0) v2++;
                if (v2 == v1) v2++;
                int v3 = 6-v0-v1-v2;

                IVec<3> ind = { 0, 0, 0 };
                int dirvop = vdir[vop];
                int dirv2 = vdir[v2];
                int dirv3 = vdir[v3];
                int kk = ii;
                for (int k = 1; k < ndGR-1; k++)
                  for (int l = 1; l < ndGR-1; l++, kk++)
                    {
                      ind[dirvop] = ndGR-1;
                      ind[dirv2] = k;
                      ind[dirv3] = l;
                      assign(kk, ind);
                    }
              }
            
            ii += sqr(ndGR-2);
          }
    }

    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> _dshape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      auto dshape = _dshape.AddSize (ndof, 3);
      dshape = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = Trans(Inv(F2*F));

      int ndGR = GaussRadauIRplus.Size();
      
      ArrayMem<AutoDiff<1>, 20> polxi(ndGR), poleta(ndGR), polzeta(ndGR);
      LagrangePolynomials(AutoDiff<1>(xi(0),0), GaussRadauIRplus, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), GaussRadauIRplus, poleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), GaussRadauIRplus, polzeta);

      auto assign =  [&](int nr, IVec<3> i)
        {
          Vec<3> grad (polxi[i[0]].DValue(0)*poleta[i[1]].Value()*polzeta[i[2]].Value(),
                       polxi[i[0]].Value()*poleta[i[1]].DValue(0)*polzeta[i[2]].Value(), 
                       polxi[i[0]].Value()*poleta[i[1]].Value()*polzeta[i[2]].DValue(0));
          dshape.Row(nr) = trafo * grad;
        };

      for (int i = 0; i < 4; i++)
        {
          IVec<3> ind = { 0, 0, 0 };
          if (i == maxlam)          
            assign(i, ind);
        }

      int ii = 4;
      for (int i = 0; i < 6; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);

          // shape at edge mid-point
          for (int j = 0; j < 2; j++)
            if (e[j] == maxlam)
              {
                IVec<3> ind = { 0, 0, 0 };
                int dirv2 = vdir[e[1-j]];
                ind[dirv2] = ndGR-1;
                assign(ii, ind);
              }
          ii++;

          // inner shapes on half-edges
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<3> ind = { 0, 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndGR-2; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      for (int i = 0; i < 4; i++)
        {
          IVec<4> f = GetVertexOrientedFace(i);
          for (int j = 0; j < 3; j++)
            {
              if (f[j] == maxlam)
                {
                  int v1 = f[(j+1)%3];
                  int v2 = f[(j+2)%3];

                  IVec<3> ind = { 0, 0, 0 };
                  int dirv1 = vdir[v1];
                  int dirv2 = vdir[v2];
                  int kk = ii;
                  for (int l = 1; l < ndGR-1; l++)
                    for (int k = 1; k < ndGR-1; k++)
                      {
                        ind[dirv1] = k;
                        ind[dirv2] = l;
                        assign(kk++, ind);
                      }
                }
              ii += sqr(ndGR-2);
            }

          // edge in face between f[j] and f[j+1]
          
          for (int j = 0; j < 3; j++)
            {
              int v0 = f[j];
              int v1 = f[(j+1)%3];
              int v2 = f[(j+2)%3];

              if (v0 == maxlam || v1 == maxlam)
                {
                  int ve2 = (v0 == maxlam) ? v1 : v0;
                  int vop = v2;
                  
                  IVec<3> ind = { 0, 0, 0 };
                  int dirvop = vdir[vop];
                  int dirve2 = vdir[ve2];
                  for (int k = 1; k < ndGR-1; k++)
                    {
                      ind[dirvop] = k;
                      ind[dirve2] = ndGR-1;
                      assign(ii+k-1, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      int iiT = ii+maxlam*sqr(ndGR-2)*(ndGR-2);
      for (int i = 1; i < ndGR-1; i++)
        for (int j = 1; j < ndGR-1; j++)
          for (int k = 1; k < ndGR-1; k++)
            assign(iiT++, { i, j, k } );
      ii += 4*sqr(ndGR-2)*(ndGR-2);

      // faces inside tet
      for (int v0 = 0; v0 < 3; v0++)
        for (int v1 = v0+1; v1 < 4; v1++)
          {
            if (maxlam == v0 || maxlam == v1)
              {
                int vop = (maxlam==v0) ? v1 : v0;
                int v2 = 0;
                if (v2 == v0) v2++;
                if (v2 == v1) v2++;
                int v3 = 6-v0-v1-v2;

                IVec<3> ind = { 0, 0, 0 };
                int dirvop = vdir[vop];
                int dirv2 = vdir[v2];
                int dirv3 = vdir[v3];
                int kk = ii;
                for (int k = 1; k < ndGR-1; k++)
                  for (int l = 1; l < ndGR-1; l++, kk++)
                    {
                      ind[dirvop] = ndGR-1;
                      ind[dirv2] = k;
                      ind[dirv3] = l;
                      assign(kk, ind);
                    }
              }
            
            ii += sqr(ndGR-2);
          }
    }

  };





  class HCurlDualCellPotentialTrig : public ScalarFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIRplus;
    bool include_central;
  public:
    HCurlDualCellPotentialTrig (const IntegrationRule & _GaussRadauIRplus, bool inc_cent)
      : ScalarFiniteElement<2> (3 + 3         // vertex + edge MP
                                + 6*(_GaussRadauIRplus.Size()-2)   // edges
                                + 3*sqr(_GaussRadauIRplus.Size()-2) // face-quad inner
                                + 3*(_GaussRadauIRplus.Size()-2)    // edges inside face
                                + (inc_cent ? 1 : 0)
                                , _GaussRadauIRplus.Size()-1),
      GaussRadauIRplus(_GaussRadauIRplus), include_central(inc_cent)
    { ; }
    
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const override
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.Range(ndof) = 0;

      int minvi = (maxlam+1)%3;
      int maxvi = minvi;
      for (int i = 0; i < 3; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }

      int vdir[3];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[maxvi] = 1;
      
      Vec<2> x(lam[minvi], lam[maxvi]);
      Vec<2> xi = MapTrig2Quad (x);

      Mat<2,2> F = DMapQuad2Trig(xi);
      Mat<2,2> F2;    // trafo from vertex permutation
      Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[maxvi]-verts[maxlam];
      
      Mat<2> trafo = Trans(Inv(F2*F));
      int ndGR = GaussRadauIRplus.Size();
      
      ArrayMem<double, 20> polxi(ndGR), poleta(ndGR);
      LagrangePolynomials(xi(0), GaussRadauIRplus, polxi);
      LagrangePolynomials(xi(1), GaussRadauIRplus, poleta);

      auto assign =  [&](int nr, IVec<2> i)
        {
          shape(nr) = polxi[i[0]]*poleta[i[1]];
        };

      for (int i = 0; i < 3; i++)
        {
          IVec<2> ind = { 0, 0 };
          if (i == maxlam)          
            assign(i, ind);
        }

      int ii = 3;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);

          // shape at edge mid-point
          for (int j = 0; j < 2; j++)
            if (e[j] == maxlam)
              {
                IVec<2> ind = { 0, 0 };
                int dirv2 = vdir[e[1-j]];
                ind[dirv2] = ndGR-1;
                assign(ii, ind);
              }
          ii++;

          // inner shapes on half-edges
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<2> ind = { 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndGR-2; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      {
        IVec<4> f = GetVertexOrientedFace(0);
        for (int j = 0; j < 3; j++)
          {
            if (f[j] == maxlam)
              {
                int v1 = f[(j+1)%3];
                int v2 = f[(j+2)%3];
                
                IVec<2> ind = { 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 1; l < ndGR-1; l++)
                  for (int k = 1; k < ndGR-1; k++)
                    {
                      ind[dirv1] = k;
                      ind[dirv2] = l;
                      assign(kk++, ind);
                    }
              }
            ii += sqr(ndGR-2);
          }

        // edge in face between f[j] and f[j+1]
        
        for (int j = 0; j < 3; j++)
          {
            int v0 = f[j];
            int v1 = f[(j+1)%3];
            int v2 = f[(j+2)%3];
            
            if (v0 == maxlam || v1 == maxlam)
              {
                int ve2 = (v0 == maxlam) ? v1 : v0;
                int vop = v2;
                
                IVec<2> ind = { 0, 0 };
                int dirvop = vdir[vop];
                int dirve2 = vdir[ve2];
                for (int k = 1; k < ndGR-1; k++)
                  {
                    ind[dirvop] = k;
                    ind[dirve2] = ndGR-1;
                    assign(ii+k-1, ind);
                  }
              }
            ii += ndGR-2;
          }

        if (include_central)
          assign(ii++, IVec<2> (ndGR-1, ndGR-1));
      }
    }

    
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> _dshape) const override
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      auto dshape = _dshape.AddSize (ndof, 2);
      dshape = 0;

      int minvi = (maxlam+1)%3;
      int maxvi = minvi;
      for (int i = 0; i < 3; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }

      int vdir[3];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[maxvi] = 1;
      
      Vec<2> x(lam[minvi], lam[maxvi]);
      Vec<2> xi = MapTrig2Quad (x);
      Mat<2,2> F = DMapQuad2Trig(xi);
      
      Mat<2,2> F2;    // trafo from vertex permutation
      Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[maxvi]-verts[maxlam];
      
      Mat<2> trafo = Trans(Inv(F2*F));

      int ndGR = GaussRadauIRplus.Size();
      
      ArrayMem<AutoDiff<1>, 20> polxi(ndGR), poleta(ndGR);
      LagrangePolynomials(AutoDiff<1>(xi(0),0), GaussRadauIRplus, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), GaussRadauIRplus, poleta);

      auto assign =  [&](int nr, IVec<2> i)
        {
          Vec<2> grad (polxi[i[0]].DValue(0)*poleta[i[1]].Value(), 
                       polxi[i[0]].Value()*poleta[i[1]].DValue(0));
          dshape.Row(nr) = trafo * grad;
        };

      for (int i = 0; i < 3; i++)
        {
          IVec<2> ind = { 0, 0 };
          if (i == maxlam)          
            assign(i, ind);
        }

      int ii = 3;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          // shape at edge mid-point
          for (int j = 0; j < 2; j++)
            if (e[j] == maxlam)
              {
                IVec<2> ind = { 0, 0 };
                int dirv2 = vdir[e[1-j]];
                ind[dirv2] = ndGR-1;
                assign(ii, ind);
              }
          ii++;
          // inner shapes on half-edges
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  IVec<2> ind = { 0, 0 };
                  int dirv2 = vdir[e[1-j]];
                  for (int k = 0; k < ndGR-2; k++)
                    {
                      ind[dirv2] = k+1;
                      assign(ii+k, ind);
                    }
                }
              ii += ndGR-2;
            }
        }

      {
        IVec<4> f = GetVertexOrientedFace(0);
        for (int j = 0; j < 3; j++)
          {
            if (f[j] == maxlam)
              {
                int v1 = f[(j+1)%3];
                int v2 = f[(j+2)%3];
                
                IVec<2> ind = { 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 1; l < ndGR-1; l++)
                  for (int k = 1; k < ndGR-1; k++)
                    {
                      ind[dirv1] = k;
                      ind[dirv2] = l;
                      assign(kk++, ind);
                    }
              }
            ii += sqr(ndGR-2);
          }

        // edge in face between f[j] and f[j+1]
        
        for (int j = 0; j < 3; j++)
          {
            int v0 = f[j];
            int v1 = f[(j+1)%3];
            int v2 = f[(j+2)%3];
            
            if (v0 == maxlam || v1 == maxlam)
              {
                int ve2 = (v0 == maxlam) ? v1 : v0;
                int vop = v2;
                
                IVec<2> ind = { 0, 0 };
                int dirvop = vdir[vop];
                int dirve2 = vdir[ve2];
                for (int k = 1; k < ndGR-1; k++)
                  {
                    ind[dirvop] = k;
                    ind[dirve2] = ndGR-1;
                    assign(ii+k-1, ind);
                  }
              }
            ii += ndGR-2;
          }
      }
      if (include_central)
        assign(ii++, IVec<2> (ndGR-1, ndGR-1));              
    }
    
    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override
    {
      auto [pnts, inds, nump] = GetMicroCellPoints<ET_TRIG> (order, true);
      Matrix mat(this->ndof, this->ndof);
      Vector rhs(this->ndof);
      Vec<1> fi;
      for (int i = 0; i < this->ndof; i++)
        {
          Vec<2> x = pnts[i];
          IntegrationPoint ip(x(0), x(1), 0, 0);
          MappedIntegrationPoint<2,2> mip(ip, trafo);
          func.Evaluate (mip, fi);
          rhs(i) = fi(0);
          this->CalcShape(ip, mat.Row(i));
        }

      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }      
    
  };
    








  
  

  HCurlDualCells::
  HCurlDualCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {
    collocated = !flags.GetDefineFlagX("collocated").IsFalse();
    // cout << "HCurlDual ctor, collocated = " << collocated << endl;

    if (ma->GetDimension()==2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<2>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShape<2>>> ());
      }
    else
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShape<3>>> ());
        additional_evaluators.Set ("Piolashape", make_shared<T_DifferentialOperator<DiffOpPiolaShape<3>>> ());
        additional_evaluators.Set ("Hodge", make_shared<T_DifferentialOperator<DiffOpHodge<3>>> ());        
      }
      
    
    
    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      {
        GaussRadauIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        tangentialIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      }
  }


  void HCurlDualCells :: Update()
  {
    FESpace::Update();
    size_t ndof = 0;

    int ndt = tangentialIR.Size();
    int ndGR = GaussRadauIR.Size();

    first_edge_dofs.SetSize(ma->GetNEdges()+1);
    first_edge_dofs[0] = ndof;
    for (auto i : Range(ma->GetNEdges()))
      first_edge_dofs[i+1] = ndof += 2*ndt;

    first_face_dofs.SetSize(ma->GetNFaces()+1);
    first_face_dofs[0] = ndof;
    for (auto i : Range(ma->GetNFaces()))
      first_face_dofs[i+1] = ndof += 3*2*(ndGR-1)*ndt;

    if (ma->GetDimension() == 3)
      {
        first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
        first_cell_dofs[0] = ndof;
        for (auto i : Range(ma->GetNE(VOL)))
          first_cell_dofs[i+1] = ndof += 4*3*(ndGR-1)*(ndGR-1)*ndt;
      }
    SetNDof(ndof);
  }
    
  void HCurlDualCells ::
  GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    auto el = ma->GetElement(ei);
    for (auto ed : el.Edges())
      dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
    for (auto fa : el.Faces())
      dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
    if (ma->GetDimension()==3 && ei.VB() == VOL)
      dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
    //cout << dnums << endl;
  }


  FiniteElement & HCurlDualCells ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) HCurlDualCellTet(GaussRadauIR, GaussRadauIR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      case ET_TRIG:
        {
          auto trig = new (alloc) HCurlDualCellTrig(GaussRadauIR, GaussRadauIR);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      default:
        throw Exception("element not implemented");
      }
  }



  shared_ptr<BaseMatrix> HCurlDualCells ::
  GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
  {
    if (!collocated)
      throw Exception("only for collocated spaces");

    static Timer t("GetMassOperator"); RegionTimer reg(t);
    static Timer tint("integrate");
    static Timer tsp("MakeSparse");

      
    HeapReset hr(lh);

    int dim = ma->GetDimension();
    switch (dim)
      {
      case 2:
        {
          auto irs = this->GetIntegrationRules();
          IntegrationRule ir = std::move(irs[ET_TRIG]);
          auto & felref = dynamic_cast<const HCurlDualCellTrig&> (GetFE(ElementId(VOL,0), lh));
          Matrix shapes(felref.GetNDof(), 2*ir.Size());
          Matrix shapes_trans(2*ir.Size(), felref.GetNDof());
      
          Array<int> rowind, colind;
          Array<double> values;
      
          Array<short> classnr(ma->GetNE(VOL));
          ma->IterateElements
            (VOL, lh, [&] (auto el, LocalHeap & llh)
             {
               classnr[el.Nr()] = 
                 SwitchET<ET_SEGM, ET_TRIG,ET_TET>
                 (el.GetType(),
                  [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
             });

          TableCreator<size_t> creator;
          for ( ; !creator.Done(); creator++)
            for (auto i : Range(classnr))
              creator.Add (classnr[i], i);
          Table<size_t> table = creator.MoveTable();

      
          for (auto elclass_inds : table)
            {
              if (elclass_inds.Size() == 0) continue;

              HeapReset hr(lh);
          

              size_t nr = elclass_inds[0];
              auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
              auto & felref = dynamic_cast<const HCurlDualCellTrig&> (GetFE(ElementId(VOL,nr), lh));
          
              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Cols(2*i, 2*i+2));
              for (int i = 0; i < shapes.Height(); i++)
                for (int j = 0; j < shapes.Width(); j++)
                  if (fabs(shapes(i,j)) < 1e-8)
                    shapes(i,j) = 0;
              SmallSparseMatrix spshapes(shapes, 1e-8);
              shapes_trans = Trans(shapes);
          
              tint.Start();

              MyMutex COOmutex;
              ParallelForRange (elclass_inds.Size(), [&](auto myrange)
                                {
                                  auto mylh = lh.Split();
                                  Array<int> myrowind, mycolind;
                                  Array<double> myvalues;
                                  Array<DofId> dofs(felref.GetNDof());
                                  Matrix rhoi_shapes_trans(2*ir.Size(), felref.GetNDof());
                                  Matrix<> elmat(felref.GetNDof());
                                  Matrix rhomat(ir.Size(), 4);
                                  Matrix rhoscal(ir.Size(), 1);
              
                                  for (auto i : myrange)
                                    {
                                      HeapReset hr(mylh);
                                      auto nr = elclass_inds[i];

                                      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                                      MappedIntegrationRule<2,2> mir(ir, trafo, mylh);
                                      GetDofNrs(ElementId(VOL,nr), dofs);

                  
                                      if (rho)
                                        {
                                          if (rho->Dimension() == 1)
                                            rho->Evaluate(mir, rhoscal);
                                          else
                                            rho->Evaluate(mir, rhomat);
                                        }
                  
                                      for (size_t i = 0; i < mir.Size(); i++)
                                        {
                                          Mat<2,2> rhoi = Id<2>();
                                          if (rho)
                                            {
                                              if (rho->Dimension() == 1)
                                                rhoi *= rhoscal(i, 0);
                                              else
                                                rhoi = rhomat.Row(i).AsMatrix(2,2);
                                            }
                      
                                          Mat<2,2> Finv = mir[i].GetJacobianInverse();
                                          rhoi = Finv * rhoi * Trans(Finv);
                                          rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                                          rhoi_shapes_trans.Rows(2*i, 2*i+2) = rhoi * shapes_trans.Rows(2*i, 2*i+2);
                                        }
                  
                                      spshapes.Mult (rhoi_shapes_trans, elmat);
                  
                                      for (int i = 0; i < dofs.Size(); i++)
                                        for (int j = 0; j < dofs.Size(); j++)
                                          if (fabs(elmat(i,j)) > 1e-10)
                                            {
                                              myrowind.Append(dofs[i]);
                                              mycolind.Append(dofs[j]);
                                              myvalues.Append(elmat(i,j));
                                            }
                                    }

                                  COOmutex.lock();
                                  rowind += myrowind;
                                  colind += mycolind;
                                  values += myvalues;
                                  COOmutex.unlock();
                                });

              tint.Stop();
            }
          tsp.Start();
          auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
          tsp.Stop();
      
          return make_shared<MySuperSparseMatrix> (move(*spmat));
        }
      case 3:
        {
          auto irs = this->GetIntegrationRules();
          IntegrationRule ir = move(irs[ET_TET]);
          auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,0), lh));
          Matrix shapes(felref.GetNDof(), 3*ir.Size());
          Matrix shapes_trans(3*ir.Size(), felref.GetNDof());
      
          Array<int> rowind, colind;
          Array<double> values;
      
          Array<short> classnr(ma->GetNE(VOL));
          ma->IterateElements
            (VOL, lh, [&] (auto el, LocalHeap & llh)
             {
               classnr[el.Nr()] = 
                 SwitchET<ET_SEGM, ET_TRIG,ET_TET>
                 (el.GetType(),
                  [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
             });

          TableCreator<size_t> creator;
          for ( ; !creator.Done(); creator++)
            for (auto i : Range(classnr))
              creator.Add (classnr[i], i);
          Table<size_t> table = creator.MoveTable();

      
          for (auto elclass_inds : table)
            {
              if (elclass_inds.Size() == 0) continue;

              HeapReset hr(lh);
          

              size_t nr = elclass_inds[0];
              auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
              auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,nr), lh));
          
              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
              for (int i = 0; i < shapes.Height(); i++)
                for (int j = 0; j < shapes.Width(); j++)
                  if (fabs(shapes(i,j)) < 1e-8)
                    shapes(i,j) = 0;
              SmallSparseMatrix spshapes(shapes, 1e-8);
              shapes_trans = Trans(shapes);
          
              tint.Start();

              MyMutex COOmutex;
              ParallelForRange (elclass_inds.Size(), [&](auto myrange)
                                {
                                  auto mylh = lh.Split();
                                  Array<int> myrowind, mycolind;
                                  Array<double> myvalues;
                                  Array<DofId> dofs(felref.GetNDof());
                                  Matrix rhoi_shapes_trans(3*ir.Size(), felref.GetNDof());
                                  Matrix<> elmat(felref.GetNDof());
                                  Matrix rhomat(ir.Size(), 9);
                                  Matrix rhoscal(ir.Size(), 1);
              
                                  for (auto i : myrange)
                                    {
                                      HeapReset hr(mylh);
                                      auto nr = elclass_inds[i];

                                      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                                      MappedIntegrationRule<3,3> mir(ir, trafo, mylh);
                                      GetDofNrs(ElementId(VOL,nr), dofs);

                  
                                      if (rho)
                                        {
                                          if (rho->Dimension() == 1)
                                            rho->Evaluate(mir, rhoscal);
                                          else
                                            rho->Evaluate(mir, rhomat);
                                        }
                  
                                      for (size_t i = 0; i < mir.Size(); i++)
                                        {
                                          Mat<3,3> rhoi = Id<3>();
                                          if (rho)
                                            {
                                              if (rho->Dimension() == 1)
                                                rhoi *= rhoscal(i, 0);
                                              else
                                                rhoi = rhomat.Row(i).AsMatrix(3,3);
                                            }
                      
                                          Mat<3,3> Finv = mir[i].GetJacobianInverse();
                                          rhoi = Finv * rhoi * Trans(Finv);
                                          rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                                          rhoi_shapes_trans.Rows(3*i, 3*i+3) = rhoi * shapes_trans.Rows(3*i, 3*i+3);
                                        }
                  
                                      spshapes.Mult (rhoi_shapes_trans, elmat);
                  
                                      for (int i = 0; i < dofs.Size(); i++)
                                        for (int j = 0; j < dofs.Size(); j++)
                                          if (fabs(elmat(i,j)) > 1e-10)
                                            {
                                              myrowind.Append(dofs[i]);
                                              mycolind.Append(dofs[j]);
                                              myvalues.Append(elmat(i,j));
                                            }
                                    }

                                  COOmutex.lock();
                                  rowind += myrowind;
                                  colind += mycolind;
                                  values += myvalues;
                                  COOmutex.unlock();
                                });
              tint.Stop();
            }
          tsp.Start();
          auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
          tsp.Stop();
      
          return make_shared<MySuperSparseMatrix> (move(*spmat));
        }
      default:
        throw Exception("MassOperator not implemented for dim!=2,3");
      }

    /*
      for (size_t nr : Range(ma->GetNE()))
      {
      HeapReset hr(lh);

      auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
      auto & felref = dynamic_cast<const HCurlDualCellTet&> (GetFE(ElementId(VOL,nr), lh));

      for (int i = 0; i < ir.Size(); i++)
      felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
      for (int i = 0; i < shapes.Height(); i++)
      for (int j = 0; j < shapes.Width(); j++)
      if (fabs(shapes(i,j)) < 1e-8)
      shapes(i,j) = 0;
        
        
      MappedIntegrationRule<3,3> mir(ir, trafo, lh);
      GetDofNrs(ElementId(VOL,nr), dofs);
        
      for (size_t i = 0; i < mir.Size(); i++)
      {
      Mat<3,3> rhoi = Id<3>();
      Mat<3,3> Finv = mir[i].GetJacobianInverse();
      rhoi = Finv * rhoi * Trans(Finv);
      rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
      rhoi_shapes.Cols(3*i, 3*i+3) = shapes.Cols(3*i, 3*i+3) * rhoi;
      }
        
      elmat = rhoi_shapes * Trans(shapes);
        
      for (int i = 0; i < dofs.Size(); i++)
      for (int j = 0; j < dofs.Size(); j++)
      if (fabs(elmat(i,j)) > 1e-10)
      {
      rowind.Append(dofs[i]);
      colind.Append(dofs[j]);
      values.Append(elmat(i,j));
      }
      }
    */

  }

  

  
  std::map<ELEMENT_TYPE, IntegrationRule> HCurlDualCells::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    
    IntegrationRule irtet;
    irtet.SetDim(3);
    int i = 0;
    for (int v = 0; v < 4; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          for (auto ipz : GaussRadauIR)
            {
              Vec<3> xi(ipx(0), ipy(0), ipz(0));
              Vec<3> x = MapHex2Tet (xi, v);
              Mat<3,3> jac = DMapHex2Tet(xi, v);
              double w = ipx.Weight()*ipy.Weight()*ipz.Weight()*fabs(Det(jac));
              IntegrationPoint ip(x(0), x(1), x(2), w);
              ip.SetNr(i);
              irtet += ip;
              i++;
            }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtet)
        ip.SetWeight(1.0/24);   // fix consistency for lowest order



    if (GaussRadauIR.Size()==-2)
      {
        // integrate gradients precisely 
        Array<int> classify = { 0, 1, 1, 2, 1, 2, 2, 3,
                                0, 1, 1, 2, 1, 2, 2, 3,
                                0, 1, 1, 2, 1, 2, 2, 3,
                                0, 1, 1, 2, 1, 2, 2, 3 };
                                
        Array<double> xi, wi;
        ComputeGaussRadauRule (order+1, xi, wi);
        IntegrationRule GaussRadauIRplus;
        for (auto i : Range(xi))
          GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        GaussRadauIRplus.Append (IntegrationPoint (1-1e-10, 0, 0, 0));

        HCurlDualCellPotentialTet tetpot(GaussRadauIRplus);
        tetpot.SetVertexNumbers (Array<int> { 1, 2, 3, 4 } );
        
        Matrix dshape(tetpot.GetNDof(), 3);
        Matrix sumdshape(tetpot.GetNDof(), 3);
        sumdshape = 0.0;

        Matrix mat(4,4);
        Vector rhs(4);
        mat = 0.0;
        rhs = 0.0;
        // sum_j w_j fi(x_j) = \int fi
        for (auto nr : Range(irtet))
          {
            auto ip = irtet[nr];
            tetpot.CalcDShape(ip, dshape);

            mat(0, classify[nr]) += dshape(5,0);
            mat(1, classify[nr]) += dshape(46,0); // first inner
            mat(2, classify[nr]) += dshape(55,0); // last one
            mat(3, classify[nr]) += dshape(4,0);
          }
        // cout << "mat = " << endl << mat << endl;
        CalcInverse (mat);
        // cout << "inv = " << endl << mat << endl;
        
        auto irs = ngcomp::GetIntegrationRules(8);
        IntegrationRule irref = move(irs[ET_TET]);

        for (auto ip : irref)
          {
            tetpot.CalcDShape(ip, dshape);
            // sumdshape += ip.Weight() * dshape;
            
            rhs(0) += ip.Weight()*dshape(5,0); // dshape(0,0);
            rhs(1) += ip.Weight()*dshape(46,0); // first inner
            rhs(2) += ip.Weight()*dshape(55,0); // last one
            rhs(3) += ip.Weight()*dshape(4,0);
          }

        /*
        cout << "exact integrals, vert+edge  " << endl
             << sumdshape.Rows(0,22) << endl;
        cout << "faces" << endl
             << sumdshape.Rows(22,22+24) << endl;
        cout << "cells" << endl
             << sumdshape.Rows(46,sumdshape.Height()) << endl;
        */
        /*
        Vec<3> rhs2 = 0.0;
        for (auto ip : irtrig)
          {
            rhs2(0) += ip.Weight() * f0(ip(0), ip(1));
            rhs2(1) += ip.Weight() * f1(ip(0), ip(1));
            rhs2(2) += ip.Weight() * f2(ip(0), ip(1));
          }
        cout << "GR - integrals =  " << rhs2 << endl;        
        */
        Vector sol(4);
        sol = mat * rhs;

        for (auto i : Range(irtet))
          irtet[i].SetWeight(sol(classify[i]));
      }


    
    rules[ET_TET] = move(irtet);


    IntegrationRule irtrig;
    irtrig.SetDim(2);
    for (int v = 0; v < 3; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          {
            Vec<2> xi(ipx(0), ipy(0));
            Vec<2> x = MapQuad2Trig (xi, v);
            Mat<2,2> jac = DMapQuad2Trig (xi, v);
            double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
            irtrig += IntegrationPoint(x(0), x(1), 0, w);
          }

    if (GaussRadauIR.Size()==1)
      for (auto & ip : irtrig)
        ip.SetWeight(1.0/6);   // fix consistency for lowest order

    if (GaussRadauIR.Size()==-2)
      {
        // A .. corner, B .. on edge, C .. in cell
        // Array<int> classify = { 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1 };
        Array<int> classify = { 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2 };
        
        auto f0 = [](double x, double y) { return 1 ; };
        auto f1 = [](double x, double y)
          {
            // return x*x ;
            double lams[3] = { x, y, 1-x-y };
            if (lams[0] > lams[1] && lams[0] > lams[2])
              {
                return x*y;
              }
            return 0.0;
          };
        auto f2 = [](double x, double y)
          {
            // return x*y*(1-x-y) ;
            double lams[3] = { x, y, 1-x-y };
            if (lams[0] > lams[1] && lams[0] > lams[2])
              {
                return x;
              }
            return 0.0;
          };

        Mat<3,3> mat = 0.0;
        Vec<3> rhs = 0.0;
        // sum_j w_j fi(x_j) = \int fi
        for (auto nr : Range(irtrig))
          {
            auto ip = irtrig[nr];
            mat(0, classify[nr]) += f0(ip(0), ip(1));
            mat(1, classify[nr]) += f1(ip(0), ip(1));
            mat(2, classify[nr]) += f2(ip(0), ip(1));
          }

        auto irs = ngcomp::GetIntegrationRules(5);
        IntegrationRule irref = move(irs[ET_TRIG]);
        // IntegrationRule irref(ET_TRIG, 5);
        for (auto ip : irref)
          {
            rhs(0) += ip.Weight() * f0(ip(0), ip(1));
            rhs(1) += ip.Weight() * f1(ip(0), ip(1));
            rhs(2) += ip.Weight() * f2(ip(0), ip(1));
          }

        cout << "exact integrals =  " << rhs << endl;
        Vec<3> rhs2 = 0.0;
        for (auto ip : irtrig)
          {
            rhs2(0) += ip.Weight() * f0(ip(0), ip(1));
            rhs2(1) += ip.Weight() * f1(ip(0), ip(1));
            rhs2(2) += ip.Weight() * f2(ip(0), ip(1));
          }
        cout << "GR - integrals =  " << rhs2 << endl;        
        
        Vec<3> sol = Inv(mat) * rhs;
        cout << "sol = " << sol << endl;
        // cout << "ir = " << irtrig << endl;

        for (auto i : Range(irtrig))
          irtrig[i].SetWeight(sol(classify[i]));
      }



    if (GaussRadauIR.Size()==-2)
      {
        // integrate gradients precisely 
        Array<int> classify = { 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2 };

        Array<double> xi, wi;
        ComputeGaussRadauRule (order+1, xi, wi);
        IntegrationRule GaussRadauIRplus;
        for (auto i : Range(xi))
          GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
        GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));

        HCurlDualCellPotentialTrig trigpot(GaussRadauIRplus, false);
        trigpot.SetVertexNumbers (Array<int> { 1, 2, 3 } );
        
        Vector shape(trigpot.GetNDof());
        Matrix dshape(trigpot.GetNDof(), 2);
        Vector sumshape(trigpot.GetNDof());
        Matrix sumdshape(trigpot.GetNDof(), 2);

        sumshape = 0.0;
        sumdshape = 0.0;
        
        // sum_j w_j fi(x_j) = \int fi
        for (auto nr : Range(irtrig))
          {
            auto ip = irtrig[nr];
            // trigpot.CalcShape(ip, shape);
            trigpot.CalcDShape(ip, dshape);
            /*
            sumshape += 
            mat(0, classify[nr]) += f0(ip(0), ip(1));
            mat(1, classify[nr]) += f1(ip(0), ip(1));
            mat(2, classify[nr]) += f2(ip(0), ip(1));
            */
          }

        auto irs = ngcomp::GetIntegrationRules(5);
        IntegrationRule irref = move(irs[ET_TRIG]);
        // IntegrationRule irref(ET_TRIG, 5);
        for (auto ip : irref)
          {
            trigpot.CalcDShape(ip, dshape);
            sumdshape += ip.Weight() * dshape;
          }
        
        // cout << "exact integrals =  " << endl << sumdshape << endl;
        /*
        Vec<3> rhs2 = 0.0;
        for (auto ip : irtrig)
          {
            rhs2(0) += ip.Weight() * f0(ip(0), ip(1));
            rhs2(1) += ip.Weight() * f1(ip(0), ip(1));
            rhs2(2) += ip.Weight() * f2(ip(0), ip(1));
          }
        cout << "GR - integrals =  " << rhs2 << endl;        
        
        Vec<3> sol = Inv(mat) * rhs;
        cout << "sol = " << sol << endl;
        // cout << "ir = " << irtrig << endl;

        for (auto i : Range(irtrig))
          irtrig[i].SetWeight(sol(classify[i]));
        */
      }


    
    rules[ET_TRIG] = move(irtrig);

    

    
    return rules;
  }





  /* ******************** Potential space **************************** */


  

  
  
  HCurlDualCellsPotential3D ::
  HCurlDualCellsPotential3D (shared_ptr<MeshAccess> ama, const Flags & flags)
    : PotentialFESpace (ama, flags)            
    {
      if (ma->GetDimension() == 2)
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        }
      else
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        }

      if (flags.GetDefineFlag("include_central"))
        include_central = true;
      
      Array<double> xi, wi;
      ComputeGaussRadauRule (order, xi, wi);
      for (auto i : Range(xi))
        GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));
      /*
        cout << "HCurlDual - Potential ctor" << endl;
        cout << "order = " << order << endl;
        cout << "IR = " << GaussRadauIRplus << endl;
      */
    }

  string HCurlDualCellsPotential3D :: GetClassName () const
  { return "hcurldualcellspotential3d"; }
  
  void HCurlDualCellsPotential3D :: Update() 
    {
      FESpace::Update();
      size_t ndof = 0;
      
      ndof += ma->GetNV();

      first_edge_dofs.SetSize(ma->GetNEdges()+1);
      first_edge_dofs[0] = ndof;
      for (auto i : Range(ma->GetNEdges()))
        first_edge_dofs[i+1] = ndof += 2*order-1;

      first_face_dofs.SetSize(ma->GetNFaces()+1);
      first_face_dofs[0] = ndof;
      
      int nf =  3*sqr(order-1) + 3*(order-1);
      if (include_central)
        {
          if (ma->GetDimension()==2)
            nf += 1;
        }
      for (auto i : Range(ma->GetNFaces()))
        first_face_dofs[i+1] = ndof += nf;

      if (ma->GetDimension() == 3)
        {
          first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
          first_cell_dofs[0] = ndof;
          for (auto i : Range(ma->GetNE(VOL)))
            first_cell_dofs[i+1] = ndof += 4*sqr(order-1)*(order-1) + 6*sqr(order-1);
        }
      
      SetNDof(ndof);
    }

    void HCurlDualCellsPotential3D :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const 
    {
      dnums.SetSize0();
      switch (ni.GetType())
        {
        case NT_VERTEX: dnums.Append(ni.GetNr()); break;
        case NT_EDGE:   dnums += Range(first_edge_dofs[ni.GetNr()],
                                       first_edge_dofs[ni.GetNr()+1]); break;
        case NT_FACE:   dnums += Range(first_face_dofs[ni.GetNr()],
                                       first_face_dofs[ni.GetNr()+1]); break;
        case NT_CELL:   dnums += Range(first_cell_dofs[ni.GetNr()],
                                       first_cell_dofs[ni.GetNr()+1]); break;
        default:
          throw Exception("undefined nodetype");
        }      
    }

    
    void HCurlDualCellsPotential3D :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const 
    {
      dnums.SetSize0();
      auto el = ma->GetElement(ei);
      dnums += el.Vertices();
      for (auto ed : el.Edges())
        dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
      for (auto fa : el.Faces())
        dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
      if (ma->GetDimension()==3)
        if (ei.VB() == VOL)
          dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
    }
    
    BitArray HCurlDualCellsPotential3D :: GetInnerDofs() const 
    {
      BitArray idofs(GetNDof());
      idofs.Clear();
      for (auto v : Range(ma->GetNV()))
        idofs.SetBit(v);
      return idofs;
    }
    
    FiniteElement & HCurlDualCellsPotential3D :: GetFE (ElementId ei, Allocator & alloc) const 
    {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
        {
        case ET_TET:
          {
            auto tet = new (alloc) HCurlDualCellPotentialTet(GaussRadauIRplus);
            tet->SetVertexNumbers (ngel.vertices);
            return *tet;
          }
        case ET_TRIG:
          {
            auto trig = new (alloc) HCurlDualCellPotentialTrig(GaussRadauIRplus, include_central);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }

        default:
          throw Exception("element not implemented");
        }
    }


  shared_ptr<FESpace> HCurlDualCells::GetPotentialSpace(bool include_central) const
  {
    Flags flags;
    flags.SetFlag ("order", order+1);
    flags.SetFlag ("dirichlet", this->flags.GetNumListFlag("dirichlet"));
    if (include_central)
      flags.SetFlag("include_central");
    auto pot = make_shared<HCurlDualCellsPotential3D>(ma, flags);
    pot -> Update();
    pot -> FinalizeUpdate();
    return pot;
  }

  // p-version dofs
  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator3D(bool dual) const
  {
    LocalHeap lh(10*1000*1000);
    shared_ptr<BaseMatrix> sum;
    
    Flags hcflags;
    hcflags.SetFlag ("order", order);
    auto fescurl = make_shared<HCurlPrimalCells> (ma, hcflags);
    fescurl->Update();
    fescurl->FinalizeUpdate();
        
    Array<short> classnr(ma->GetNE(VOL));
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        // size_t nr = classnr[elclass_inds[0]];
        ElementId ei(VOL,elclass_inds[0]);

        // here we go ....
        
        auto & fel = dynamic_cast<const HCurlCellFiniteElement<3>&> (GetFE (ei, lh));
        auto & felc = dynamic_cast<const HCurlCellFiniteElement<3>&> (fescurl->GetFE (ei, lh));
        auto & trafo = ma->GetTrafo(ei, lh);

        IntegrationRule irGR2;
        for (auto ip : GaussRadauIR)
          irGR2.Append (IntegrationPoint (1-1e-32-ip(0), 0, 0, 0));
        
        H1PrimalCellTrig trig(order, irGR2);
        int ndoftrig = trig.GetNDof();
        // total dofs = 4 * (order+1) * ndoftrig

        Matrix<> curlmat(felc.GetNDof(), fel.GetNDof());
        Matrix<> mixedmass(felc.GetNDof(), felc.GetNDof());

        Matrix<> shape(fel.GetNDof(), 3);
        Matrix<> shapec(felc.GetNDof(), 3);
        Matrix<> transshape(fel.GetNDof(), 3);
        Matrix<> transshapec(felc.GetNDof(), 3);
        
        Vector<> shapetrig(trig.GetNDof());
        Matrix<> dshapetrig(trig.GetNDof(),2);
        
        curlmat = 0.;
        mixedmass = 0.;

        IntegrationRule irsegm_broken = std::move(ngcomp::GetIntegrationRules(2*order+8)[ET_SEGM]);
        
        IntegrationRule irsegm(ET_SEGM, 2*order+8);
        IntegrationRule irquad(ET_QUAD, 2*order+8);
        L2HighOrderFE<ET_SEGM> segm(order-1);
        Vector<> shapez(segm.GetNDof());
        // Matrix<> dshapez(segm.GetNDof(),1);  

        Matrix<> testshape1(ndoftrig,3);
        Matrix<> transtestshape1(ndoftrig,3);
        Matrix<> testshape(ndoftrig*segm.GetNDof(),3);
        Matrix<> transtestshape(ndoftrig*segm.GetNDof(),3);

        // bool PiolashapesE = false;
        // bool altshapesE = true; // does not matter

        // one of POLYNOMIAL, COVARIANT, PIOLA:
        auto mappingtypeE = POLYNOMIAL; 
        auto mappingtypeH = POLYNOMIAL;
        
        for (int f = 0; f < 4; f++)
          {
            int ii = f*ndoftrig*(order+1);
            IntRange range_face(ii, ii+ndoftrig);
            IntRange range_el(ii+ndoftrig, (f+1)*ndoftrig*(order+1));
            
            IVec<4> fverts = ET_trait<ET_TET>::GetFace(f);
            int vop = 6 - fverts[0]-fverts[1]-fverts[2];
            Vec<3> points2d[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2>(0,0) };
            Vec<3> points3d[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3>(0,0,1), Vec<3>(0,0,0) };
            
            // distributional edge terms for curlmat:
            for (int fe = 0; fe < 3; fe++) // face edges
              {
                Vec<2> p12d = points2d[fe];
                Vec<2> p22d = points2d[(fe+1)%3];
                Vec<3> p13d = points3d[fverts[fe]];
                Vec<3> p23d = points3d[fverts[(fe+1)%3]];
                
                //  int_e shape * tau * shapeface ds
                for (auto ip : irsegm_broken)
                  {
                    Vec<3> tau3d = p23d-p13d;
                    Vec<2> tau2d = p22d-p12d;
                    IntegrationPoint ip3d(p13d+ip(0)*tau3d);
                    IntegrationPoint ip2d(p12d+ip(0)*tau2d);

                    fel.CalcShape(ip3d, shape, mappingtypeE);
                    trig.CalcShape(ip2d, shapetrig); 
                    
                    curlmat.Rows(range_face) += ip.Weight()*shapetrig * Trans(shape * tau3d);
                  }
              }


            
            // int_f shape \cross n dshapeelement  on other faces
            for (int fv = 0; fv < 3; fv++)
              for (int side = 0; side < 2; side++)
                {
                  for (auto ipe : irsegm)
                    for (auto ipz : irsegm)
                      {
                        Vec<2> ipquad = (side==1) ? Vec<2>(ipe(0), 0) : Vec<2> (0, ipe(0));
                        Vec<2> ptrig = MapQuad2Trig (ipquad, fv);
                        Mat<2,2> Dptrig = DMapQuad2Trig (Vec<2> (ipquad(0), ipquad(1)), fv);
                        Vec<3> ptet = MapHex2Tet (Vec<3> (ipquad(0), ipquad(1), ipz(0)), fverts[fv], vop); 
                        Mat<3,3> Dptet = DMapHex2Tet (Vec<3> (ipquad(0), ipquad(1), ipz(0)), fverts[fv], vop);
                        
                        Mat<3,3> Dpprism = { Dptrig(0,0), Dptrig(0,1), 0,
                                             Dptrig(1,0), Dptrig(1,1), 0,
                                             0, 0, 1 };
                        Mat<3,3> Dprism2tet = Dptet*Inv(Dpprism);
                        
                        segm.CalcShape (ipz, shapez);
                        trig.CalcShape (IntegrationPoint(ptrig), shapetrig);
                        fel.CalcShape (IntegrationPoint(ptet), shape, mappingtypeE);
                          
                        testshape=0;
                        Vec<2> tquad = (side==1) ? Vec<2>(1,0) : Vec<2> (0,-1);
                        Vec<2> ttrig = Dptrig * tquad;
                        for (int k = 0, ll=0; k < shapez.Size(); k++)
                          for (int l = 0; l < shapetrig.Size(); l++, ll++)
                            {
                              testshape(ll,0) = shapez(k)*shapetrig(l) * ttrig(0);
                              testshape(ll,1) = shapez(k)*shapetrig(l) * ttrig(1);
                            }
                        
                        transshape = shape * Dprism2tet;

                        curlmat.Rows(range_el) +=
                          ipe.Weight()*ipz.Weight()*  testshape * Trans(transshape);
                      }
                }
            
            for (int fv = 0; fv < 3; fv++) 
              { // now integrate over quad and hex above quad

                for (auto ipquad : irquad)
                  {
                    Vec<2> ptrig = MapQuad2Trig (Vec<2> (ipquad(0), ipquad(1)), fv);
                    Mat<2,2> Dptrig = DMapQuad2Trig (Vec<2> (ipquad(0), ipquad(1)), fv);
                    Vec<3> ptet = MapHex2Tet (Vec<3> (ipquad(0), ipquad(1), 0), fverts[fv], vop); 
                    Mat<3,3> Dptet = DMapHex2Tet (Vec<3> (ipquad(0), ipquad(1), 0), fverts[fv], vop);

                    Mat<3,3> Dpprism = { Dptrig(0,0), Dptrig(0,1), 0,
                                         Dptrig(1,0), Dptrig(1,1), 0,
                                         0, 0, 1 };
                    Mat<3,3> Dprism2tet = Dptet*Inv(Dpprism);

                    fel.CalcShape (IntegrationPoint(ptet), shape, mappingtypeE);
                    felc.CalcShape (IntegrationPoint(ptet), shapec, mappingtypeH);
                    trig.CalcShape (IntegrationPoint(ptrig), shapetrig);
                    *testout << "ptrig = " << endl << ptrig << endl;
                    *testout << "dshapetrig,before = " << endl << dshapetrig << endl;
                    trig.CalcDShape (IntegrationPoint(ptrig), dshapetrig);

                    testshape1 = 0;
                    testshape1.Col(2) = shapetrig;
                    
                    transshapec = shapec * Cof(Dprism2tet);
                    mixedmass.Rows(range_face) +=
                      ipquad.Weight()*Det(Dptrig) * testshape1 * Trans(transshapec);
                    
                    // surface curl:
                    testshape1 = 0;
                    testshape1.Col(0) = dshapetrig.Col(1);
                    testshape1.Col(1) = -dshapetrig.Col(0);

                    transshape = shape * Dprism2tet;
                    curlmat.Rows(range_face) +=
                      ipquad.Weight()*Det(Dptrig) * testshape1 * Trans(transshape);

                    
                    
                    for (auto ipz : irsegm)
                      {
                        Vec<2> ptrig = MapQuad2Trig (Vec<2> (ipquad(0), ipquad(1)), fv);
                        Mat<2,2> Dptrig = DMapQuad2Trig (Vec<2> (ipquad(0), ipquad(1)), fv);
                        Vec<3> ptet = MapHex2Tet (Vec<3> (ipquad(0), ipquad(1), ipz(0)), fverts[fv], vop); 
                        Mat<3,3> Dptet = DMapHex2Tet (Vec<3> (ipquad(0), ipquad(1), ipz(0)), fverts[fv], vop);
                        
                        Mat<3,3> Dpprism = { Dptrig(0,0), Dptrig(0,1), 0,
                                             Dptrig(1,0), Dptrig(1,1), 0,
                                             0, 0, 1 };
                        Mat<3,3> Dprism2tet = Dptet*Inv(Dpprism);
                        
                        segm.CalcShape (ipz, shapez);
                        // shapez *= 1-ipz(0);
                        // segm.CalcDShape (ipz, dshapez);                        
                        trig.CalcShape (IntegrationPoint(ptrig), shapetrig);
                        trig.CalcDShape (IntegrationPoint(ptrig), dshapetrig);
                        fel.CalcShape (IntegrationPoint(ptet), shape, mappingtypeE);
                        felc.CalcShape (IntegrationPoint(ptet), shapec, mappingtypeH);
                        
                        testshape=0;
                        for (int k = 0, ll=0; k < shapez.Size(); k++)
                          for (int l = 0; l < shapetrig.Size(); l++, ll++)
                            testshape(ll,2) = shapez(k)*shapetrig(l);

                        transshapec = shapec * Cof(Dprism2tet);
                        mixedmass.Rows(range_el) +=
                          ipquad.Weight()*ipz.Weight()*Det(Dptrig) *  testshape * Trans(transshapec);

                        
                        testshape=0;
                        for (int k = 0, ll=0; k < shapez.Size(); k++)
                          for (int l = 0; l < shapetrig.Size(); l++, ll++)
                            {
                              testshape(ll,0) = shapez(k)*dshapetrig(l,1);
                              testshape(ll,1) = -shapez(k)*dshapetrig(l,0);
                              testshape(ll,2) = 0;
                            }

                        transshape = shape * Dprism2tet;
                        curlmat.Rows(range_el) +=
                          ipquad.Weight()*ipz.Weight()*Det(Dptrig) *  testshape * Trans(transshape);
                      }
                  }
              }
          }

        // checking for nans:
        // cout << "norm curlmat = " << L2Norm(curlmat) << endl;
        // cout << "norm mixedmat = " << L2Norm(mixedmass) << endl;
        // *testout << "curlmat = " << endl << curlmat << endl;
        // *testout << "mixedmass = " << endl << mixedmass << endl;
        
        CalcInverse (mixedmass);
        // cout << "norm inv mixedmat = " << L2Norm(mixedmass) << endl;        
        Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
        tmp = mixedmass * curlmat;


        Matrix<> finalmat(felc.GetNDof(), fel.GetNDof());
        if (dual)
        {
          // multiply finalmat with lumpedmass
          finalmat = tmp;
          cout << "dual not implemented for 3D" << endl;
        }
        else
          finalmat = tmp;

        Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
          ydofs(elclass_inds.Size(), felc.GetNDof());
        
        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            GetDofNrs(ei, dnumsx);
            fescurl->GetDofNrs(ei, dnumsy);
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }
        
        auto mat = make_shared<ConstantElementByElementMatrix>
          (fescurl->GetNDof(), GetNDof(),
           finalmat, std::move(ydofs), std::move(xdofs));
        
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }
    
    return sum;
  }

  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator3DNano(bool dual) const
  {
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;
        Flags h1flags;
        h1flags.SetFlag ("order", order);
        auto fescurl = make_shared<HCurlPrimalCells> (ma, h1flags);
        fescurl->Update();
        fescurl->FinalizeUpdate();
        
        Array<short> classnr(ma->GetNE(VOL));
        ma->IterateElements
          (VOL, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
    
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();

        IntegrationRule irquad(ET_QUAD, 2*order+4);
        IntegrationRule irsegm(ET_SEGM, 2*order+4);

        // bool altshapesH = false;
        auto mappingtypeB = POLYNOMIAL;
        auto mappingtypeH = POLYNOMIAL;        
        
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            ElementId ei(VOL,elclass_inds[0]);
        
            auto & fel = dynamic_cast<const HCurlCellFiniteElement<3>&> (GetFE (ei, lh));
            auto & felc = dynamic_cast<const HCurlCellFiniteElement<3>&> (fescurl->GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);

        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TET> (order+1);
            auto faces = GetMicroCellFaces<ET_TET> (order+1);
            auto facebnds = GetMicroCellFaceBoundaries<ET_TET> (order+1);

            Matrix<> curlmat(faces.Size(), fel.GetNDof());
            Matrix<> mixedmass(faces.Size(), felc.GetNDof());
            Matrix<> lumpedmass(felc.GetNDof(), felc.GetNDof());

            curlmat = 0.;
            mixedmass = 0.;
            lumpedmass = 0.;
            Matrix shape(fel.GetNDof(), 3);
            Matrix shapecurl(felc.GetNDof(), 3);

            // curlmat ...
            for (auto i : Range(facebnds))
              for (auto fb : facebnds[i])
                {
                  Vec<3> micropts[2] = {
                    pnts[fb[0]],
                    pnts[fb[1]]
                  };
                  Vec<3> tau = micropts[1]-micropts[0];
                  for (auto ip : irsegm)
                    {
                      Vec<3> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), p(2), 0);
                      fel.CalcShape(ipseg, shape);    // same as CalcAltShapes
                      curlmat.Row(i) += ip.Weight() * shape*tau;
                    }
                }
        
            // mixedmass
            for (auto i : Range(faces))
              for (auto f : faces[i])
                {
                  Vec<3> micropts[4] =
                    { pnts[f[0]], pnts[f[1]], pnts[f[2]], pnts[f[3]] };
              
                  for (auto ip : irquad)
                    {
                      AutoDiff<2> adx(ip(0), 0);
                      AutoDiff<2> ady(ip(1), 1);
                      AutoDiff<2> adlams[] =
                        { (1-adx)*(1-ady), adx*(1-ady), adx*ady, (1-adx)*ady };
                  
                      Vec<3, AutoDiff<2>> adp = AutoDiff<2>(0.0);
                      for (int j = 0; j < 4; j++)
                        adp += adlams[j] * micropts[j];

                      Mat<3,2> jac;                  
                      for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 2; l++)
                          jac(k,l) = adp(k).DValue(l);
                  
                      IntegrationPoint iptrig(adp(0).Value(), adp(1).Value(), adp(2).Value(), 0);
                      felc.CalcShape(iptrig, shapecurl, mappingtypeB);                      
                      mixedmass.Row(i) += ip.Weight()* shapecurl * Cross(jac.Col(0), jac.Col(1));
                    }
                }

            // cout << "mixedmass = " << mixedmass << endl;
        
            // integrate
            // lumpedmass ...

            // lumped integration is the same (except order=0 fixup)
            // auto irs = fescurl->GetIntegrationRules();
            // IntegrationRule ir = move(irs[ET_TET]);
            
            auto irs = ngcomp::GetIntegrationRules(2*order+4);
            IntegrationRule ir = move(irs[ET_TET]);            
            
            Matrix shapes(felc.GetNDof(), ir.Size()*3);
            Matrix Piolashapes(felc.GetNDof(), ir.Size()*3);
        
            for (int i = 0; i < ir.Size(); i++)
              {
                felc.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3), mappingtypeH);
                felc.CalcShape (ir[i], Piolashapes.Cols(3*i, 3*i+3), mappingtypeB);
                Piolashapes.Cols(3*i, 3*i+3) *= ir[i].Weight();
              }

            lumpedmass = shapes * Trans (Piolashapes);

            /*
             *testout << "curlmat = " << curlmat << endl;
             *testout << "mixedmass = " << mixedmass << endl;
             *testout << "lumpedmass = " << lumpedmass << endl;
             */
            
            CalcInverse (mixedmass);
            *testout << "inv mixedmass = " << mixedmass << endl;            
            Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
            tmp = mixedmass * curlmat;

            if (dual)
              curlmat = lumpedmass * tmp;
            else
              curlmat = tmp;              

            // *testout << "final mat = " << curlmat << endl;
        
            Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
              ydofs(elclass_inds.Size(), felc.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsx);
                fescurl->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            auto mat = make_shared<ConstantElementByElementMatrix>
              (fescurl->GetNDof(), GetNDof(),
               curlmat, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        return sum;
  }

  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator2D(bool dual) const
  {
    throw Exception("not implemented for 2D");
  }
  
  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator2DNano(bool dual) const
  { 
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;
        Flags h1flags;
        bool printmats = false;
        h1flags.SetFlag ("order", order);
        auto fescurl = make_shared<H1PrimalCells> (ma, h1flags);
        fescurl->Update();
        fescurl->FinalizeUpdate();
        
        Array<short> classnr(ma->GetNE(VOL));
        ma->IterateElements
          (VOL, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
    
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();

        IntegrationRule irquad(ET_QUAD, 2*order+4);
        IntegrationRule irsegm(ET_SEGM, 2*order+4);
    
        //MAPPINGTYPE mappingE = COVARIANT; // POLYNOMIAL doesn't make any difference
        MAPPINGTYPE mappingE = POLYNOMIAL;
        //MAPPINGTYPEH1 mappingB = L2; // L2
        MAPPINGTYPEH1 mappingB = L2POLYNOMIAL; // 
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(VOL,elclass_inds[0]);
        
            auto & fel = dynamic_cast<const HCurlCellFiniteElement<2>&> (GetFE (ei, lh));
            auto & felc = dynamic_cast<const H1PrimalCellTrig&> (fescurl->GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);

            //cout << "fel.ndof = " << fel.GetNDof() << endl;
            //cout << "felc.ndof = " << felc.GetNDof() << endl;
        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TRIG> (order+1);
            auto faces = GetMicroCellFaces<ET_TRIG> (order+1);
            auto facebnds = GetMicroCellFaceBoundaries<ET_TRIG> (order+1);

            Matrix<> curlmat(faces.Size(), fel.GetNDof());
            Matrix<> mixedmass(faces.Size(), felc.GetNDof());
            Matrix<> lumpedmass(felc.GetNDof(), felc.GetNDof());

            curlmat = 0.;
            mixedmass = 0.;
            lumpedmass = 0.;
            Matrix shape(fel.GetNDof(), 2);
            Vector shapecurl(felc.GetNDof());



            // curlmat ...
            for (auto i : Range(facebnds))
              for (auto fb : facebnds[i])
                {
                  Vec<2> micropts[2] = {
                    pnts[fb[0]],
                    pnts[fb[1]]
                  };
                  Vec<2> tau = micropts[1]-micropts[0];
                  for (auto ip : irsegm)
                    {
                      Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), 0, 0);
                      fel.CalcShape(ipseg, shape,mappingE);
                      //fel.CalcAltShape(ipseg, shape);  
                      curlmat.Row(i) += ip.Weight() * shape*tau;
                    }
                }
            if (printmats)
            cout << "curlmat = " << endl << curlmat << endl;

        
            // mixedmass
            for (auto i : Range(faces))
              for (auto f : faces[i])
                {
                  Vec<2> micropts[4] =
                    { pnts[f[0]], pnts[f[1]], pnts[f[2]], pnts[f[3]] };
              
                  for (auto ip : irquad)
                    {
                      AutoDiff<2> adx(ip(0), 0);
                      AutoDiff<2> ady(ip(1), 1);
                      AutoDiff<2> adlams[] =
                        { (1-adx)*(1-ady), adx*(1-ady), adx*ady, (1-adx)*ady };
                  
                      Vec<2, AutoDiff<2>> adp = AutoDiff<2>(0.0);
                      for (int j = 0; j < 4; j++)
                        adp += adlams[j] * micropts[j];

                      Mat<2,2> jac;                  
                      for (int k = 0; k < 2; k++)
                        for (int l = 0; l < 2; l++)
                          jac(k,l) = adp(k).DValue(l);
                  
                      IntegrationPoint iptrig(adp(0).Value(), adp(1).Value(), 0, 0);
                      felc.CalcShape(iptrig, shapecurl,mappingB);
                  
                      mixedmass.Row(i) += ip.Weight()*Det(jac) * shapecurl;
                    }
                }

            if (printmats)
            cout << "mixedmass = " << mixedmass << endl;
        
            // integrate
            // lumpedmass ...
            auto irs = fescurl->GetIntegrationRules();
            //auto irs = ngcomp::GetIntegrationRules(2*order+6);

            IntegrationRule ir = move(irs[ET_TRIG]);

            Matrix shapes(felc.GetNDof(), ir.Size());
            Matrix shapesB(felc.GetNDof(), ir.Size());
        
            for (int i = 0; i < ir.Size(); i++)
            {
              felc.CalcShape (ir[i], shapes.Col(i));
              felc.CalcShape (ir[i], shapesB.Col(i),mappingB);
              shapesB.Col(i) *= ir[i].Weight();
            }

            // cout << "shapes = " << endl << shapes << endl;
                  
            lumpedmass = shapes * Trans(shapesB);
            if (printmats)
            cout << "lumpedmass = " << endl << lumpedmass << endl;
            // cout << "mixedmass = " << endl << mixedmass << endl;

            CalcInverse (mixedmass);
            // cout << "mixedmass^-1 = " << mixedmass << endl;
            Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
            tmp = mixedmass * curlmat;

            if (dual)
              curlmat = lumpedmass * tmp;
            else
              curlmat = tmp;              

            // cout << "final mat = " << curlmat << endl;
        
            Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
              ydofs(elclass_inds.Size(), felc.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsx);
                fescurl->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            /*
              cout << "elmat = " << curlmat << endl;
              cout << "xdofs = " << xdofs << endl;
              cout << "ydofs = " << ydofs << endl;
            */
        
            auto mat = make_shared<ConstantElementByElementMatrix>
              (fescurl->GetNDof(), GetNDof(),
               curlmat, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        return sum;
  }


  shared_ptr<BaseMatrix> HCurlDualCells::GetCurlOperator2DKronecker(bool dual) const
  {
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;
        Flags h1flags;
        h1flags.SetFlag ("order", order);
        auto fescurl = make_shared<H1PrimalCells> (ma, h1flags);
        fescurl->Update();
        fescurl->FinalizeUpdate();
        
        Array<short> classnr(ma->GetNE(VOL));
        ma->IterateElements
          (VOL, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
    
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();

        IntegrationRule irquad(ET_QUAD, 2*order+4);
        IntegrationRule irsegm(ET_SEGM, 2*order+4);
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(VOL,elclass_inds[0]);
        
            auto & fel = dynamic_cast<const HCurlFiniteElement<2>&> (GetFE (ei, lh));
            auto & felc = dynamic_cast<const ScalarFiniteElement<2>&> (fescurl->GetFE (ei, lh));

            auto & trafo = ma->GetTrafo(ei, lh);

            Matrix curlmat(felc.GetNDof(), fel.GetNDof());
            curlmat = 0.0;
            // with boundary point 0 or 1 included
            Array<double> xi, wi;
            IntegrationRule GRruleLeft, GRruleRight;
            ComputeGaussRadauRule (order+1, xi, wi);
            for (auto i : Range(xi))
              {
                GRruleLeft.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
                GRruleRight.Append (IntegrationPoint (1-xi[i], 0, 0, wi[i]));                
              }
            
            H1PrimalCellSegm felh1d(order, GRruleRight);

            Matrix shapeE(fel.GetNDof(), 2);
            Matrix grad1d(felh1d.GetNDof(), 1);
            Vector shape1d(felh1d.GetNDof());
            Vector shapec(felc.GetNDof());            
            Matrix grad1delmat(felh1d.GetNDof(), fel.GetNDof());
            Matrix mass1delmat(felh1d.GetNDof());

            for (int frustum = 0; frustum < 3; frustum++)
              {
                auto e2v = ElementTopology::GetEdges(ET_TRIG)[frustum];
                int vop = 3-e2v[0]-e2v[1];
                for (auto ipeta : GRruleLeft)
                  {
                    grad1delmat = 0;
                    mass1delmat = 0;
                    for (int side = 0; side < 2; side++)
                      {
                        for (auto ipxi : GRruleLeft)
                          {
                            double xi = ipxi(0);
                            double eta = ipeta(0);
                            double x1d = side ? 0.5*xi : 1-0.5*xi;
                            Vec<2> x2d = MapQuad2Trig(Vec<2>(xi,eta), e2v[side], vop);
                            Mat<2,2> jac = DMapQuad2Trig (Vec<2>(xi,eta), e2v[side], vop);

                            IntegrationPoint ip2d(x2d);
                            IntegrationPoint ip1d(x1d);
                            fel.CalcShape (ip2d, shapeE);
                            felh1d.CalcDShape (ip1d, grad1d);
                            // if (side==1) grad1d*=-1;
                            felh1d.CalcShape (ip1d, shape1d);                            
                            
                            grad1delmat += ipxi.Weight() * grad1d * Trans( (shapeE*jac.Col(1) ) );
                            mass1delmat += ipxi.Weight() * shape1d*Trans(shape1d);
                          }
                        
                        double xi = 0;
                        double eta = ipeta(0);
                        double x1d = side ? 0.5*xi : 1-0.5*xi;
                        Vec<2> x2d = MapQuad2Trig(Vec<2>(xi,eta), e2v[side], vop);
                        Mat<2,2> jac = DMapQuad2Trig (Vec<2>(xi,eta), e2v[side], vop);

                        IntegrationPoint ip2d(x2d);
                        IntegrationPoint ip1d(x1d);

                        double normal = (side==0) ? -2 : 2;

                        fel.CalcShape (ip2d, shapeE);
                        felh1d.CalcShape (ip1d, shape1d);                            
                        
                        // cout << "shape1d = " << shape1d << endl;
                        // cout << "shapeEeta = " << shapeE*jac.Col(1) << endl;                        
                        grad1delmat += normal * shape1d * Trans( (shapeE*jac.Col(1) ) );
                      }
                    CalcInverse(mass1delmat);
                    Matrix E2HJ1d = mass1delmat * grad1delmat;
                    // cout << "diff1d = " << endl << Truncate(E2HJ1d,1e-9) << endl;
                    

                    for (int side = 0; side < 2; side++)
                      {
                        for (auto ipxi : GRruleRight)
                          {
                            double xi = ipxi(0);
                            double eta = ipeta(0);
                            double x1d = side ? 0.5*xi : 1-0.5*xi;
                            Vec<2> x2d = MapQuad2Trig(Vec<2>(xi,eta), e2v[side], vop);
                            Mat<2,2> jac = DMapQuad2Trig (Vec<2>(xi,eta), e2v[side], vop);

                            IntegrationPoint ip2d(x2d);
                            IntegrationPoint ip1d(x1d);

                            felh1d.CalcShape (ip1d, shape1d);
                            felc.CalcShape(ip2d, shapec);

                            curlmat += ipeta.Weight()*ipxi.Weight() *  shapec * Trans(Vector (Trans(E2HJ1d)*shape1d));
                          }
                      }
                  }
              }

            // cout << "curlmat = " << curlmat << endl;
            
#ifdef OLD            
            
            //cout << "fel.ndof = " << fel.GetNDof() << endl;
            //cout << "felc.ndof = " << felc.GetNDof() << endl;
        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TRIG> (order+1);
            auto faces = GetMicroCellFaces<ET_TRIG> (order+1);
            auto facebnds = GetMicroCellFaceBoundaries<ET_TRIG> (order+1);

            Matrix<> curlmat(faces.Size(), fel.GetNDof());
            Matrix<> mixedmass(faces.Size(), felc.GetNDof());
            Matrix<> lumpedmass(felc.GetNDof(), felc.GetNDof());

            curlmat = 0.;
            mixedmass = 0.;
            lumpedmass = 0.;
            Matrix shape(fel.GetNDof(), 2);
            Vector shapecurl(felc.GetNDof());

            // curlmat ...
            for (auto i : Range(facebnds))
              for (auto fb : facebnds[i])
                {
                  Vec<2> micropts[2] = {
                    pnts[fb[0]],
                    pnts[fb[1]]
                  };
                  Vec<2> tau = micropts[1]-micropts[0];
                  for (auto ip : irsegm)
                    {
                      Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), 0, 0);
                      fel.CalcShape(ipseg, shape);
                      curlmat.Row(i) += ip.Weight() * shape*tau;
                    }
                }
            // cout << "curlmat = " << endl << curlmat << endl;

        
            // mixedmass
            for (auto i : Range(faces))
              for (auto f : faces[i])
                {
                  Vec<2> micropts[4] =
                    { pnts[f[0]], pnts[f[1]], pnts[f[2]], pnts[f[3]] };
              
                  for (auto ip : irquad)
                    {
                      AutoDiff<2> adx(ip(0), 0);
                      AutoDiff<2> ady(ip(1), 1);
                      AutoDiff<2> adlams[] =
                        { (1-adx)*(1-ady), adx*(1-ady), adx*ady, (1-adx)*ady };
                  
                      Vec<2, AutoDiff<2>> adp = AutoDiff<2>(0.0);
                      for (int j = 0; j < 4; j++)
                        adp += adlams[j] * micropts[j];

                      Mat<2,2> jac;                  
                      for (int k = 0; k < 2; k++)
                        for (int l = 0; l < 2; l++)
                          jac(k,l) = adp(k).DValue(l);
                  
                      IntegrationPoint iptrig(adp(0).Value(), adp(1).Value(), 0, 0);
                      felc.CalcShape(iptrig, shapecurl);
                  
                      mixedmass.Row(i) += ip.Weight()*Det(jac) * shapecurl;
                    }
                }

            // cout << "mixedmass = " << mixedmass << endl;
        
            // integrate
            // lumpedmass ...
            auto irs = fescurl->GetIntegrationRules();
            // auto irs = ngcomp::GetIntegrationRules(2*order+6);
            IntegrationRule ir = move(irs[ET_TRIG]);
            Matrix shapes(felc.GetNDof(), ir.Size());
            Matrix rhoi_shapes_trans(ir.Size(), felc.GetNDof());
        
            for (int i = 0; i < ir.Size(); i++)
              felc.CalcShape (ir[i], shapes.Col(i));

            // cout << "shapes = " << endl << shapes << endl;
        
            MappedIntegrationRule<2,2> mir(ir, trafo, lh);
            for (size_t i = 0; i < mir.Size(); i++)
              rhoi_shapes_trans.Row(i) = mir[i].IP().Weight() * shapes.Col(i);
                  
            lumpedmass = shapes * rhoi_shapes_trans;
            // cout << "lumpedmass = " << endl << lumpedmass << endl;
        
            // cout << "mixedmass = " << endl << mixedmass << endl;

            CalcInverse (mixedmass);
            // cout << "mixedmass^-1 = " << mixedmass << endl;
            Matrix<> tmp(felc.GetNDof(), fel.GetNDof());
            tmp = mixedmass * curlmat;

            if (dual)
              curlmat = lumpedmass * tmp;
            else
              curlmat = tmp;              
#endif
            
            // cout << "final mat = " << curlmat << endl;
        
            Table<DofId> xdofs(elclass_inds.Size(), fel.GetNDof()),
              ydofs(elclass_inds.Size(), felc.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsx);
                fescurl->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            /*
              cout << "elmat = " << curlmat << endl;
              cout << "xdofs = " << xdofs << endl;
              cout << "ydofs = " << ydofs << endl;
            */
        
            auto mat = make_shared<ConstantElementByElementMatrix>
              (fescurl->GetNDof(), GetNDof(),
               curlmat, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        return sum;
  }



  
  shared_ptr<BaseMatrix> HCurlDualCells::
  GetCurlOperator(bool dual, bool nanocells, bool Kronecker) const
  {
    if (ma -> GetDimension() == 2)
      {
        if (Kronecker)
          return GetCurlOperator2DKronecker(dual);
        else if (nanocells)
          return GetCurlOperator2DNano(dual);
        else
          return GetCurlOperator2D(dual);
    }
    else // dim == 3
    {
      if (nanocells)
          return GetCurlOperator3DNano(dual);
      else
          return GetCurlOperator3D(dual);
    }
  }

  

  shared_ptr<BaseMatrix> HCurlDualCells::GetGradientOperator3D () const
  {
    LocalHeap lh(10*1000*1000);
    shared_ptr<BaseMatrix> sum;
    
    auto fespot = GetPotentialSpace(false);
    
    Array<short> classnr(ma->GetNE(VOL));
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    
    IntegrationRule irsegm(ET_SEGM, 2*order+10);
    
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        
        auto & felpot = dynamic_cast<const ScalarFiniteElement<3>&> (fespot->GetFE (ei, lh));
        auto & fel = dynamic_cast<const HCurlCellFiniteElement<3>&> (GetFE (ei, lh));
        auto & trafo = ma->GetTrafo(ei, lh);
        L2HighOrderFE<ET_SEGM> fel1d(order);
        L2HighOrderFE<ET_SEGM> fel1dxy(order-1);
        IntegrationRule irsegm(ET_SEGM, 2*order+6);
        
        
        Matrix<> gradmat(fel.GetNDof(), felpot.GetNDof());
        Matrix<> mixedmass(fel.GetNDof(), fel.GetNDof());
        
        gradmat = 0.;
        mixedmass = 0.;

        Matrix shape(fel.GetNDof(), 3);        
        Vector shape1d(fel1d.GetNDof());
        Vector shape1dx(fel1dxy.GetNDof());
        Vector shape1dy(fel1dxy.GetNDof());
        Vector shape2d(fel1d.GetNDof()*fel1dxy.GetNDof());
        Vector shape3d(fel1d.GetNDof()*sqr(fel1dxy.GetNDof()));
        
        Vector shapepot(felpot.GetNDof());
        Matrix dshapepot(felpot.GetNDof(), 3);        
        
        bool altshapesE = true;
        
        for (int e = 0; e < 6; e++)
          for (int dir = 0; dir < 2; dir++)
            {
              int ii = (e*2+dir) * (order+1)*(order+1)*(order+1);
              
              IntRange range_edge(ii, ii+(order+1));            
              IntRange range_face(ii+(order+1), ii+(order+1)+2*order*(order+1));
              IntRange range_el(ii+(order+1)+2*order*(order+1), ii+(order+1)*(order+1)*(order+1));
              // cout << "ranges edge/face/el = " << range_edge << range_face << range_el << endl;
              
              auto [v0,v1] = ET_trait<ET_TET>::GetEdge(e);
              if (dir == 1) Swap (v0,v1);

              // edge integrals:
              
              Vec<3> points3d[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3>(0,0,1), Vec<3>(0,0,0) };
              Vec<3> p13d = points3d[v0];
              Vec<3> p23d = points3d[v1];
              
              for (auto ip : irsegm)
                {
                  Vec<3> tau3d = p23d-p13d;
                  IntegrationPoint ip3d(p13d+0.5*ip(0)*tau3d);
                  
                  fel.CalcShape(ip3d, shape, altshapesE);
                  felpot.CalcDShape(ip3d, dshapepot); 
                  fel1d.CalcShape (ip, shape1d);
                  
                  gradmat.Rows(range_edge) += ip.Weight()*shape1d * Trans(dshapepot * tau3d);
                  mixedmass.Rows(range_edge) += ip.Weight()*shape1d * Trans(shape * tau3d);
                }

              // face integrals
              for (int side = 0; side < 2; side++)
                for (auto ipx : irsegm)
                  for (auto ipz : irsegm)
                    {
                      Vec<3> p3d = (side==1) ? Vec<3>(ipx(0), 0, ipz(0)) : Vec<3>(0, ipx(0), ipz(0));

                      Vec<3> ptet = MapHex2Tet (p3d, v0, v1);
                      Mat<3,3> Dptet = DMapHex2Tet (p3d, v0, v1);

                      fel1dxy.CalcShape(ipx, shape1dx);
                      fel1d.CalcShape(ipz, shape1d);
                      for (int k = 0, kk=0; k < shape1d.Size(); k++)
                        for (int l = 0; l < shape1dx.Size(); l++, kk++)
                          shape2d(kk) = shape1d(k)*shape1dx(l);

                      fel.CalcShape (IntegrationPoint(ptet), shape, altshapesE);
                      felpot.CalcDShape (IntegrationPoint(ptet), dshapepot);

                      Vec<3> tau3d = Dptet.Col(2);
                      double w = ipx.Weight()*ipz.Weight();
                      gradmat.Rows(range_face+side*order*(order+1)) += w*shape2d * Trans(dshapepot * tau3d);
                      mixedmass.Rows(range_face+side*order*(order+1)) += w*shape2d * Trans(shape * tau3d);
                    }
                for (auto ipx : irsegm)
                  for (auto ipy : irsegm)
                    for (auto ipz : irsegm)
                      {
                        Vec<3> p3d(ipx(0), ipy(0), ipz(0));

                        Vec<3> ptet = MapHex2Tet (p3d, v0, v1);
                        Mat<3,3> Dptet = DMapHex2Tet (p3d, v0, v1);
                        
                        fel1dxy.CalcShape(ipx, shape1dx);
                        fel1dxy.CalcShape(ipy, shape1dy);
                        fel1d.CalcShape(ipz, shape1d);
                        for (int k = 0, kk=0; k < shape1d.Size(); k++)
                          for (int l = 0; l < shape1dx.Size(); l++)
                            for (int m = 0; m < shape1dy.Size(); m++, kk++)
                              shape3d(kk) = shape1d(k)*shape1dx(l)*shape1dy(m);
                        
                        fel.CalcShape (IntegrationPoint(ptet), shape, altshapesE);
                        felpot.CalcDShape (IntegrationPoint(ptet), dshapepot);
                        
                        Vec<3> tau3d = Dptet.Col(2);
                        double w = ipx.Weight()*ipy.Weight()*ipz.Weight();
                        gradmat.Rows(range_el) += w*shape3d * Trans(dshapepot * tau3d);
                        mixedmass.Rows(range_el) += w*shape3d * Trans(shape * tau3d);
                      }
            } 
        
        // cout << "gradmat = " << gradmat.Rows(0, gradmat.Height()/12) << endl;
        // cout << "mixedmass = " << mixedmass.Rows(0, gradmat.Height()/12) << endl;
        // continue;
        
        CalcInverse (mixedmass);
        Matrix<> finalmat(fel.GetNDof(), felpot.GetNDof());
        finalmat = mixedmass * gradmat;
        // cout << "finalmat = " << finalmat << endl;

        
        Table<DofId> xdofs(elclass_inds.Size(), felpot.GetNDof()),
          ydofs(elclass_inds.Size(), fel.GetNDof());
        
        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            GetDofNrs(ei, dnumsy);
            fespot->GetDofNrs(ei, dnumsx);
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }
        
        auto mat = make_shared<ConstantElementByElementMatrix>
          (this->GetNDof(), fespot->GetNDof(),
           finalmat, std::move(ydofs), std::move(xdofs));
        
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }
    
    auto diag = make_shared<VVector<double>> (this->GetNDof());
    FlatVector<> fd = diag->FV();
    fd = 0.0;
    Array<DofId> dnums;
    for (size_t i = 0; i < ma->GetNE(VOL); i++)
      {
        GetDofNrs( { VOL, i }, dnums);
        for (auto d : dnums)
          fd[d] += 1;
      }
    for (auto & val : fd)
      val = 1/val;
    auto diagmat = make_shared<DiagonalMatrix<double>>(diag);
    return ComposeOperators(diagmat, sum);
  }
  

  
  shared_ptr<BaseMatrix> HCurlDualCells::GetGradientOperator(bool nanocells) const
  {
    if (ma -> GetDimension() == 2)
      {
        cout << "GetGradient called" << endl;
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;

        auto fespot = GetPotentialSpace(false);
        
        Array<short> classnr(ma->GetNE(VOL));
        ma->IterateElements
          (VOL, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
    
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();

        IntegrationRule irquad(ET_QUAD, 2*order+6);
        IntegrationRule irsegm(ET_SEGM, 2*order+6);
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(VOL,elclass_inds[0]);
                                     
            auto & felpot = dynamic_cast<const ScalarFiniteElement<2>&> (fespot->GetFE (ei, lh));
            auto & fel = dynamic_cast<const HCurlCellFiniteElement<2>&> (GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);

            // cout << "fel.ndof = " << fel.GetNDof() << endl;
            // cout << "felc.ndof = " << felc.GetNDof() << endl;
        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TRIG> (order+1);
            auto [edges,nedge] = GetMicroCellEdges<ET_TRIG> (order+1);            

            // cout << "edges.size = " << edges.Size() << ", nedge = " << nedge << endl;
            // cout << "felpot.ndof = " << felpot.GetNDof() << endl;
            // cout << "fel.ndof = " << fel.GetNDof() << endl;
            Matrix<> gradmat(nedge, felpot.GetNDof());
            Matrix<> mixedmass(nedge, fel.GetNDof());

            gradmat = 0.;
            mixedmass = 0.;
            Matrix shape(fel.GetNDof(), 2);
            Vector shapepot(felpot.GetNDof());
              

            // gradmat ...
            for (auto i : Range(nedge))
              {
                auto [p0,p1] = edges[i];
                
                // gradmat(i, p0) = -1;
                // gradmat(i, p1) = 1;
                
                Vec<2> micropts[2] = { pnts[p0], pnts[p1] };

                IntegrationPoint ip0(micropts[0](0), micropts[0](1), 0, 0);                
                felpot.CalcShape (ip0, shapepot);
                gradmat.Row(i) = shapepot;
                IntegrationPoint ip1(micropts[1](0), micropts[1](1), 0, 0);                
                felpot.CalcShape (ip1, shapepot);
                gradmat.Row(i) -= shapepot;                
                
                Vec<2> tau = micropts[1]-micropts[0];
                for (auto ip : irsegm)
                  {
                    Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                    IntegrationPoint ipseg(p(0), p(1), 0, 0);
                    fel.CalcAltShape(ipseg, shape);
                    mixedmass.Row(i) += ip.Weight() * shape*tau;
                  }
              }

            // cout << "gradmat = " << gradmat << endl;
            // cout << "mixedmass = " << mixedmass << endl;
            // continue;

            CalcInverse (mixedmass);
            Matrix<> tmp(fel.GetNDof(), felpot.GetNDof());
            tmp = mixedmass * gradmat;

            Table<DofId> xdofs(elclass_inds.Size(), felpot.GetNDof()),
              ydofs(elclass_inds.Size(), fel.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsy);
                fespot->GetDofNrs(ei, dnumsx);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            auto mat = make_shared<ConstantElementByElementMatrix>
              (this->GetNDof(), fespot->GetNDof(),
               tmp, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        auto diag = make_shared<VVector<double>> (this->GetNDof());
        FlatVector<> fd = diag->FV();
        fd = 0.0;
        Array<DofId> dnums;
        for (size_t i = 0; i < ma->GetNE(VOL); i++)
          {
            GetDofNrs( { VOL, i }, dnums);
            for (auto d : dnums)
              fd[d] += 1;
          }
        for (auto & val : fd)
          val = 1/val;
        auto diagmat = make_shared<DiagonalMatrix<double>>(diag);
        return ComposeOperators(diagmat, sum);
      }

    if (ma -> GetDimension() == 3)
      {
        if (nanocells == false)
          {
            // if (dual == true)
            // throw Exception("dual available only for nanocells");
            return GetGradientOperator3D();
          }


        
        LocalHeap lh(10*1000*1000);
        shared_ptr<BaseMatrix> sum;

        auto fespot = GetPotentialSpace(false);
        
        Array<short> classnr(ma->GetNE(VOL));
        ma->IterateElements
          (VOL, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
    
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();

        IntegrationRule irsegm(ET_SEGM, 2*order+6);
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            ElementId ei(VOL,elclass_inds[0]);
                                     
            auto & felpot = dynamic_cast<const ScalarFiniteElement<3>&> (fespot->GetFE (ei, lh));
            auto & fel = dynamic_cast<const HCurlCellFiniteElement<3>&> (GetFE (ei, lh));
            auto & trafo = ma->GetTrafo(ei, lh);
        
            auto [pnts,ind,num] = GetMicroCellPoints<ET_TET> (order+1);
            auto [edges,nedge] = GetMicroCellEdges<ET_TET> (order+1);            
            Matrix<> gradmat(nedge, felpot.GetNDof());
            Matrix<> mixedmass(nedge, fel.GetNDof());

            gradmat = 0.;
            mixedmass = 0.;
            Matrix shape(fel.GetNDof(), 3);
            Vector shapepot(felpot.GetNDof());
              
            // gradmat ...
            for (auto i : Range(nedge))
              {
                auto [p0,p1] = edges[i];
                
                // gradmat(i, p0) = -1;
                // gradmat(i, p1) = 1;
                
                Vec<3> micropts[2] = { pnts[p0], pnts[p1] };

                IntegrationPoint ip0(micropts[0](0), micropts[0](1), micropts[0](2), 0);                
                felpot.CalcShape (ip0, shapepot);
                gradmat.Row(i) = shapepot;
                IntegrationPoint ip1(micropts[1](0), micropts[1](1), micropts[1](2), 0);                
                felpot.CalcShape (ip1, shapepot);
                gradmat.Row(i) -= shapepot;                
                
                Vec<3> tau = micropts[1]-micropts[0];
                for (auto ip : irsegm)
                  {
                    Vec<3> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                    IntegrationPoint ipseg(p(0), p(1), p(2), 0);
                    fel.CalcAltShape(ipseg, shape);
                    mixedmass.Row(i) += ip.Weight() * shape*tau;
                  }
              }

            // cout << "gradmat = " << gradmat << endl;
            // cout << "mixedmass = " << mixedmass << endl;
            // continue;

            CalcInverse (mixedmass);
            Matrix<> tmp(fel.GetNDof(), felpot.GetNDof());
            tmp = mixedmass * gradmat;

            Table<DofId> xdofs(elclass_inds.Size(), felpot.GetNDof()),
              ydofs(elclass_inds.Size(), fel.GetNDof());
        
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                GetDofNrs(ei, dnumsy);
                fespot->GetDofNrs(ei, dnumsx);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }

            auto mat = make_shared<ConstantElementByElementMatrix>
              (this->GetNDof(), fespot->GetNDof(),
               tmp, std::move(ydofs), std::move(xdofs));
        
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }

        auto diag = make_shared<VVector<double>> (this->GetNDof());
        FlatVector<> fd = diag->FV();
        fd = 0.0;
        Array<DofId> dnums;
        for (size_t i = 0; i < ma->GetNE(VOL); i++)
          {
            GetDofNrs( { VOL, i }, dnums);
            for (auto d : dnums)
              fd[d] += 1;
          }
        for (auto & val : fd)
          val = 1/val;
        auto diagmat = make_shared<DiagonalMatrix<double>>(diag);
        return ComposeOperators(diagmat, sum);
      }
    return nullptr;
  }
}
