#include <comp.hpp>
#include "hcurlprimalcells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"


namespace ngcomp
{
  class HCurlPrimalCellTet : public HCurlCellFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIR;
    const IntegrationRule & GaussRadauIRinv;
    const IntegrationRule & tangentialIR;
  public:
    HCurlPrimalCellTet (const IntegrationRule & _GaussRadauIR,
                        const IntegrationRule & _GaussRadauIRinv,
                        const IntegrationRule & _tangentialIR)
      : HCurlCellFiniteElement<3> (_tangentialIR.Size() * (4+6*2*(_GaussRadauIR.Size()-1)+4*3*sqr(_GaussRadauIR.Size()-1)),
                                   _GaussRadauIR.Size()-1),
      GaussRadauIR(_GaussRadauIR), GaussRadauIRinv(_GaussRadauIRinv), tangentialIR(_tangentialIR)
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
      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
        {
          if (v != maxlam)
            {
              IVec<3> ind = { 0, 0, 0 };
              int dirv2 = vdir[v];
              for (int k = 0; k < ndt; k++)
                {
                  ind[dirv2] = k;
                  assign(ii+k, ind, dirv2);
                }
            }
          ii += ndt;
        }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
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
                                BareSliceMatrix<> curlshape) const override
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
      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
        {
          if (v != maxlam)
            {
              IVec<3> ind = { 0, 0, 0 };
              int dirv2 = vdir[v];
              for (int k = 0; k < ndt; k++)
                {
                  ind[dirv2] = k;
                  assign(ii+k, ind, dirv2);
                }
            }
          ii += ndt;
        }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
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


    
    virtual void CalcAltCurlShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> curlshape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //curlshape = 0;
      curlshape.AddSize(ndof,3) = 0;

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
      
      Mat<3> trafo = Inv(F2*F);

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
          Vec<3> xinode(GaussRadauIR[i[0]](0), GaussRadauIR[i[1]](0), GaussRadauIR[i[2]](0));
          Mat<3,3> Fnode = DMapHex2Tet(xinode);          
          Mat<3,3> trafonode = Trans(Inv(F2*Fnode));

          Mat<3,3> Dshapeorig = 0.;
          switch (dir)
            {
            case 0:
            {
            Dshapeorig.Col(0) = Trans(trafo)*Vec<3>(polxitang[i[0]].DValue(0)*poleta[i[1]].Value()*polzeta[i[2]].Value(),
              polxitang[i[0]].Value()*poleta[i[1]].DValue(0)*polzeta[i[2]].Value(),
              polxitang[i[0]].Value()*poleta[i[1]].Value()*polzeta[i[2]].DValue(0));
              break;
            }
            case 1:
            {
            Dshapeorig.Col(1) = Trans(trafo)*Vec<3>(polxi[i[0]].DValue(0)*poletatang[i[1]].Value()*polzeta[i[2]].Value(),
              polxi[i[0]].Value()*poletatang[i[1]].DValue(0)*polzeta[i[2]].Value(),
              polxi[i[0]].Value()*poletatang[i[1]].Value()*polzeta[i[2]].DValue(0));
              break;
            }
            case 2:
            {
            Dshapeorig.Col(2) = Trans(trafo)*Vec<3>(polxi[i[0]].DValue(0)*poleta[i[1]].Value()*polzetatang[i[2]].Value(),
              polxi[i[0]].Value()*poleta[i[1]].DValue(0)*polzetatang[i[2]].Value(),
              polxi[i[0]].Value()*poleta[i[1]].Value()*polzetatang[i[2]].DValue(0));
              break;
            }
            break;
            }
          Mat<3,3> Dshape = trafonode*Trans(Dshapeorig);
          curlshape.Row(nr) = Vec<3>(Dshape(2,1)-Dshape(1,2),-Dshape(2,0)+Dshape(0,2),Dshape(1,0)-Dshape(0,1));

        };


        

      
      int ii = 0;
      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
        {
          if (v != maxlam)
            {
              IVec<3> ind = { 0, 0, 0 };
              int dirv2 = vdir[v];
              for (int k = 0; k < ndt; k++)
                {
                  ind[dirv2] = k;
                  assign(ii+k, ind, dirv2);
                }
            }
          ii += ndt;
        }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
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
                            BareSliceMatrix<> shape) const override
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
      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
        {
          if (v != maxlam)
            {
              IVec<3> ind = { 0, 0, 0 };
              int dirv2 = vdir[v];
              for (int k = 0; k < ndt; k++)
                {
                  ind[dirv2] = k;
                  assign(ii+k, ind, dirv2);
                }
            }
          ii += ndt;
        }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
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
                                 BareSliceMatrix<> shape) const override
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
          // Mat<3> trafo = Trans(Inv(F2*F));
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

      /*
      ArrayMem<double, 20> polxiinv(ndGR), poletainv(ndGR), polzetainv(ndGR);
      LagrangePolynomials(xi(0), GaussRadauIRinv, polxiinv);
      LagrangePolynomials(xi(1), GaussRadauIRinv, poletainv);
      LagrangePolynomials(xi(2), GaussRadauIRinv, polzetainv);

      auto assign =  [&](int nr, IVec<3> i, int dir)
        {
          Vec<3> xinode(dir == 0 ? GaussRadauIRinv[i[0]](0) : GaussRadauIR[i[0]](0),
                        dir == 1 ? GaussRadauIRinv[i[1]](0) : GaussRadauIR[i[1]](0),
                        dir == 2 ? GaussRadauIRinv[i[2]](0) : GaussRadauIR[i[2]](0));

          Mat<3,3> F = DMapHex2Tet(xinode);          
          // Mat<3> trafo = Trans(Inv(F2*F));
          Mat<3> trafo1 = Det(F2*F) *Inv(F2*F) * Trans(Inv(F2*F));
          Mat<3> trafo = trafo * trafo1;

          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafo*Vec<3>(polxiinv[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafo*Vec<3>(0, polxi[i[0]]*poletainv[i[1]]*polzeta[i[2]], 0);
              break;
            }
            case 2:
            {
              shape.Row(nr) = trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzetainv[i[2]]);
              break;
            }
            break;
            }
        };
      */  
      
      int ii = 0;
      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
        {
          if (v != maxlam)
            {
              IVec<3> ind = { 0, 0, 0 };
              int dirv2 = vdir[v];
              for (int k = 0; k < ndt; k++)
                {
                  ind[dirv2] = k;
                  assign(ii+k, ind, dirv2);
                }
            }
          ii += ndt;
        }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
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

    
  };

#define GONETOHEADER
#ifdef GONETOHEADER
  class DiffOpAltShape : public DiffOp<DiffOpAltShape>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      static_cast<const HCurlPrimalCellTet&>(fel).CalcAltShape (mip.IP(), Trans(mat));
      Mat<3,3> F = mip.GetJacobian();
      Mat<3,3> trafo = Trans(Inv(F));
      for (int i = 0; i < mat.Width(); i++)
        {
          Vec<3> shape = mat.Col(i);
          mat.Col(i) = trafo * shape;
        }
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
      Mat<D,D> F = static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobian();
      mat = 1/Det(F) * F;   
    }
    
  };

#endif
  

  

  HCurlPrimalCells::
  HCurlPrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {
    collocated = !flags.GetDefineFlagX("collocated").IsFalse();
    // cout << "HCurlPrimal ctor, collocated = " << collocated << endl;
    
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
    additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHCurl<3>>> ());    
    additional_evaluators.Set ("Hodge", make_shared<T_DifferentialOperator<DiffOpHodgeHCurl<3>>> ());        
    additional_evaluators.Set ("Piolashape", make_shared<T_DifferentialOperator<DiffOpPiolaShapeHCurl<3>>> ());        
    additional_evaluators.Set ("altcurl", make_shared<T_DifferentialOperator<DiffOpAltCurlHCurl<3>>> ());        
    
    
    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      {
        GaussRadauIR.Append (IntegrationPoint (1-xi[i]-1e-10, 0, 0, wi[i]));
        GaussRadauIRinv.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));        
        tangentialIR.Append (IntegrationPoint (1-xi[i]-1e-10, 0, 0, wi[i]));
      }
  }


  void HCurlPrimalCells :: Update()
  {
    FESpace::Update();
    size_t ndof = 0;

    int ndt = tangentialIR.Size();
    int ndGR = GaussRadauIR.Size();

    first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
    first_cell_dofs[0] = ndof;
    for (auto i : Range(ma->GetNE(VOL)))
      first_cell_dofs[i+1] = ndof += ndt*(4+6*2*(ndGR-1)+4*3*sqr(ndGR-1));

    SetNDof(ndof);
  }
    
  void HCurlPrimalCells ::
  GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    if (ei.VB() != VOL) return;
    dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
  }


  FiniteElement & HCurlPrimalCells ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) HCurlPrimalCellTet(GaussRadauIR, GaussRadauIRinv, GaussRadauIR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      default:
        throw Exception("element not implemented");
      }
  }

  shared_ptr<BaseMatrix> HCurlPrimalCells ::
  GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
  {
    // cout << "get hcurl mass" << endl;
    if (!collocated)
      throw Exception("only for collocated spaces");


   static Timer t("GetMassOperator"); RegionTimer reg(t);
    static Timer tint("integrate");
    static Timer tsp("MakeSparse");
      
    HeapReset hr(lh);
    auto irs = this->GetIntegrationRules();
    IntegrationRule ir = std::move(irs[ET_TET]);

    auto & felref = dynamic_cast<const HCurlPrimalCellTet&> (GetFE(ElementId(VOL,0), lh));
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
        auto & felref = dynamic_cast<const HCurlPrimalCellTet&> (GetFE(ElementId(VOL,nr), lh));
        
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
            
            for (auto i : myrange)
              {
                HeapReset hr(mylh);
                auto nr = elclass_inds[i];

                auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                MappedIntegrationRule<3,3> mir(ir, trafo, mylh);
                GetDofNrs(ElementId(VOL,nr), dofs);
                
                for (size_t i = 0; i < mir.Size(); i++)
                  {
                    Mat<3,3> rhoi = Id<3>();
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

    
    /*
    static Timer t("GetMassOperator"); RegionTimer reg(t);
    static Timer tint("integrate");
    static Timer tsp("MakeSparse");
      
    HeapReset hr(lh);
    auto irs = this->GetIntegrationRules();
    IntegrationRule ir = std::move(irs[ET_TET]);

    auto & felref = dynamic_cast<const HCurlPrimalCellTet&> (GetFE(ElementId(VOL,0), lh));
    Matrix shapes(felref.GetNDof(), 3*ir.Size());
    Matrix rhoi_shapes(felref.GetNDof(), 3*ir.Size());

    for (int i = 0; i < ir.Size(); i++)
      felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));

    for (int i = 0; i < shapes.Height(); i++)
      for (int j = 0; j < shapes.Width(); j++)
        if (fabs(shapes(i,j)) < 1e-10)
          shapes(i,j) = 0;
    rhoi_shapes = shapes;
    
    Array<int> rowind, colind;
    Array<double> values;
    Matrix<> elmat(felref.GetNDof());
    Array<DofId> dofs(felref.GetNDof());

    tint.Start();
    for (size_t nr : Range(ma->GetNE()))
      {
        HeapReset hr(lh);
        
        
        auto & felref = dynamic_cast<const HCurlPrimalCellTet&> (GetFE(ElementId(VOL,nr), lh));
        Matrix shapes(felref.GetNDof(), 3*ir.Size());
        Matrix rhoi_shapes(felref.GetNDof(), 3*ir.Size());
        
        for (int i = 0; i < ir.Size(); i++)
          felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
        
        for (int i = 0; i < shapes.Height(); i++)
          for (int j = 0; j < shapes.Width(); j++)
            if (fabs(shapes(i,j)) < 1e-10)
              shapes(i,j) = 0;

        SmallSparseMatrix spshapes(shapes, 1e-8);        
        rhoi_shapes = shapes;

        
        
        MappedIntegrationRule<3,3> mir(ir, ma->GetTrafo(ElementId(VOL,nr), lh), lh);
        GetDofNrs(ElementId(VOL,nr), dofs);
        
        for (size_t i = 0; i < mir.Size(); i++)
          {
            Mat<3,3> rhoi = Id<3>();
            Mat<3,3> Finv = mir[i].GetJacobianInverse();
            rhoi = Finv * rhoi * Trans(Finv);
            rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
            rhoi_shapes.Cols(3*i, 3*i+3) = shapes.Cols(3*i, 3*i+3) * rhoi;
          }
        // elmat = rhoi_shapes * Trans(shapes);
        spshapes.Mult (Trans(rhoi_shapes), Trans(elmat));
        
        // *testout << "elmat = " << nr << endl << elmat << endl;
        for (int i = 0; i < dofs.Size(); i++)
          for (int j = 0; j < dofs.Size(); j++)
            if (fabs(elmat(i,j)) > 1e-10)
              {
                rowind.Append(dofs[i]);
                colind.Append(dofs[j]);
                values.Append(elmat(i,j));
              }
        }
      tint.Stop();

      // for (int i = 0; i < rowind.Size(); i++)
      // *testout << rowind[i] << " " << colind[i] << " " << values[i] << endl;
      */

    
      tsp.Start();
      auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
      tsp.Stop();
      
      return make_shared<MySuperSparseMatrix> (std::move(*spmat));
    }


  
  std::map<ELEMENT_TYPE, IntegrationRule> HCurlPrimalCells::GetIntegrationRules(bool fix_lo) const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    
    IntegrationRule irtet;
    irtet.SetDim(3);
    int i=0;
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

    if (GaussRadauIR.Size()==1 && fix_lo)
      for (auto & ip : irtet)
        ip.SetWeight(1.0/24);   // fix consistency for lowest order
    
    rules[ET_TET] = std::move(irtet);


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
    if (GaussRadauIR.Size()==1 && fix_lo)
      for (auto & ip : irtrig)
        ip.SetWeight(1.0/6);   // fix consistency for lowest order
    
    rules[ET_TRIG] = std::move(irtrig);

    



    
    return rules;
  }    
}
