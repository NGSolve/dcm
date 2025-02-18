#include <comp.hpp>
#include "h1primalcells.hpp"
#include "hcurlcells.hpp"
#include "hcurldualcells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"

namespace ngcomp
{

  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;

  H1PrimalCellSegm :: H1PrimalCellSegm (int order, const IntegrationRule & _GaussRadauIR)
    : ScalarFiniteElement<1> (2*order+1, order), GaussRadauIR(_GaussRadauIR)
  { ; }

  void H1PrimalCellSegm :: CalcShape (const IntegrationPoint & ip, 
      BareSliceVector<> shape) const
  {
    ArrayMem<double, 20> polx(order+1);

    double lam[] = { ip(0), 1-ip(0) };
    int maxlam = 0;
    if (lam[1] > lam[0])
      maxlam = 1;

    shape.Range(GetNDof()) = 0;

    double lx = lam[(maxlam+1)%2];
    double xi = 2*lx;

    LagrangePolynomials(xi, GaussRadauIR, polx);
    shape(0) = polx[0];

    int ii = 1;
    IVec<2> e = GetVertexOrientedEdge(0);
    for (int j = 0; j < 2; j++)
    {
      if (e[j] == maxlam)
      {
        for (int k = 0; k < order; k++)
          shape(ii+k) = polx[k+1];
      }
      ii+=order;
    }
  }

  void H1PrimalCellSegm :: CalcDShape (const IntegrationPoint & ip, 
      BareSliceMatrix<> baredshape) const
  {
    ArrayMem<AutoDiff<1>, 20> polx(order+1);

    auto dshape = baredshape.AddSize(GetNDof(), 1);
    dshape = 0;

    double lam[] = { ip(0), 1-ip(0) };
    int maxlam = 0;
    if (lam[1] > lam[0])
      maxlam = 1;

    double lx = lam[(maxlam+1)%2];
    double xi = 2*lx;

    AutoDiff<1> adxi(xi, 0);
    LagrangePolynomials(adxi, GaussRadauIR, polx);
    double sig = 1;
    if (maxlam == 0)
      sig = -1;
    dshape(0,0) = 2*sig*polx[0].DValue(0);

    int ii = 1;
    IVec<2> e = GetVertexOrientedEdge(0);
    for (int j = 0; j < 2; j++)
    {
      if (e[j] == maxlam)
      {
        for (int k = 0; k < order; k++)
          dshape(ii+k,0) = 2*sig*polx[k+1].DValue(0);
      }
      ii+=order;
    }
  }


  void H1PrimalCellTrig ::
    CalcShape (const IntegrationPoint & ip, 
        BareSliceVector<> shape) const 
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.Range(GetNDof()) = 0;

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

      ArrayMem<double, 20> polxi(order+1), poleta(order+1);
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      auto assign =  [&](int nr, IVec<2> ind)
      {
        shape(nr) = polxi[ind[0]]*poleta[ind[1]];
      };

      assign(0, { 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 3; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
      }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 3; v1++)
        for (int v2 = 0; v2 < v1; v2++)
        {
          if (v1 != maxlam && v2 != maxlam)
          {
            IVec<2> ind = { 0, 0 };
            int dirv1 = vdir[v1];
            int dirv2 = vdir[v2];
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }
    }

  void H1PrimalCellTrig ::
    CalcL2Shape (const IntegrationPoint & ip, 
        BareSliceVector<> shape) const 
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.Range(GetNDof()) = 0;

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

      double trafo = 1./fabs(Det(F2*F));
      ArrayMem<double, 20> polxi(order+1), poleta(order+1);
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      auto assign =  [&](int nr, IVec<2> ind)
      {
        shape(nr) = trafo*polxi[ind[0]]*poleta[ind[1]];
      };

      assign(0, { 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 3; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
      }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 3; v1++)
        for (int v2 = 0; v2 < v1; v2++)
        {
          if (v1 != maxlam && v2 != maxlam)
          {
            IVec<2> ind = { 0, 0 };
            int dirv1 = vdir[v1];
            int dirv2 = vdir[v2];
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }
    }

  void H1PrimalCellTrig ::
    CalcL2AltShape (const IntegrationPoint & ip, 
        BareSliceVector<> shape) const 
    {
      //just L2 mapping for now!!
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.Range(GetNDof()) = 0;

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

      Mat<2,2> F2;    // trafo from vertex permutation
      Vec<2> verts[] = { Vec<2>(1,0), Vec<2>(0,1), Vec<2> (0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[maxvi]-verts[maxlam];


      ArrayMem<double, 20> polxi(order+1), poleta(order+1);
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);

      auto assign =  [&](int nr, IVec<2> ind)
      {
        Vec<2> xinode(GaussRadauIR[ind[0]](0), GaussRadauIR[ind[1]](0));
        Mat<2,2> F = DMapQuad2Trig(xinode);
        double trafo = 1/fabs(Det(F2*F));
        shape(nr) = trafo*polxi[ind[0]]*poleta[ind[1]];
      };

      assign(0, { 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 3; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
      }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 3; v1++)
        for (int v2 = 0; v2 < v1; v2++)
        {
          if (v1 != maxlam && v2 != maxlam)
          {
            IVec<2> ind = { 0, 0 };
            int dirv1 = vdir[v1];
            int dirv2 = vdir[v2];
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }
    }

  void H1PrimalCellTrig ::
    CalcDShape (const IntegrationPoint & ip, 
        BareSliceMatrix<> baredshape) const 
    {
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      auto dshape = baredshape.AddSize(GetNDof(), 2);
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

      ArrayMem<AutoDiff<1>, 20> polxi(order+1), poleta(order+1);
      LagrangePolynomials(AutoDiff<1>(xi(0),0), GaussRadauIR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), GaussRadauIR, poleta);

      auto assign =  [&](int nr, IVec<2> i)
      {
        dshape.Row(nr) =
          trafo*Vec<2>(polxi[i[0]].DValue(0)*poleta[i[1]].Value(),
              polxi[i[0]].Value()*poleta[i[1]].DValue(0));
      };


      assign(0, { 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 3; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
      }

      // internal faces into direction of vertices v1, v2
      for (int v1 = 0; v1 < 3; v1++)
        for (int v2 = 0; v2 < v1; v2++)
        {
          if (v1 != maxlam && v2 != maxlam)
          {
            IVec<2> ind = { 0, 0 };
            int dirv1 = vdir[v1];
            int dirv2 = vdir[v2];
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }
    }

  void H1PrimalCellTrig ::
    Interpolate (const ElementTransformation & trafo, 
        const class CoefficientFunction & func, SliceMatrix<> coefs,
        LocalHeap & lh) const 
    {
      // auto [pnts, inds, nump] = GetMicroCellPoints<ET_TRIG> (order, true);

      IntegrationRule irtrig;
      for (int v = 0; v < 3; v++)
        for (auto ipx : GaussRadauIR)
          for (auto ipy : GaussRadauIR)
          {
            if (ipx(0) > 1-1e-13) continue;
            Vec<2> xi(ipx(0), ipy(0));
            Vec<2> x = MapQuad2Trig (xi, v);
            Mat<2,2> jac = DMapQuad2Trig(xi, v);
            double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
            irtrig += IntegrationPoint(x(0), x(1), 0, w);
          }
      irtrig += IntegrationPoint(1.0/3, 1.0/3, 0, 0);

      Matrix mat(this->ndof, this->ndof);
      Vector rhs(this->ndof);
      Vec<1> fi;
      for (int i = 0; i < this->ndof; i++)
      {
        // Vec<2> x = pnts[i];
        // IntegrationPoint ip(x(0), x(1), 0, 0);
        MappedIntegrationPoint<2,2> mip(irtrig[i], trafo);
        func.Evaluate (mip, fi);
        rhs(i) = fi(0);
        this->CalcShape(irtrig[i], mat.Row(i));
      }
      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }          



  class H1PrimalCellTet : public ScalarFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & GaussRadauIR;
    public:
    H1PrimalCellTet (int order, const IntegrationRule & _GaussRadauIR)
      : ScalarFiniteElement<3> (1+4*order+6*sqr(order)+4*sqr(order)*order, order),
      GaussRadauIR(_GaussRadauIR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }


    virtual void CalcShape (const IntegrationPoint & ip, 
        BareSliceVector<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      shape.Range(GetNDof()) = 0;

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

      ArrayMem<double, 20> polxi(order+1), poleta(order+1), polzeta(order+1);   
      LagrangePolynomials(xi(0), GaussRadauIR, polxi);
      LagrangePolynomials(xi(1), GaussRadauIR, poleta);
      LagrangePolynomials(xi(2), GaussRadauIR, polzeta);

      auto assign =  [&](int nr, IVec<3> ind)
      {
        shape(nr) = polxi[ind[0]]*poleta[ind[1]]*polzeta[ind[2]];
      };


      assign(0, { 0, 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
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
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }

      ii += maxlam*order*order*order;
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          for (int k = 0; k < order; k++)
            assign(ii++, { i+1, j+1, k+1 });
    }


    virtual void CalcDShape (const IntegrationPoint & ip, 
        BareSliceMatrix<> baredshape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      auto dshape = baredshape.AddSize(GetNDof(), 3);
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


      ArrayMem<AutoDiff<1>, 20> polxi(order+1), poleta(order+1), polzeta(order+1);   
      LagrangePolynomials(AutoDiff<1>(xi(0),0), GaussRadauIR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), GaussRadauIR, poleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), GaussRadauIR, polzeta);

      auto assign =  [&](int nr, IVec<3> i)
      {
        dshape.Row(nr) =
          trafo*Vec<3>(polxi[i[0]].DValue(0)*poleta[i[1]].Value()*polzeta[i[2]].Value(),
              polxi[i[0]].Value()*poleta[i[1]].DValue(0)*polzeta[i[2]].Value(),
              polxi[i[0]].Value()*poleta[i[1]].Value()*polzeta[i[2]].DValue(0));
      };



      assign(0, { 0, 0, 0 });

      int ii = 1;

      // internal edges into direction of vertex v
      for (int v = 0; v < 4; v++)
      {
        if (v != maxlam)
        {
          IVec<3> ind = { 0, 0, 0 };
          int dirv2 = vdir[v];
          for (int k = 0; k < order; k++)
          {
            ind[dirv2] = k+1;
            assign(ii+k, ind);
          }
        }
        ii += order;
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
            for (int k = 0, kk=ii; k < order; k++)
              for (int l = 0; l < order; l++, kk++)
              {
                ind[dirv1] = k+1;
                ind[dirv2] = l+1;
                assign(kk, ind);
              }
          }
          ii += sqr(order);
        }

      ii += maxlam*order*order*order;
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          for (int k = 0; k < order; k++)
            assign(ii++, { i+1, j+1, k+1 });
    } 
  };


  H1PrimalCells::
    H1PrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
    {
      switch (ma->GetDimension())
      {
        case 1:
          {
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();

            break;
          }
        case 2:
          {
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();

            additional_evaluators.Set ("dual_mapping", make_shared<T_DifferentialOperator<DiffOpIdDual<2,2>>> ());
            additional_evaluators.Set ("l2shape", make_shared<T_DifferentialOperator<DiffOpL2Shape2D>> ());
            break;
          }
        case 3:
          {
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
            additional_evaluators.Set ("dual_mapping", make_shared<T_DifferentialOperator<DiffOpIdDual<3,3>>> ());        
            break;
          }
      }


      Array<double> xi, wi;
      ComputeGaussRadauRule (order+1, xi, wi);
      for (auto i : Range(xi))
        GaussRadauIR.Append (IntegrationPoint (1-xi[i]-1e-14, 0, 0, wi[i]));
    }


  void H1PrimalCells :: Update()
  {
    FESpace::Update();
    size_t ndof = 0;

    first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
    first_cell_dofs[0] = ndof;
    switch (ma->GetDimension())
    {
      case 1:
        {
          for (auto i : Range(ma->GetNE(VOL)))
            first_cell_dofs[i+1] = ndof += 1+2*order;
          break;
        }
      case 2:
        {
          for (auto i : Range(ma->GetNE(VOL)))
            first_cell_dofs[i+1] = ndof += 1+3*order+3*sqr(order);
          break;
        }
      case 3:
        {
          for (auto i : Range(ma->GetNE(VOL)))
            first_cell_dofs[i+1] = ndof += 1+4*order+6*sqr(order)+4*sqr(order)*order;
          break;
        }
      default:
        throw Exception("invalid dimension");
    }

    SetNDof(ndof);
  }

  void H1PrimalCells ::
    GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
      dnums.SetSize0();
      if (ei.VB() != VOL) return;
      dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
    }


  FiniteElement & H1PrimalCells ::
    GetFE (ElementId ei, Allocator & alloc) const
    {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
      {
        case ET_SEGM:
          {
            auto segm = new (alloc) H1PrimalCellSegm(order, GaussRadauIR);
            segm->SetVertexNumbers (ngel.vertices);
            return *segm;
          }
        case ET_TRIG:
          {
            auto trig = new (alloc) H1PrimalCellTrig(order, GaussRadauIR);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }
        case ET_TET:
          {
            auto tet = new (alloc) H1PrimalCellTet(order, GaussRadauIR);
            tet->SetVertexNumbers (ngel.vertices);
            return *tet;
          }
        default:
          throw Exception("element not implemented");
      }
    }


  shared_ptr<BaseMatrix> H1PrimalCells::GetRotOperator2DNano(bool dual) const
  {
    cout << "careful: only implemented for mesh consisting solely of ref element!!" << endl;
    bool printmats = false;
    LocalHeap lh(10*1000*1000);
    shared_ptr<BaseMatrix> sum;
    Flags hcflags;
    hcflags.SetFlag ("order", order);
    auto feshcurl = make_shared<HCurlDualCells> (ma, hcflags);
    feshcurl->Update();
    feshcurl->FinalizeUpdate();

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

    IntegrationRule irsegm(ET_SEGM, 2*order+4);

    //MAPPINGTYPE mappingE = COVARIANT; // POLYNOMIAL doesn't make any difference
    MAPPINGTYPE mappingD = PIOLAPOLYNOMIAL;
    MAPPINGTYPE mappingE = POLYNOMIAL; // 

    auto [pnts,ind,num] = GetNanoPoints<ET_TRIG> (order+1,false);
    auto [lines,line_ends] = GetNanoLines<ET_TRIG> (order+1,false);

    for (auto elclass_inds : table)
    {
      if (elclass_inds.Size() == 0) continue;

      ElementId ei(VOL,elclass_inds[0]);

      auto & felh1 = dynamic_cast<const ScalarFiniteElement<2>&> (GetFE (ei, lh));
      auto & felhcurl = dynamic_cast<const HCurlCellFiniteElement<2>&> (feshcurl->GetFE (ei, lh));
      auto & trafo = ma->GetTrafo(ei, lh);

      //cout << felhcurl.GetNDof() << endl; 
      Matrix<> curlmat(felhcurl.GetNDof(), felh1.GetNDof());
      Matrix<> mixedmass(felhcurl.GetNDof(), felhcurl.GetNDof());
      Matrix<> lumpedmass(felhcurl.GetNDof(), felhcurl.GetNDof());

      curlmat = 0.;
      mixedmass = 0.;
      lumpedmass = 0.;

      Matrix shapehcurl(felhcurl.GetNDof(), 2);        
      Vector shapeh1(felh1.GetNDof());

      // curlmat
      for (auto i : Range(line_ends))
      {
        if (line_ends[i][0]<num)
        {
          IntegrationPoint p(pnts[line_ends[i][0]][0], pnts[line_ends[i][0]][1],0);
          felh1.CalcShape(p, shapeh1);    
          curlmat.Row(i) += -shapeh1;
        }
        if (line_ends[i][1]<num)
        {
          IntegrationPoint q(pnts[line_ends[i][1]][0], pnts[line_ends[i][1]][1],0);
          felh1.CalcShape(q, shapeh1);   
          curlmat.Row(i) += shapeh1;
        }

      }
      if (printmats)
        cout << "curlmat = " << curlmat << endl;
      //throw Exception("test");

      // mixedmass
      for (auto i : Range(lines))
      {
        //cout << i << endl;
        for (auto e : lines[i])
        {
          //cout << e[0] << endl;
          Vec<2> micropts[2] =
          { pnts[e[0]], pnts[e[1]]};
          //cout << micropts[1] << micropts[0]<< endl;
          Vec<2> tau = micropts[1]-micropts[0];
          Vec<2> taurot = Vec<2>{tau[1],-tau[0]};

          for (auto ip : irsegm)
          {
            Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
            IntegrationPoint ipseg(p(0), p(1), 0);
            felhcurl.CalcShape(ipseg, shapehcurl,mappingD);    
            mixedmass.Row(i) += ip.Weight()*(shapehcurl*taurot);
            //mixedmass.Row(i) += ip.Weight()*(shapehcurl*tau);
          }
        }
      }
      if (printmats)
        cout << "mixedmass = " << mixedmass << endl;
      CalcInverse (mixedmass);

      //cout << "mixedmassinv = " << mixedmass << endl;
      // integrate
      // lumpedmass ...
      //lumped integration is the same (except order=0 fixup)
      //auto irs = feshcurl->GetIntegrationRules();
      //auto irs = GetIntegrationRules(2*order+6);
      //IntegrationRule ir = std::move(irs[ET_TRIG]);

      //auto irs = ngcomp::GetIntegrationRules(2*order+6);
      auto irs = feshcurl->GetIntegrationRules();
      IntegrationRule ir = std::move(irs[ET_TRIG]);            

      Matrix shapes(felhcurl.GetNDof(), ir.Size()*2);
      Matrix Piolashapes(felhcurl.GetNDof(), ir.Size()*2);

      for (int i = 0; i < ir.Size(); i++)
      {
        felhcurl.CalcShape (ir[i], shapes.Cols(2*i, 2*i+2),mappingE);
        felhcurl.CalcShape (ir[i], Piolashapes.Cols(2*i, 2*i+2),mappingD);
        Piolashapes.Cols(2*i, 2*i+2) *= ir[i].Weight();
      }

      lumpedmass = shapes * Trans (Piolashapes);

      if (printmats) 
        cout << "lumpedmass = " << lumpedmass << endl;

      Matrix<> tmpmat(felhcurl.GetNDof(), felh1.GetNDof());
      tmpmat = mixedmass * curlmat;

      if (dual)
        curlmat = lumpedmass * tmpmat;
      else
        curlmat = tmpmat;
      //cout << "finalmat = " << finalmat << endl;

      Table<DofId> xdofs(elclass_inds.Size(), felh1.GetNDof()),
        ydofs(elclass_inds.Size(), felhcurl.GetNDof());

      Array<DofId> dnumsx, dnumsy;
      for (auto i : Range(elclass_inds))
      {
        ElementId ei(VOL, elclass_inds[i]);
        feshcurl->GetDofNrs(ei, dnumsy);
        GetDofNrs(ei, dnumsx);
        xdofs[i] = dnumsx;
        ydofs[i] = dnumsy;
      }
      auto mat = make_shared<T_ConstEBE>
        ( feshcurl->GetNDof(),GetNDof(),
          curlmat, std::move(ydofs),std::move(xdofs));
      if (sum)
        sum = make_shared<SumMatrix>(sum, mat);
      else
        sum = mat;
    }
    return sum;

  }


  shared_ptr<BaseMatrix> H1PrimalCells ::
    GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
    {
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

            auto & felref = dynamic_cast<const H1PrimalCellTrig&> (GetFE(ElementId(VOL,0), lh));
            Matrix shapes(felref.GetNDof(), ir.Size());
            Matrix shapes_trans(ir.Size(), felref.GetNDof());

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
              auto & felref = dynamic_cast<const H1PrimalCellTrig&> (GetFE(ElementId(VOL,nr), lh));

              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Col(i));
              for (int i = 0; i < shapes.Height(); i++)
                for (int j = 0; j < shapes.Width(); j++)
                  if (fabs(shapes(i,j)) < 1e-14)
                    shapes(i,j) = 0;
              SmallSparseMatrix spshapes(shapes, 1e-14);
              shapes_trans = Trans(shapes);

              tint.Start();

              MyMutex COOmutex;
              ParallelForRange (elclass_inds.Size(), [&](auto myrange)
                  {
                  auto mylh = lh.Split();
                  Array<int> myrowind, mycolind;
                  Array<double> myvalues;
                  Array<DofId> dofs(felref.GetNDof());
                  Matrix rhoi_shapes_trans(ir.Size(), felref.GetNDof());
                  Matrix<> elmat(felref.GetNDof());

                  for (auto i : myrange)
                  {
                  HeapReset hr(mylh);
                  auto nr = elclass_inds[i];

                  auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                  MappedIntegrationRule<2,2> mir(ir, trafo, mylh);
                  GetDofNrs(ElementId(VOL,nr), dofs);

                  for (size_t i = 0; i < mir.Size(); i++)
                  {
                  double rhoi = 1.;  //coefficient is always 1 right now
                  rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                  rhoi_shapes_trans.Row(i) = rhoi * shapes_trans.Row(i);
                  }

                  spshapes.Mult (rhoi_shapes_trans, elmat);

                  for (int i = 0; i < dofs.Size(); i++)
                    for (int j = 0; j < dofs.Size(); j++)
                      if (fabs(elmat(i,j)) > 1e-14)
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

            return make_shared<MySuperSparseMatrix> (std::move(*spmat));
          }
        case 3:
          {
            auto irs = this->GetIntegrationRules();
            IntegrationRule ir = std::move(irs[ET_TET]);

            auto & felref = dynamic_cast<const H1PrimalCellTet&> (GetFE(ElementId(VOL,0), lh));
            Matrix shapes(felref.GetNDof(), ir.Size());
            Matrix shapes_trans(ir.Size(), felref.GetNDof());

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
              auto & felref = dynamic_cast<const H1PrimalCellTet&> (GetFE(ElementId(VOL,nr), lh));

              for (int i = 0; i < ir.Size(); i++)
                felref.CalcShape (ir[i], shapes.Col(i));
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
                  Matrix rhoi_shapes_trans(ir.Size(), felref.GetNDof());
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
                  double rhoi = 1.;  //coefficient is always 1 right now
                  rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
                  rhoi_shapes_trans.Row(i) = rhoi * shapes_trans.Row(i);
                  }

                  spshapes.Mult (rhoi_shapes_trans, elmat);

                  for (int i = 0; i < dofs.Size(); i++)
                    for (int j = 0; j < dofs.Size(); j++)
                      if (fabs(elmat(i,j)) > 1e-14)
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

            return make_shared<MySuperSparseMatrix> (std::move(*spmat));
          }
        default:
          throw Exception("MassOperator not implemented for dim!=2,3");
      }
    }




  std::map<ELEMENT_TYPE, IntegrationRule> H1PrimalCells::GetIntegrationRules(bool fix_lo) const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;

    IntegrationRule irtrig;
    for (int v = 0; v < 3; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
        {
          Vec<2> xi(ipx(0), ipy(0));
          Vec<2> x = MapQuad2Trig (xi, v);
          Mat<2,2> jac = DMapQuad2Trig(xi, v);
          double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
          irtrig += IntegrationPoint(x(0), x(1), 0, w);
        }

    if (GaussRadauIR.Size()==1 && fix_lo)
      for (auto & ip : irtrig)
        ip.SetWeight(1.0/6);   // fix consistency for lowest order


    if (GaussRadauIR.Size()==-2)
    { // problem: get negative weights

      // A .. corner, B .. on edge, C .. in cell
      Array<int> classify = { 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2 };
      auto f0 = [](double x, double y) { return 1 ; };
      auto f1 = [](double x, double y)
      {
        return x*x ;
        /*
           double lams[3] = { x, y, 1-x-y };
           if (lams[0] > lams[1] && lams[0] > lams[2])
           {
           return x*y;
           }
           return 0.0;
           */
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
      cout << "mat = " << endl << mat << endl;
      cout << "Inv(mat) = " << endl << Inv(mat) << endl;
      auto irs = ngcomp::GetIntegrationRules(5);
      IntegrationRule irref = std::move(irs[ET_TRIG]);
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

      cout.precision(16);
      Vec<3> sol = Inv(mat) * rhs;
      cout << "sol = " << sol << endl;
      cout << "ir = " << irtrig << endl;

      Vec<3> rsol { irtrig[0].Weight(), irtrig[1].Weight(), irtrig[3].Weight() };
      cout << "residual = " << mat*rsol << "-" << rhs << " = " << mat*rsol-rhs << endl;

      for (auto i : Range(irtrig))
        irtrig[i].SetWeight(sol(classify[i]));
    }



    rules[ET_TRIG] = std::move(irtrig);

    IntegrationRule irtet;
    for (int v = 0; v < 4; v++)
      for (auto ipx : GaussRadauIR)
        for (auto ipy : GaussRadauIR)
          for (auto ipz : GaussRadauIR)
          {
            Vec<3> xi(ipx(0), ipy(0), ipz(0));
            Vec<3> x = MapHex2Tet (xi, v);
            Mat<3,3> jac = DMapHex2Tet(xi, v);
            double w = ipx.Weight()*ipy.Weight()*ipz.Weight()*fabs(Det(jac));
            irtet += IntegrationPoint(x(0), x(1), x(2), w);
          }

    if (GaussRadauIR.Size()==1 && fix_lo)
      for (auto & ip : irtet)
        ip.SetWeight(1.0/24);   // fix consistency for lowest order





    rules[ET_TET] = std::move(irtet);

    IntegrationRule irseg;
    for (auto & ip : GaussRadauIR)
    {
      double x = ip(0), w = ip.Weight();
      irseg.Append ( IntegrationPoint(x/2, 0, 0, 0.5*w));
      irseg.Append ( IntegrationPoint((x+1)/2, 0, 0, 0.5*w));
    }
    rules[ET_SEGM] = std::move(irseg);

    return rules;
  }    
}
