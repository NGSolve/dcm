#include <comp.hpp>
#include "hdivdualcells.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"
#include "intrules.hpp"


namespace ngcomp
{
  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;
  
  class HDivDualCellTrig : public HDivCellFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & IR;
  public:
    HDivDualCellTrig (const IntegrationRule & _IR)
      : HDivCellFiniteElement<2> (3*2*_IR.Size()*_IR.Size(), _IR.Size()-1),      
      IR(_IR)
    { ; }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    {
      //cout << "HDivDualCellTrig.CalcShape called" << endl ;

      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.AddSize(ndof, 2) = 0;
      //shape = 0;
      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;

      if (vnums[minv]>vnums[maxv])
      {
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

      Mat<2> trafo = 1.0/Det(F2*F)*(F2*F);


      int nd = IR.Size();
      ArrayMem<double, 20> polxi(nd), poleta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);

      auto assign_shape = [&](int nr, int ix, int iy,int dir)
      {
        switch (dir)
        {
          case 0:
            {
              shape.Row(nr) = trafo*Vec<2>(polxi[ix]*poleta[iy],0);
              break;
            }
          case 1:
            {
              shape.Row(nr) = trafo*Vec<2>(0,-polxi[ix]*poleta[iy]);
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
              for (int k = 0; k < nd; k++)
                assign_shape(ii+k, k, 0, 1);
            }
            else
            {
              for (int k = 0; k < nd; k++)
                assign_shape(ii+k, 0, k, 0);
            }
          }
          ii+=nd;
        }
      }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
      {

        if (f[i] == maxlam)
        {
          int v1 = f[(i+1)%3];
          int kk = ii;
          for (int l = 1; l < nd; l++)
          {
            for (int k = 0; k < nd; k++)
            {
              if (v1 == minv)
              {
                assign_shape(kk++, k, l,1);
                assign_shape(kk++, l, k,0);
              }
              else 
              {
                assign_shape(kk++, l, k,0);
                assign_shape(kk++, k, l,1);
              }

            }
          }
        }
        ii += 2*(nd-1)*nd;
      }
      //cout << "CalcShape done, shape = " << shape << endl;
    }

    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> divshape) const override
    {  

      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1)};
      int maxlam = PosMax(lam);

      //divshape = 0;

      divshape.AddSize(ndof, 1) = 0;

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

      double trafo = 1.0/Det(F2*F); 

      int nd = IR.Size();

      ArrayMem<AutoDiff<1>, 20> polxi(nd), poleta(nd);

      LagrangePolynomials(AutoDiff<1>(xi(0),0), IR, polxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), IR, poleta);


      auto assign =  [&](int nr, int ix, int iy, int dir)
      {
        switch (dir)
        {
          case 0:
            {
              divshape.Row(nr) = trafo*Vec<1>(polxi[ix].DValue(0)*poleta[iy].Value());
              break;
            }
          case 1:
            {
              divshape.Row(nr) = trafo*Vec<1>(-polxi[ix].Value()*poleta[iy].DValue(0));
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
              for (int k = 0; k < nd; k++)
                assign(ii+k, k, 0, 1);
            }
            else
            {
              for (int k = 0; k < nd; k++)
                assign(ii+k, 0, k, 0);
            }
          }
          ii+=nd;
        }
      }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
      {

        if (f[i] == maxlam)
        {
          int v1 = f[(i+1)%3];
          int kk = ii;
          for (int l = 1; l < nd; l++)
          {
            for (int k = 0; k < nd; k++)
            {
              if (v1 == minv)
              {
                assign(kk++, k, l,1);
                assign(kk++, l, k,0);
              }
              else 
              {
                assign(kk++, l, k,0);
                assign(kk++, k, l,1);
              }

            }
          }
        }
        ii += 2*(nd-1)*nd;
      }
    }

    virtual void CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    { 
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = PosMax(lam);

      shape.AddSize(ndof, 2) = 0;
      //shape = 0;
      int minv = (maxlam+1)%3;
      int maxv = (maxlam+2)%3;

      if (vnums[minv]>vnums[maxv])
      {
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

      Mat<2> trafo = 1.0/Det(F2*F)*(F2*F);


      int nd = IR.Size();
      ArrayMem<double, 20> polxi(nd), poleta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);


      auto assign_shape = [&](int nr, int ix, int iy,int dir)
      {
        Vec<2> xinode(IR[ix](0), IR[iy](0));

        Mat<2,2> F = DMapQuad2Trig(xinode);

        Mat<2> trafo = 1.0/Det(F2*F)*(F2*F);
        switch (dir)
        {
          case 0:
            {
              shape.Row(nr) = trafo*Vec<2>(polxi[ix]*poleta[iy],0);
              break;
            }
          case 1:
            {
              shape.Row(nr) = trafo*Vec<2>(0,-polxi[ix]*poleta[iy]);
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
              for (int k = 0; k < nd; k++)
                assign_shape(ii+k, k, 0, 1);
            }
            else
            {
              for (int k = 0; k < nd; k++)
                assign_shape(ii+k, 0, k, 0);
            }
          }
          ii+=nd;
        }
      }
      IVec<4> f = GetVertexOrientedFace(0);
      for (int i = 0; i < 3; i++)
      {

        if (f[i] == maxlam)
        {
          int v1 = f[(i+1)%3];
          int kk = ii;
          for (int l = 1; l < nd; l++)
          {
            for (int k = 0; k < nd; k++)
            {
              if (v1 == minv)
              {
                assign_shape(kk++, k, l,1);
                assign_shape(kk++, l, k,0);
              }
              else 
              {
                assign_shape(kk++, l, k,0);
                assign_shape(kk++, k, l,1);
              }

            }
          }
        }
        ii += 2*(nd-1)*nd;
      }
    }

    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    { ; }
  };


  class HDivDualCellTet : public HDivCellFiniteElement<3>, public VertexOrientedFE<ET_TET>
  {
    const IntegrationRule & IR;
  public:
    HDivDualCellTet (const IntegrationRule & _IR)
      : HDivCellFiniteElement<3> (4*3*_IR.Size()*_IR.Size()*_IR.Size()-6*_IR.Size()*_IR.Size(), _IR.Size()-1),
      IR(_IR)
    { ; }
    using VertexOrientedFE<ET_TET>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      shape.AddSize(ndof, 3) = 0;
      //shape = 0;

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

      Mat<3> trafo = 1./Det(F2*F)*(F2*F);

      int nd = IR.Size();


      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);


      auto assign =  [&](int nr, IVec<3> i, int dir, bool flipsign = false)
      {
        double sig = 1;
        if (flipsign)
          sig = -1;
        switch (dir)
        {
          case 0:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
          case 1:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0);
              break;
            }
          case 2:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]]);
              break;
            }
            break;
        }
      };




      int ii = 0;
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
            for (int l = 0; l < nd; l++)
              for (int k = 0; k < nd; k++)
              {
                ind[dirv1] = k;
                ind[dirv2] = l;
                assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
              }
          }
          ii += nd*nd;
        }
      }

      ii += maxlam*3*sqr(nd)*(nd-1);
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
          {
            assign(ii++, { k, i, j }, 0);
            assign(ii++, { i, k, j }, 1);
            assign(ii++, { i, j, k }, 2);
          }
    }

    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> divshape) const override
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //divshape = 0;
      divshape.AddSize(ndof, 3) = 0;

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
      
      double trafo = 1.0/fabs(Det(F2*F));

      int nd = IR.Size();

      ArrayMem<AutoDiff<1>, 20> Dpolxi(nd), Dpoleta(nd), Dpolzeta(nd);
      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);         
      LagrangePolynomials(AutoDiff<1>(xi(0),0), IR, Dpolxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), IR, Dpoleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), IR, Dpolzeta);

      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);

      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          switch (dir)
            {
            case 0:
            {
              divshape(nr) = trafo*Dpolxi[i[0]].DValue(0)*poleta[i[1]]*polzeta[i[2]];
              break;
            }
            case 1:
            {
              divshape(nr) = trafo*polxi[i[0]]*Dpoleta[i[1]].DValue(0)*polzeta[i[2]];
              break;
            }
            case 2:
            {
              divshape(nr) = trafo*polxi[i[0]]*poleta[i[1]]*Dpolzeta[i[2]].DValue(0);
              break;
            }
            break;
            }
          if (flipsign)
            divshape(nr)*=-1.;
        };

      int ii = 0;
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
            for (int l = 0; l < nd; l++)
              for (int k = 0; k < nd; k++)
              {
                ind[dirv1] = k;
                ind[dirv2] = l;
                assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
              }
          }
          ii += nd*nd;
        }
      }

      ii += maxlam*3*sqr(nd)*(nd-1);
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
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

      shape.AddSize(ndof, 3) = 0;
      //shape = 0;

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

      Mat<3> trafo = 1./Det(F2*F)*(F2*F);

      int nd = IR.Size();


      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);


      auto assign =  [&](int nr, IVec<3> i, int dir, bool flipsign = false)
      {
        Vec<3> xinode(IR[i[0]](0), IR[i[1]](0), IR[i[2]](0));
        Mat<3,3> Fnode = DMapHex2Tet(xinode);          
        Mat<3> trafonode = 1.0/fabs(Det(F2*Fnode))*(F2*Fnode);
        if (flipsign)
          trafonode *= -1;
        switch (dir)
        {
          case 0:
            {
              shape.Row(nr) = trafonode*Vec<3>(polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
          case 1:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0);
              break;
            }
          case 2:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]]);
              break;
            }
            break;
        }
      };




      int ii = 0;
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
            for (int l = 0; l < nd; l++)
              for (int k = 0; k < nd; k++)
              {
                ind[dirv1] = k;
                ind[dirv2] = l;
                assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
              }
          }
          ii += nd*nd;
        }
      }

      ii += maxlam*3*sqr(nd)*(nd-1);
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
          {
            assign(ii++, { k, i, j }, 0);
            assign(ii++, { i, k, j }, 1);
            assign(ii++, { i, j, k }, 2);
          }
    }




    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    { ; }
  };




  HDivDualCells::
  HDivDualCells (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)      
  {

    if (ma->GetDimension()==2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<2>>> ());
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<3>>> ());
      }
      
    
    
    Array<double> xi, wi;
    ComputeGaussRadauRule (order+1, xi, wi);
    for (auto i : Range(xi))
      {
        IR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      }
  }


  void HDivDualCells :: Update()
  {
    //cout << "HDivDualCells.Update called" << endl;
    FESpace::Update();
    size_t ndof = 0;

    int nd = IR.Size();

    if (ma->GetDimension() == 2)
    {
      first_edge_dofs.SetSize(ma->GetNEdges()+1);
      first_edge_dofs[0] = ndof;
      for (auto i : Range(ma->GetNEdges()))
        first_edge_dofs[i+1] = ndof += 2*nd;


      first_face_dofs.SetSize(ma->GetNFaces()+1);
      first_face_dofs[0] = ndof;
      for (auto i : Range(ma->GetNFaces()))
        first_face_dofs[i+1] = ndof += 3*2*(nd-1)*nd;
    }
    if (ma->GetDimension() == 3)
    {
      first_face_dofs.SetSize(ma->GetNFaces()+1);
      first_face_dofs[0] = ndof;
      for (auto i : Range(ma->GetNFaces()))
        first_face_dofs[i+1] = ndof += 3*nd*nd;

      first_cell_dofs.SetSize(ma->GetNE(VOL)+1);
      first_cell_dofs[0] = ndof;
      for (auto i : Range(ma->GetNE(VOL)))
        first_cell_dofs[i+1] = ndof += 4*3*nd*(nd-1)*(nd-1);
    }
    SetNDof(ndof);
    //cout << "HDivDualCells.Update done" << endl;
  }

  void HDivDualCells ::
    GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
      //cout << "HDivDualCells.GetDofNrs called for ElementId " << ei << endl;
      dnums.SetSize0();
      auto el = ma->GetElement(ei);
      if (ma->GetDimension() == 2)
      {
        for (auto ed : el.Edges())
          dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
        for (auto fa : el.Faces())
          dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
      }
      if (ma->GetDimension() == 3)
      {
        for (auto fa : el.Faces())
          dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
        if (ei.VB() == VOL)
          dnums += Range(first_cell_dofs[ei.Nr()], first_cell_dofs[ei.Nr()+1]);
      }
      //cout << "HDivDualCells.GetDofNrs done" << endl;
    }


  FiniteElement & HDivDualCells ::
  GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TET:
        {
          auto tet = new (alloc) HDivDualCellTet(IR);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      case ET_TRIG:
        {
          auto trig = new (alloc) HDivDualCellTrig(IR);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      default:
        throw Exception("element not implemented");
      }
  }



  shared_ptr<BaseMatrix> HDivDualCells ::
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
          IntegrationRule ir = move(irs[ET_TRIG]);
          auto & felref = dynamic_cast<const HDivDualCellTrig&> (GetFE(ElementId(VOL,0), lh));
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
              auto & felref = dynamic_cast<const HDivDualCellTrig&> (GetFE(ElementId(VOL,nr), lh));
          
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
                      
                                          Mat<2,2> F = mir[i].GetJacobian();
                                          rhoi = Trans(F) * rhoi * F;
                                          rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
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
          auto & felref = dynamic_cast<const HDivDualCellTet&> (GetFE(ElementId(VOL,0), lh));
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
              auto & felref = dynamic_cast<const HDivDualCellTet&> (GetFE(ElementId(VOL,nr), lh));
          
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
                      
                                          Mat<3,3> F = mir[i].GetJacobian();
                                          rhoi = Trans(F) * rhoi * F;
                                          rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
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

  }

  

  
  std::map<ELEMENT_TYPE, IntegrationRule> HDivDualCells::GetIntegrationRules(bool fix_lo) const
  {
      std::map<ELEMENT_TYPE, IntegrationRule> rules;

      rules[ET_TRIG] = PrimalCellIR(IR);
      rules[ET_TET] = PrimalVolIR(IR);
      rules[ET_SEGM] = PrimalSegmIR(IR);
      if (IR.Size()==1 && fix_lo)
      {
        for (auto & ip : rules[ET_TRIG])
          ip.SetWeight(1.0/6);
        for (auto & ip : rules[ET_TET])
          ip.SetWeight(1.0/24);
      }

      return rules;
  }

};
