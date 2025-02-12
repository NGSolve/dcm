#include <comp.hpp>
#include "dcs_common.hpp"

namespace ngcomp
{


  void ComputeSortedGaussRadauPoints(int order,Array<double> &xi,bool ascending)
  {
    Array<double> wi,xitmp;
    xi.SetSize(0);
    ComputeGaussRadauRule(order,xitmp,wi);
    BubbleSort(xitmp);
    if (ascending)
      for (auto i : Range(xitmp))
        xi.Append(xitmp[i]);
    else
      for (auto i : Range(xitmp))
        xi.Append(xitmp[order-i-1]);
  }


  template <typename T>
  void LagrangePolynomials (T x, const IntegrationRule & ir, FlatArray<T> shapes)
  {
    for (int i = 0; i < shapes.Size(); i++)
      {
        T prod = 1;
        for (int j = 0; j < ir.Size(); j++)
          if (j != i)
            prod *= (x-ir[j](0)) / (ir[i](0)-ir[j](0));
        shapes[i] = prod;
      }
  }

  template void LagrangePolynomials<double> (double x, const IntegrationRule & ir, FlatArray<double> shapes);
  template void LagrangePolynomials<AutoDiff<1>> (AutoDiff<1> x, const IntegrationRule & ir, FlatArray<AutoDiff<1>> shapes);
  template void LagrangePolynomials<AutoDiff<2>> (AutoDiff<2> x, const IntegrationRule & ir, FlatArray<AutoDiff<2>> shapes);  

  
  IntegrationRule PrimalVolIR(const IntegrationRule & IR, bool mirroredIR)
  {
        IntegrationRule ircube;
        for (auto ipx : IR)
          for (auto ipy : IR)
            for (auto ipz : IR)
          {
            if (mirroredIR)
              ircube += IntegrationPoint(1-ipx(0), 1-ipy(0), 1-ipz(0), ipx.Weight()*ipy.Weight()*ipz.Weight());
            else 
              ircube += IntegrationPoint(ipx(0), ipy(0), ipz(0), ipx.Weight()*ipy.Weight()*ipz.Weight());
          }

        IntegrationRule ir;
        ir.SetDim(3);
        for (int i = 0; i < 4; i++)
          {
            for (auto & ip : ircube)
              {
                Vec<3> p = MapHex2Tet(Vec<3> ({ip(0),ip(1),ip(2)}), i);
                Mat<3> F = DMapHex2Tet(Vec<3> ({ip(0),ip(1),ip(2)}), i);
                ir.Append (IntegrationPoint( p(0), p(1), p(2), ip.Weight()*fabs(Det(F))));
              }
          }
    return ir;
  }

  IntegrationRule PrimalCellIR(const IntegrationRule & IR, bool mirroredIR)
  {
        IntegrationRule irquad;
        for (auto ipx : IR)
          for (auto ipy : IR)
          {
            if (mirroredIR)
              irquad += IntegrationPoint(1-ipx(0)-1e-12, 1-ipy(0)-1e-12, 0, ipx.Weight()*ipy.Weight());
            else 
            irquad += IntegrationPoint(ipx(0), ipy(0), 0, ipx.Weight()*ipy.Weight());
          }

        IntegrationRule ir;
        Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };
        for (int i = 0; i < 3; i++)
          {
            Vec<2> p0 = verts[i];
            Vec<2> p1 = 0.5 * (verts[i]+verts[(i+1)%3]);
            Vec<2> p2 = { 1.0/3, 1.0/3 };
            Vec<2> p3 = 0.5 * (verts[i]+verts[(i+2)%3]);

            for (auto & ip : irquad)
              {
                double x = ip(0), y = ip(1);
                Vec<2> p = (1-x)*(1-y) * p0 + x*(1-y) * p1 + x*y*p2+(1-x)*y*p3;
                Vec<2> dpdx = y*(p2-p3)+(1-y)*(p1-p0);
                Vec<2> dpdy = x*(p2-p1)+(1-x)*(p3-p0);
                Mat<2> F;
                F.Col(0) = dpdx;
                F.Col(1) = dpdy;
                ir.Append (IntegrationPoint( p(0), p(1), 0, ip.Weight()*Det(F)));
              }
          }
    return ir;
  }

  IntegrationRule PrimalSegmIR(const IntegrationRule & irseg)
  {
        IntegrationRule ir;
        for (auto & ip : irseg)
          {
            double x = ip(0), w = ip.Weight();
            ir.Append ( IntegrationPoint(x/2, 0, 0, 0.5*w));
            ir.Append ( IntegrationPoint((x+1)/2, 0, 0, 0.5*w));
          }
        return ir;


  }



  std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(int order) 
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    IntegrationRule ir(ET_SEGM, order);

    IntegrationRule irsegm;
    for (int v = 0; v < 2; v++)
      for (auto ipx : ir)
        {
          double xi = ipx(0);
          double x = v==0 ? 0.5*xi : 1-0.5*xi;
          double w = 0.5*ipx.Weight();
          irsegm += IntegrationPoint(x, 0, 0, w);
        }
    
    rules[ET_SEGM] = std::move(irsegm);

    

    IntegrationRule irtrig;
    for (int v = 0; v < 3; v++)
      for (auto ipx : ir)
        for (auto ipy : ir)
            {
              Vec<2> xi(ipx(0), ipy(0));
              Vec<2> x = MapQuad2Trig (xi, v);
              Mat<2,2> jac = DMapQuad2Trig (xi, v);
              double w = ipx.Weight()*ipy.Weight()*fabs(Det(jac));
              irtrig += IntegrationPoint(x(0), x(1), 0, w);
            }
    
    rules[ET_TRIG] = std::move(irtrig);

    
    IntegrationRule irtet;
    for (int v = 0; v < 4; v++)
      for (auto ipx : ir)
        for (auto ipy : ir)
          for (auto ipz : ir)
            {
              Vec<3> xi(ipx(0), ipy(0), ipz(0));
              Vec<3> x = MapHex2Tet (xi, v);
              Mat<3,3> jac = DMapHex2Tet(xi, v);
              double w = ipx.Weight()*ipy.Weight()*ipz.Weight()*fabs(Det(jac));
              irtet += IntegrationPoint(x(0), x(1), x(2), w);
            }
    
    rules[ET_TET] = std::move(irtet);

    return rules;
  }



  std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRulesInnerFacets(int intorder)
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;

    IntegrationRule ir(ET_SEGM, intorder);

    IntegrationRule irtrig;
    for (int i = 0; i < 3; i++)
      {
        for (auto ip : ir)
          {
            AutoDiff<1> s(ip(0),0);
            AutoDiff<1> lami[3];
            lami[i] = 1.0/3 * s;
            lami[(i+1)%3] = 0.5 - 1.0/6 * s;
            lami[(i+2)%3] = 0.5 - 1.0/6 * s;

            AutoDiff<1> x = lami[0];
            AutoDiff<1> y = lami[1];

            // irtrig += IntegrationPoint(x.Value(), y.Value(), 0, sqrt(sqr(x.DValue(0))+sqr(y.DValue(0))));
            irtrig += IntegrationPoint(x.Value(), y.Value(), 0, ip.Weight()/6);
          }
      }
    rules[ET_TRIG] = std::move(irtrig);

    return rules;
  }


  
  // i ... 2D coordinates inm corner
  // returns 3D coordinates on manifold
  //
  //inline constexpr array<int,3> C2M (int corner, IVec<2> i, int order)
  inline array<int,3> C2M (int corner, IVec<2> i, int order)
  {
    if (corner == 0) return { order, i[0], i[1] };
    if (corner == 1) return { i[1], order, i[0] };
    return { i[0], i[1], order };
  }


  template<>
  tuple<Array<array<int,2>>,int> GetMicroCellEdges<ET_TRIG> (int order, bool dual)
  {
    Array<array<int,2>> edges;
    
    auto [pnts, _inds, np] = GetMicroCellPoints<ET_TRIG> (order, dual);
    auto inds = std::move(_inds);
    
    auto addedge = [&] (int corner, IVec<2> p0, IVec<2> p1)
      {
        int i0 = inds.Pos (C2M (corner, p0, order));
        int i1 = inds.Pos (C2M (corner, p1, order));
        edges.Append ( array<int,2>{min(i0, i1),max(i0,i1)});
      };
    
    for (int c = 0; c < 3; c++)
      for (int k = 0; k < order; k++)
        for (int l = 0; l < order; l++)
          {
            addedge ( c, { k, l }, { k+1,l } );
            addedge ( c, { k, l }, { k,l+1 } );
          }
    
    int nreal = edges.Size();
    
    for (int c = 0; c < 3; c++)
      for (int k = 0; k < order; k++)
        addedge ( c, { k, order }, { k+1, order } );


    return make_tuple(edges, nreal);
  } 


  template<>
  Table<array<int,4>> GetMicroCellFaces<ET_TRIG> (int order, bool dual)
  {
    // cout << "Called GetMicroCellFaces" << endl;
    auto [pnts, _inds, np] = GetMicroCellPoints<ET_TRIG> (order, dual);
    auto inds = std::move(_inds);
    TableCreator<array<int,4>> creator;

    auto addface = [&] (int i, int corner, IVec<2> p0, IVec<2> p2)
      {
        creator.Add (i,
                     array<int,4> {
                       int(inds.Pos (C2M ( corner, p0, order ))),
                         int(inds.Pos (C2M ( corner, { p0[0], p2[1] }, order ))),
                         int(inds.Pos (C2M ( corner, p2, order ))),
                         int(inds.Pos (C2M ( corner, { p2[0], p0[1] }, order ))),
                         });
      };
    
    for ( ; !creator.Done(); creator++)
      {
        int i = 0;
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < order-1; k++)
              for (int l = 0; l < order-1; l++)
                {
                  // creator.Add(i,addface (c, {k,l},{k+1,l+1}));
                  addface (i, c, {k,l},{k+1,l+1});
                  i++;
                }
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < order-1; k++)
              {
                addface (i, c, {order-1,k},{order,k+1});
                addface (i, (c+1)%3, {k,order-1},{k+1,order});
                i++;
              }
          for (int c = 0; c < 3; c++)
            addface (i,c, {order-1,order-1},{order,order});
      }
    Table<array<int,4>> faces_table = creator.MoveTable();
    
    return faces_table;
  }


  template<>
  Table<array<int,2>> GetMicroCellFaceBoundaries<ET_TRIG> (int order, bool dual)
  {
    Table<array<int,4>> faces_tab = GetMicroCellFaces<ET_TRIG> (order, dual);
    TableCreator<array<int,2>> creator;
    
    for ( ; !creator.Done(); creator++)
      {
        for (auto i : Range(faces_tab))
          {
            Array<array<int,2>> tmparray;
            auto faces = faces_tab[i];
            for (auto k : Range(faces))
              for (int j = 0; j < 4; j++)
                {
                  array<int,2> tmpedge = { faces[k][j],faces[k][(j+1)%4] };
                  array<int,2> tmpedge2 = { tmpedge[1], tmpedge[0] }; 
                  if (!tmparray.Contains(tmpedge2))
                    tmparray.Append(tmpedge);
                  else
                    tmparray.DeleteElement(tmparray.Pos(tmpedge2)); 
                }
            
            for (auto ed : tmparray)
              creator.Add(i,ed);
          }
      }
    
    return creator.MoveTable();
  }






  // i ... 2D coordinates inm corner
  // returns 3D coordinates on manifold
  //
  //inline constexpr array<int,4> C2M (int corner, IVec<3> i, int order)
  inline array<int,4> C2M3D (int corner, IVec<3> i, int order)
  {
    if (corner == 0) return { order, i[0], i[1], i[2] };
    if (corner == 1) return { i[2], order, i[0], i[1] };
    if (corner == 2) return { i[1], i[2], order, i[0] };
                     return { i[0], i[1], i[2], order };
  }

  template<>
  tuple<Array<Vec<3>>,Array<array<int,4>>,int> GetMicroCellPoints<ET_TET> (int order, bool dual)
  {  // order x order squares per corner
    Array<double> xi;
    IntegrationRule GaussRadauIRplus;
    ComputeSortedGaussRadauPoints (order, xi, false);
    for (auto i : Range(xi))
    {
      GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, 0));
    }
    GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));

    Array<Vec<3>> pnts;
    Array<array<int,4>> inds;
    
    auto addpoint = [&] (int corner, IVec<3> i)
    {
      pnts.Append ( MapHex2Tet ( Vec<3> (GaussRadauIRplus[i[0]](0),
                                         GaussRadauIRplus[i[1]](0),
                                         GaussRadauIRplus[i[2]](0)),
                                 corner) );
      inds.Append (C2M3D(corner, i, order));
    };
    
    for (int c = 0; c < 4; c++)
      for (int k = 0; k < order; k++)
        for (int l = 0; l < order; l++)
          for (int m = 0; m < order; m++)
            addpoint (c, { k, l, m } );

    for (int k = 0; k < order; k++)
      for (int l = 0; l < order; l++)
        {
          addpoint (0, { k, l, order } );
          addpoint (0, { order, k, l } );
          addpoint (0, { l, order, k } );

          addpoint (1, { order, k, l } );
          addpoint (1, { l, order, k } );

          addpoint (2, { order, k, l } );          
        }

    int nreal = pnts.Size();
    
    for (int k = 0; k < order; k++)
      {
        addpoint (0, { k, order, order } );
        addpoint (0, { order, k, order } );
        addpoint (0, { order, order, k } );
        addpoint (1, { order, order, k } );
      }

    addpoint(0, { order, order, order });
    return make_tuple(pnts, inds, nreal);
  }

  template<>
  tuple<Array<array<int,2>>,int> GetMicroCellEdges<ET_TET> (int order, bool dual)
  {
    Array<array<int,2>> edges;
    
    auto [pnts, _inds, np] = GetMicroCellPoints<ET_TET> (order, dual);
    auto inds = std::move(_inds);
    
    auto addedge = [&] (int corner, IVec<3> p0, IVec<3> p1)
      {
        int i0 = inds.Pos (C2M3D (corner, p0, order));
        int i1 = inds.Pos (C2M3D (corner, p1, order));
        edges.Append ( array<int,2>{min(i0, i1),max(i0,i1)});
      };

    auto addedgedir = [&] (int corner, IVec<3> p0, int dir)
      {
        IVec<3> p1 = p0;
        p1[dir] += 1;
        addedge (corner, p0, p1);
      };
    
    
    for (int c = 0; c < 4; c++)
      for (int k = 0; k < order; k++)
        for (int l = 0; l < order; l++)
          for (int m = 0; m < order; m++)
            {
              addedgedir ( c, { m, k, l }, 0);
              addedgedir ( c, { l, m, k }, 1);
              addedgedir ( c, { k, l, m }, 2);
            }
    
    int nreal = edges.Size();

    for (int k = 0; k < order; k++)
      for (int l = 0; l < order; l++)
        {
          addedgedir (0, { k, l, order }, 0);
          addedgedir (0, { k, l, order }, 1);
          
          addedgedir (0, { order, k, l }, 1);
          addedgedir (0, { order, k, l }, 2);
          
          addedgedir (0, { l, order, k }, 0);
          addedgedir (0, { l, order, k }, 2);

          addedgedir (1, { order, k, l }, 1);
          addedgedir (1, { order, k, l }, 2);
          
          addedgedir (1, { l, order, k }, 0);
          addedgedir (1, { l, order, k }, 2);

          addedgedir (2, { order, k, l }, 1);          
          addedgedir (2, { order, k, l }, 2);          
        }
    
    for (int k = 0; k < order; k++)
      {
        addedgedir (0, { k, order, order }, 0);
        addedgedir (0, { order, k, order }, 1);
        addedgedir (0, { order, order, k }, 2);
        addedgedir (1, { order, order, k }, 2);
      }

    return make_tuple(edges, nreal);
  } 


  template<>
  Table<array<int,4>> GetMicroCellFaces<ET_TET> (int order, bool dual)
  {
    // cout << "Called GetMicroCellFaces" << endl;
    auto [pnts, _inds, np] = GetMicroCellPoints<ET_TET> (order, dual);
    auto inds = std::move(_inds);
    TableCreator<array<int,4>> creator;

    auto addface = [&] (int ind, int corner, IVec<3> p0, int dir1, int dir2)
      {
        IVec<3> p1 = p0, p2 = p0, p3 = p0;
        p1[dir1]++;
        p2[dir1]++;
        p2[dir2]++;
        p3[dir2]++;
        creator.Add (ind,
                     array<int,4> {
                       int(inds.Pos (C2M3D ( corner, p0, order ))),
                         int(inds.Pos (C2M3D ( corner, p1, order ))),
                         int(inds.Pos (C2M3D ( corner, p2, order ))),
                         int(inds.Pos (C2M3D ( corner, p3, order ))),
                         });
      };
    
    for ( ; !creator.Done(); creator++)
      {
        int i = 0;
        for (int c = 0; c < 4; c++)
          for (int k = 0; k < order-1; k++)
            for (int l = 0; l < order-1; l++)
              for (int m = 0; m < order; m++)
                {
                  addface (i, c, {k,l,m}, 0, 1);
                  i++;                  
                  addface (i, c, {m,k,l}, 1, 2);
                  i++;
                  addface (i, c, {l,m,k}, 2, 0);
                  i++;
                }
        
        for (int k = 0; k < order-1; k++)
          for (int l = 0; l < order; l++)
            { // 6 x 4
              addface (i, 0, { order-1, k, l }, 1, 0 );
              addface (i, 1, { k, l, order-1 }, 2, 0 );
              i++;
              addface (i, 0, { order-1, l, k }, 2, 0 );
              addface (i, 1, { l, k, order-1 }, 2, 1 );
              i++;

              addface (i, 0, { k, order-1, l }, 0, 1 );
              addface (i, 2, { l, order-1, k }, 1, 2 );
              i++;
              addface (i, 0, { l, order-1, k }, 2, 1 );
              addface (i, 2, { k, order-1, l }, 1, 0 );
              i++;
              
              addface (i, 0, { k, l, order-1 }, 0, 2 );
              addface (i, 3, { order-1, k, l }, 0, 1 );
              i++;
              addface (i, 0, { l, k, order-1 }, 1, 2 );
              addface (i, 3, { order-1, l, k }, 0, 2 );
              i++;

              addface (i, 1, { order-1, k, l }, 1, 0 );
              addface (i, 2, { k, l, order-1 }, 2, 0 );
              i++;
              addface (i, 1, { order-1, l, k }, 2, 0 );
              addface (i, 2, { l, k, order-1 }, 2, 1 );
              i++;
              
              addface (i, 1, { k, order-1, l }, 0, 1 );
              addface (i, 3, { l, order-1, k }, 1, 2 );
              i++;

              addface (i, 1, { l, order-1, k }, 2, 1 );
              addface (i, 3, { k, order-1, l }, 1, 0 );
              i++;

              addface (i, 2, { order-1, k, l }, 1, 0 );
              addface (i, 3, { k, l, order-1 }, 2, 0 );
              i++;
              
              addface (i, 2, { order-1, l, k }, 2, 0 );
              addface (i, 3, { l, k, order-1 }, 2, 1 );
              i++;
              
            }

        for (int k = 0; k < order; k++)
          {
            addface (i, 0, { order-1, order-1, k }, 0, 1 );
            addface (i, 1, { order-1, k, order-1 }, 0, 2 );
            addface (i, 2, { k, order-1, order-1 }, 1, 2 );
            i++;

            addface (i, 3, { k, order-1, order-1 }, 1, 2 );
            addface (i, 1, { order-1, order-1, k }, 0, 1 );
            addface (i, 2, { order-1, k, order-1 }, 0, 2 );
            i++;

            addface (i, 3, { order-1, k, order-1 }, 0, 2 );
            addface (i, 0, { k, order-1, order-1 }, 1, 2 );
            addface (i, 2, { order-1, order-1, k }, 0, 1 );
            i++;            
            
            addface (i, 3, { order-1, order-1, k }, 0, 1 );
            addface (i, 0, { order-1, k, order-1 }, 0, 2 );
            addface (i, 1, { k, order-1, order-1 }, 1, 2 );
            i++;            
          }
      }
    Table<array<int,4>> faces_table = creator.MoveTable();
    
    return faces_table;
  }


  template<>
  Table<array<int,2>> GetMicroCellFaceBoundaries<ET_TET> (int order, bool dual)
  {
    Table<array<int,4>> faces_tab = GetMicroCellFaces<ET_TET> (order, dual);
    TableCreator<array<int,2>> creator;
    
    for ( ; !creator.Done(); creator++)
      {
        for (auto i : Range(faces_tab))
          {
            Array<array<int,2>> tmparray;
            auto faces = faces_tab[i];
            for (auto k : Range(faces))
              for (int j = 0; j < 4; j++)
                {
                  array<int,2> tmpedge = { faces[k][j],faces[k][(j+1)%4] };
                  array<int,2> tmpedge2 = { tmpedge[1], tmpedge[0] }; 
                  if (!tmparray.Contains(tmpedge2))
                    tmparray.Append(tmpedge);
                  else
                    tmparray.DeleteElement(tmparray.Pos(tmpedge2)); 
                }
            
            for (auto ed : tmparray)
              creator.Add(i,ed);
          }
      }
    
    return creator.MoveTable();
  }



  template <>
  tuple<Array<Vec<2>>,Array<array<int,3>>,int> GetNanoPoints<ET_TRIG> (int order, bool dual)
  {
    // order x order squares per corner
    Array<double> xi;
    IntegrationRule GaussRadauIRplus;
    Array<Vec<2>> pnts;
    Array<array<int,3>> inds;

    int nreal;
    if (dual)
    {
      ComputeSortedGaussRadauPoints (order, xi,true);
      for (auto i : Range(xi))
        GaussRadauIRplus.Append (IntegrationPoint (xi[i]+1e-12, 0, 0, 0));
      GaussRadauIRplus.Append (IntegrationPoint (1-1e-12, 0, 0, 0));
      auto addpoint = [&] (int corner, IVec<2> i)
      {
        pnts.Append ( MapQuad2Trig ( Vec<2> (GaussRadauIRplus[i[0]](0),
                GaussRadauIRplus[i[1]](0)), corner) );
        inds.Append (C2M(corner, i, order));
      };

      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
          for (int l = 0; l < order; l++)
            addpoint (c, { k, l } );
      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
          addpoint (c, { k, order } );

      nreal = pnts.Size();

      addpoint(0, { order, order });
    }
    else
    {
      ComputeSortedGaussRadauPoints (order, xi,false);
      GaussRadauIRplus.Append (IntegrationPoint (1e-12, 0, 0, 0));
      for (auto i : Range(xi))
        GaussRadauIRplus.Append (IntegrationPoint (1-xi[i]-1e-12, 0, 0, 0));

      auto addpoint = [&] (int corner, IVec<2> i)
      {
        pnts.Append ( MapQuad2Trig ( Vec<2> (GaussRadauIRplus[i[0]](0),
                GaussRadauIRplus[i[1]](0)), corner) );
        //inds.Append (C2M(corner, i, order));
        inds.Append ({corner, i[0],i[1]});
      };

      //points only once
      /*
      for (int c = 0; c < 3; c++)
        for (int k = 1; k < order; k++)
          for (int l = 1; l < order; l++)
            addpoint (c, { k, l } );
      for (int c = 0; c < 3; c++)
        for (int k = 1; k < order; k++)
          addpoint (c, { k, order } );

      addpoint(0, { order, order });

      nreal = pnts.Size();

      for (int c = 0; c < 3; c++)
        for (int l = 0; l < order; l++)
        {
          addpoint (c, { l, 0 } );
          addpoint (c, { 0, l } );
        }
      for (int c = 0; c < 3; c++)
          addpoint (c, { 0, order } );
          */
      //micro interface points counted twice (as number of unknowns)
      for (int c = 0; c < 3; c++)
        for (int k = 1; k < order+1; k++)
          for (int l = 1; l < order+1; l++)
            addpoint (c, { k, l } );


      nreal = pnts.Size();

      for (int c = 0; c < 3; c++)
        for (int l = 1; l < order+1; l++)
        {
          addpoint (c, { l, 0 } );
          addpoint (c, { 0, l } );
          addpoint (c, { 0, 0 } );
        }

    }
    return make_tuple(pnts, inds, nreal);
  }

  template<>
  tuple<Array<Vec<2>>,Array<array<int,3>>,int> GetMicroCellPoints<ET_TRIG> (int order, bool dual)
  {  
    return GetNanoPoints<ET_TRIG>(order,dual);
  }



  template <>
  tuple<Array<array<int,2>>,int> GetNanoEdges<ET_TRIG> (int order, bool dual)
  {
    Array<array<int,2>> edges;
    
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = std::move(_inds);
    
    auto addedge = [&] (int corner, IVec<2> p0, IVec<2> p1)
      {
        int i0 = inds.Pos (C2M (corner, p0, order));
        int i1 = inds.Pos (C2M (corner, p1, order));
        edges.Append ( array<int,2>{min(i0, i1),max(i0,i1)});
      };
    
    int nreal; 
    if (dual)
    {
      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
          for (int l = 0; l < order; l++)
            {
              addedge ( c, { k, l }, { k+1,l } );
              addedge ( c, { k, l }, { k,l+1 } );
            }
      
      nreal = edges.Size();
      
      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
          addedge ( c, { k, order }, { k+1, order } );
    }
    else
    {
      for (int c = 0; c < 3; c++)
        for (int k = 1; k < order; k++)
          for (int l = 0; l < order; l++)
          {
              addedge ( c, { k, l }, { k,l+1 } );
              addedge ( c, { l, k }, { l+1,k } );
          }
      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
          addedge ( c, { k, order }, { k+1,order } );

      nreal = edges.Size();

      for (int c = 0; c < 3; c++)
        for (int k = 0; k < order; k++)
        {
          addedge ( c, { k, 0 }, { k+1, 0 } );
          addedge ( c, { 0, k }, { 0, k+1 } );
        }
    }

    return make_tuple(edges, nreal);
  }

  template <>
  tuple<Array<array<int,4>>,int> GetNanoFaces<ET_TRIG> 
      (int order, bool dual)
  {
    throw Exception ("Nano faces not implemented for ET_TRIG");
  }

  template <>
  tuple<Array<array<int,8>>,int> GetNanoCells<ET_TRIG>
     (int order, bool dual)
  {
    throw Exception ("Nano cells not implemented for ET_TRIG");
  }


  template <>
  tuple<Table<array<int,2>>,Array<array<int,2>>> 
        GetNanoLines<ET_TRIG> (int order, bool dual)
  {
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = std::move(_inds);
    TableCreator<array<int,2>> creator;
    
    if (dual)
    {
      auto addline = [&] (int i, int corner, IVec<2> p0, IVec<2> p1)
      {
        creator.Add (i,
            array<int,2> {
            int(inds.Pos (C2M ( corner, p0, order ))),
            int(inds.Pos (C2M ( corner, p1, order ))),
            });
      };
      for ( ; !creator.Done(); creator++)
      {
        int i = 0;
        for (int c = 0; c < 3; c++)
          for (int k = 0; k < order-1; k++)
            for (int l = 0; l < order; l++)
            {
              addline (i, c, {k,l},{k+1,l});
              i++;
              addline (i, c, {l,k},{l,k+1});
              i++;
            }
        for (int c = 0; c < 3; c++)
          for (int k = 0; k < order; k++)
          {
            addline (i, c, {order-1,k},{order,k});
            addline (i, (c+1)%3, {k,order},{k,order-1});
            i++;
          }
      }
    }
    else
    {
      /* lines only once
      for ( ; !creator.Done(); creator++)
      {
        int i = 0;
        for (int c = 0; c < 3; c++)
          for (int k = 1; k < order; k++)
          {
            addline (i, c, {0,k},{1,k});
            i++;
            addline (i, c, {k,0},{k,1});
            i++;
          }
        for (int c = 0; c < 3; c++)
        {
          addline (i, c, {order,0},{order,1});
          i++;
        }

        for (int c = 0; c < 3; c++)
          for (int k = 1; k < order; k++)
            for (int l = 1; l < order; l++)
            {
              addline (i, c, {k,l},{k+1,l});
              i++;
              addline (i, c, {l,k},{l,k+1});
              i++;
            }
        for (int c = 0; c < 3; c++)
          for (int k = 1; k < order; k++)
          {
            addline (i, c, {order,k},{order,k+1});
            i++;
          }
      }
      */ 
      // lines for double dofs twice
      auto addline = [&] (int i, int corner, IVec<2> p0, IVec<2> p1)
      {
        creator.Add (i,
            array<int,2> {
            int(inds.Pos ( {corner, p0[0],p0[1] })),
            int(inds.Pos ( {corner, p1[0],p1[1] })),
            });
      };
      for ( ; !creator.Done(); creator++)
      {
        int i = 0;

        for (int c = 0; c < 3; c++)
          for (int k = 0; k < order; k++)
            for (int l = 1; l < order+1; l++)
            {
              addline (i, c, {k,l},{k+1,l});
              i++;
              addline (i, c, {l,k},{l,k+1});
              i++;
            }
      }
    }
    Table<array<int,2>> lines_table = creator.MoveTable();

    Array<array<int,2>> lines_ends;

    for (auto i : Range(lines_table))
    {

       Array<int> all_pts;
       for (auto k : Range(lines_table[i]))
           {
           array<int,2> tmpedge = lines_table[i][k];
           if(!all_pts.Contains(tmpedge[0]))
             all_pts.Append(tmpedge[0]);
           else 
             all_pts.DeleteElement(all_pts.Pos(tmpedge[0]));
           if(!all_pts.Contains(tmpedge[1]))
             all_pts.Append(tmpedge[1]);
           else 
             all_pts.DeleteElement(all_pts.Pos(tmpedge[1]));
           }
      lines_ends.Append(array<int,2> {all_pts[0],all_pts[1]});
    }
    return make_tuple(lines_table,lines_ends);
  }

  template <>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces<ET_TRIG> (int order, bool dual)
  {
    // cout << "Called GetMicroCellFaces" << endl;
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = std::move(_inds);
    TableCreator<array<int,4>> creator;

    auto addface = [&] (int i, int corner, IVec<2> p0, IVec<2> p2)
      {
        creator.Add (i,
                     array<int,4> {
                       int(inds.Pos (C2M ( corner, p0, order ))),
                         int(inds.Pos (C2M ( corner, IVec<2>{ p0[0], p2[1] }, order ))),
                         int(inds.Pos (C2M ( corner, p2, order ))),
                         int(inds.Pos (C2M ( corner, IVec<2>{ p2[0], p0[1] }, order ))),
                         });
      };
    
    for ( ; !creator.Done(); creator++)
      {
        int i = 0;
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < order-1; k++)
              for (int l = 0; l < order-1; l++)
                {
                  // creator.Add(i,addface (c, {k,l},{k+1,l+1}));
                  addface (i, c, {k,l},{k+1,l+1});
                  i++;
                }
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < order-1; k++)
              {
                addface (i, c, {order-1,k},{order,k+1});
                addface (i, (c+1)%3, {k,order-1},{k+1,order});
                i++;
              }
          for (int c = 0; c < 3; c++)
            addface (i,c, {order-1,order-1},{order,order});
      }
    Table<array<int,4>> faces_tab = creator.MoveTable();
    

    TableCreator<array<int,2>> creator_bd;
    
    for ( ; !creator_bd.Done(); creator_bd++)
      {
        for (auto i : Range(faces_tab))
          {
            Array<array<int,2>> tmparray;
            auto faces = faces_tab[i];
            for (auto k : Range(faces))
              for (int j = 0; j < 4; j++)
                {
                  array<int,2> tmpedge = { faces[k][j],faces[k][(j+1)%4] };
                  array<int,2> tmpedge2 = { tmpedge[1], tmpedge[0] }; 
                  if (!tmparray.Contains(tmpedge2))
                    tmparray.Append(tmpedge);
                  else
                    tmparray.DeleteElement(tmparray.Pos(tmpedge2)); 
                }
            
            for (auto ed : tmparray)
              creator_bd.Add(i,ed);
          }
      }
    
    

    return make_tuple(faces_tab,creator_bd.MoveTable());
  }

  template <>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes<ET_TRIG> (int order, bool dual)
  {
    throw Exception ("Nano volumes not implemented for ET_TRIG");
  }

  //TET
  template <>
  tuple<Array<Vec<3>>,Array<array<int,4>>,int> 
    GetNanoPoints<ET_TET> (int order, bool dual)
  {  
    Array<double> xi;
    IntegrationRule GaussRadauIRplus;
    ComputeSortedGaussRadauPoints (order, xi, false);
    for (auto i : Range(xi))
    {
      GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, 0));
    }
    GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));

    Array<Vec<3>> pnts;
    Array<array<int,4>> inds;
    
    auto addpoint = [&] (int corner, IVec<3> i)
    {
      pnts.Append ( MapHex2Tet ( Vec<3> (GaussRadauIRplus[i[0]](0),
                                         GaussRadauIRplus[i[1]](0),
                                         GaussRadauIRplus[i[2]](0)),
                                 corner) );
      inds.Append (C2M3D(corner, i, order));
    };
    
    for (int c = 0; c < 4; c++)
      for (int k = 0; k < order; k++)
        for (int l = 0; l < order; l++)
          for (int m = 0; m < order; m++)
            addpoint (c, { k, l, m } );

    for (int k = 0; k < order; k++)
      for (int l = 0; l < order; l++)
        {
          addpoint (0, { k, l, order } );
          addpoint (0, { order, k, l } );
          addpoint (0, { l, order, k } );

          addpoint (1, { order, k, l } );
          addpoint (1, { l, order, k } );

          addpoint (2, { order, k, l } );          
        }

    int nreal = pnts.Size();
    
    for (int k = 0; k < order; k++)
      {
        addpoint (0, { k, order, order } );
        addpoint (0, { order, k, order } );
        addpoint (0, { order, order, k } );
        addpoint (1, { order, order, k } );
      }

    addpoint(0, { order, order, order });
    return make_tuple(pnts, inds, nreal);
  }

  template <>
  tuple<Array<array<int,2>>,int> GetNanoEdges<ET_TET>
   (int order, bool dual)
  {
    throw Exception ("Nano edges not implemented for ET_TET");
  }
  
  template <>
  tuple<Array<array<int,4>>,int> GetNanoFaces<ET_TET> 
      (int order, bool dual)
  {
    throw Exception ("Nano faces not implemented for ET_TET");
  }

  template <>
  tuple<Array<array<int,8>>,int> GetNanoCells<ET_TET>
     (int order, bool dual)
  {
    throw Exception ("Nano cells not implemented for ET_TET");
  }


  template <>
  tuple<Table<array<int,2>>,Array<array<int,2>>> 
        GetNanoLines<ET_TET> (int order, bool dual)
  {
    auto [pnts, _inds, np] = GetNanoPoints<ET_TET> (order, dual);
    auto inds = std::move(_inds);
    TableCreator<array<int,2>> creator;
    auto addline = [&] (int i, int corner, IVec<3> p0, IVec<3> p1)
      {
        creator.Add (i,
                     array<int,2> {
                       int(inds.Pos (C2M3D ( corner, p0, order ))),
                         int(inds.Pos (C2M3D ( corner, p1, order ))),
                         });
      };
    
    for ( ; !creator.Done(); creator++)
      {
        int i = 0;
          for (int c = 0; c < 4; c++)
            for (int k = 0; k < order-1; k++)
              for (int l = 0; l < order; l++)
                for (int m = 0; m < order; m++)
                {
                  addline (i, c, {k,l,m},{k+1,l,m});
                  i++;
                  addline (i, c, {l,m,k},{l,m,k+1});
                  i++;
                  addline (i, c, {m,k,l},{m,k+1,l});
                  i++;
                }
              for (int k = 0; k < order; k++)
                for (int l = 0; l < order; l++)
              {
                //6 inner faces
                addline (i, 0, {order-1,k,l},{order,k,l});
                addline (i, 1, {k,l,order},{k,l,order-1});
                i++;
                addline (i, 0, {l,order-1,k},{l,order,k});
                addline (i, 2, {k,order,l},{k,order-1,l});
                i++;
                addline (i, 0, {k,l,order-1},{k,l,order});
                addline (i, 3, {order,k,l},{order-1,k,l});
                i++;
                addline (i, 1, {order-1,k,l},{order,k,l});
                addline (i, 2, {k,l,order},{k,l,order-1});
                i++;
                addline (i, 1, {l,order-1,k},{l,order,k});
                addline (i, 3, {k,order,l},{k,order-1,l});
                i++;
                addline (i, 2, {order-1,k,l},{order,k,l});
                addline (i, 3, {k,l,order},{k,l,order-1});
                i++;
              }
      }
    Table<array<int,2>> lines_table = creator.MoveTable();

    Array<array<int,2>> lines_ends;

    for (auto i : Range(lines_table))
    {

       Array<int> all_pts;
       for (auto k : Range(lines_table[i]))
           {
           array<int,2> tmpedge = lines_table[i][k];
           if(!all_pts.Contains(tmpedge[0]))
             all_pts.Append(tmpedge[0]);
           else 
             all_pts.DeleteElement(all_pts.Pos(tmpedge[0]));
           if(!all_pts.Contains(tmpedge[1]))
             all_pts.Append(tmpedge[1]);
           else 
             all_pts.DeleteElement(all_pts.Pos(tmpedge[1]));
           }
      lines_ends.Append(array<int,2> {all_pts[0],all_pts[1]});
    }
    return make_tuple(lines_table,lines_ends);
  }

  template <>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces<ET_TET> (int order, bool dual)
  {
    throw Exception ("Nano surfaces not implemented for ET_TET");
  }

  template <>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes<ET_TET> (int order, bool dual)
  {
    throw Exception ("Nano volumes not implemented for ET_TET");
  }

  
  

  std::map<ELEMENT_TYPE, IntegrationRule> GetWebGuiPoints(int order)
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;

    IntegrationRule irtrig;
    double eps = 1e-10;
    switch (order)
      {
      case 1:
        for (int v = 0; v < 3; v++)
          {
            Vec<2> pnts[] =
              {
                { 0, 0}, { 1-eps, 0 }, { 1.0-eps, 1.0-eps },
                { 0, 0}, { 1-eps, 1-eps }, { 0, 1.0-eps }          
              };
            
            for (auto xi : pnts)
              {
                Vec<2> x = MapQuad2Trig (xi, v);
                irtrig += IntegrationPoint(x(0), x(1), 0, 0 );
              }
          }
        break;
      case 2:
        for (int v = 0; v < 3; v++)
          {
            Vec<2> pnts[] =
              {
                { 0, 0}, { 0.5, 0 }, { 1-eps, 0 }, { 0.5, 0.5 }, { 1.0-eps, 0.5 }, { 1.0-eps, 1.0-eps },
                { 0, 0}, { 0.5, 0.5 }, { 1-eps, 1-eps }, { 0, 0.5 }, { 0.5, 1-eps }, { 0, 1.0-eps }          
              };
            
            for (auto xi : pnts)
              {
                Vec<2> x = MapQuad2Trig (xi, v);
                irtrig += IntegrationPoint(x(0), x(1), 0, 0 );
              }
          }
        break;
      case 3:
        for (int v = 0; v < 3; v++)
          {
            Vec<2> pnts1[] =
              {
                { 0, 0}, { 1-eps, 0 }, { 1.0-eps, 1.0-eps },
                { 0, 0}, { 1-eps, 1-eps }, { 0, 1.0-eps }          
              };
            
            Vec<2> pnts[20];
            int ii = 0;
            for (int k = 0; k < 2; k++)
              for (int i = 0; i <= 3; i++)
                for (int j = 0; j <= 3-i; j++, ii++)
                  {
                    double lami = double(i)/3;
                    double lamj = double(j)/3;
                    pnts[ii] = lamj * pnts1[3*k] + lami * pnts1[3*k+1] + (1-lami-lamj) * pnts1[3*k+2];
                  }
            
            for (auto xi : pnts)
              {
                Vec<2> x = MapQuad2Trig (xi, v);
                irtrig += IntegrationPoint(x(0), x(1), 0, 0 );
              }
          }
        break;
      default:
        throw Exception ("only order 1,2, and 3 are supported");
      }
      
    rules[ET_TRIG] = std::move(irtrig);


    IntegrationRule irtet;
    switch (order)
      {
      case 1:
        for (int v = 0; v < 4; v++)
          {
            Vec<3> pnts[] =
              {
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 0, 0 }, { 1, 0, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 0, 1 }, { 0, 0, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 0, 1 }, { 0, 1, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 1, 1 }, { 0, 1, 0},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 1, 0 }, { 1, 1, 0},
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 1, 0 }, { 1, 0, 0}
              };
            
            for (auto xi : pnts)
              {
                xi *= 1-1e-10;
                Vec<3> x = MapHex2Tet (xi, v);
                irtet += IntegrationPoint(x(0), x(1), x(2), 0 );
              }
          }
        break;
      case 2:
        for (int v = 0; v < 4; v++)
          {
            Vec<3> pnts1[] =
              {
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 0, 0 }, { 1, 0, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 0, 1 }, { 0, 0, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 0, 1 }, { 0, 1, 1},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 1, 1 }, { 0, 1, 0},
                { 0, 0, 0}, { 1, 1, 1 }, { 0, 1, 0 }, { 1, 1, 0},
                { 0, 0, 0}, { 1, 1, 1 }, { 1, 1, 0 }, { 1, 0, 0}
              };

            int between[6][2] = { { 0, 3}, { 1, 3 }, { 2, 3 }, { 0, 1 }, { 0, 2 }, { 1, 2 } };
            int ii = 0;
            for (int k = 0; k < 6; k++) // 6 tets of the cube
              {
                for (int i = 0; i < 4; i++)
                  {
                    Vec<3> xi = pnts1[4*k+i];
                    xi *= 1-1e-10;
                    Vec<3> x = MapHex2Tet (xi, v);
                    irtet += IntegrationPoint(x(0), x(1), x(2), 0 );
                  }
                for (int i = 0; i < 6; i++)
                  {
                    Vec<3> xi = 0.5 * (pnts1[4*k+between[i][0]] + pnts1[4*k+between[i][1]]);
                    xi *= 1-1e-10;
                    Vec<3> x = MapHex2Tet (xi, v);
                    irtet += IntegrationPoint(x(0), x(1), x(2), 0 );
                  }
              }
          }
        break;
      default:
        throw Exception ("only order 1 is supported");
      }
    rules[ET_TET] = std::move(irtet);

    return rules;
  }
}
