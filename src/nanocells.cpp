#include <comp.hpp>
#include "dcs_common.hpp"
#include "intrules.hpp"



namespace ngcomp
{


  //TRIG
  // i ... 2D coordinates inm corner
  // returns 3D coordinates on manifold
  //
  //inline constexpr array<int,3> C2M (int corner, INT<2> i, int order)
  //
  /*
  inline array<int,3> C2M (int corner, INT<2> i, int order)
  {
    if (corner == 0) return { order, i[0], i[1] };
    if (corner == 1) return { i[1], order, i[0] };
    return { i[0], i[1], order };
  }
*/

/*
  template <>
  tuple<Array<Vec<2>>,Array<array<int,3>>,int> GetNanoPoints<ET_TRIG> (int order, bool dual)
  {
    // order x order squares per corner
    Array<double> xi, wi;
    IntegrationRule GaussRadauIRplus;
    ComputeGaussRadauRule (order, xi, wi);
    for (auto i : Range(xi))
      GaussRadauIRplus.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
    GaussRadauIRplus.Append (IntegrationPoint (1, 0, 0, 0));

    Array<Vec<2>> pnts;
    Array<array<int,3>> inds;
    
    auto addpoint = [&] (int corner, INT<2> i)
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

    int nreal = pnts.Size();

    addpoint(0, { order, order });

    return make_tuple(pnts, inds, nreal);
  }


  */
  template <>
  tuple<Array<array<int,2>>,int> GetNanoEdges<ET_TRIG> (int order, bool dual)
  {
    Array<array<int,2>> edges;
    
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = move(_inds);
    
    auto addedge = [&] (int corner, INT<2> p0, INT<2> p1)
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

  /*
  template <>
  tuple<Array<array<int,4>>,int> GetNanoFaces<ET_TRIG> 
      (int order, bool dual)
  {
    ;
  }

  template <>
  tuple<Array<array<int,8>>,int> GetNanoCells<ET_TRIG>
     (int order, bool dual)
  {
    ;
  }


  template <>
  tuple<Table<array<int,2>>,array<int,2>> 
        GetNanoLines<ET_TRIG> (int order, bool dual)
  {
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = move(_inds);
    TableCreator<array<int,2>> creator;

    auto addline = [&] (int i, int corner, INT<2> p0, INT<2> p1)
      {
        creator.Add (i,
                     array<int,2> {
                       int(inds.Pos (C2M ( corner, p0, order ))),
                         int(inds.Pos (C2M ( corner, p1, order ))),
                         });
      };
    auto [pnts, _inds, np] = GetNanoPoints<ET_TRIG> (order, dual);
    auto inds = move(_inds);
    TableCreator<array<int,2>> creator;

    auto addline = [&] (int i, int corner, INT<2> p0, INT<2> p1)
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
                addline (i, (c+1)%3, {k,order-1},{k,order});
                i++;
              }
      }
    Table<array<int,2>> lines_table = creator.MoveTable();
    
    return lines_table;
  }

  template <>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces<ET_TRIG> (int order, bool dual);
  {
    // cout << "Called GetMicroCellFaces" << endl;
    auto [pnts, _inds, np] = GetPoints<ET_TRIG> (order, dual);
    auto inds = move(_inds);
    TableCreator<array<int,4>> creator;

    auto addface = [&] (int i, int corner, INT<2> p0, INT<2> p2)
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

  template <>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes<ET_TRIG> (int order, bool dual)
  {
    ;
  }

  //TET
  template <>
  tuple<Array<Vec<3>>,Array<array<int,4>>,int> 
    GetNanoPoints<ET_TET> (int order, bool dual)
  {  ;  }

  template <>
  tuple<Array<array<int,2>>,int> GetNanoEdges<ET_TET>
   (int order, bool dual)
  {
    ;
  }
  
  template <>
  tuple<Array<array<int,4>>,int> GetNanoFaces<ET_TET> 
      (int order, bool dual)
  {
    ;
  }

  template <>
  tuple<Array<array<int,8>>,int> GetNanoCells<ET_TET>
     (int order, bool dual)
  {
    ;
  }


  template <>
  tuple<Table<array<int,2>>,array<int,2>> 
        GetNanoLines<ET_TET> (int order, bool dual)
  {
    ;
  }

  template <>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces<ET_TET> (int order, bool dual);
  {
    ;
  }

  template <>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes<ET_TET> (int order, bool dual)
  {
    ;
  }


  */
}
