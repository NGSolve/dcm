#ifndef FILE_NANOCELLS
#define FILE_NANOCELLS
#include <comp.hpp>

namespace ngcomp
{
  template <ELEMENT_TYPE ET>
  tuple<Array<Vec<Dim(ET)>>,Array<array<int,Dim(ET)+1>>,int> GetNanoPoints (int order, bool dual = true);


  template <int ELEMENT_TYPE>
  tuple<Array<array<int,2>>,int> GetNanoEdges (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,4>>,int> GetNanoFaces (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,8>>,int> GetNanoCells (int order, bool dual = true);




  template <int ELEMENT_TYPE>
  tuple<Table<array<int,2>>,array<int,2>> 
        GetNanoLines (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes (int order, bool dual = true);

}
#endif
