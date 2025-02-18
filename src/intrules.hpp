#ifndef FILE_INTRULES
#define FILE_INTRULES
#include <comp.hpp>
namespace ngcomp
{

  void ComputeSortedGaussRadauPoints(int order,Array<double> &xi,bool ascending = true);


  template <typename T>
  void LagrangePolynomials (T x, const IntegrationRule & ir, FlatArray<T> shapes);

  IntegrationRule PrimalVolIR(const IntegrationRule & IR, bool mirroredIR = false);
  
  IntegrationRule PrimalCellIR(const IntegrationRule & IR, bool mirroredIR = false);

  IntegrationRule PrimalSegmIR(const IntegrationRule & irseg, bool mirroredIR = false);


  std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(int intorder);

  std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRulesInnerFacets(int intorder);  

  std::map<ELEMENT_TYPE, IntegrationRule> GetWebGuiPoints(int intorder);



  //legacy names
  template <ELEMENT_TYPE ET>
  tuple<Array<Vec<Dim(ET)>>,Array<array<int,Dim(ET)+1>>,int> GetMicroCellPoints (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,2>>,int> GetMicroCellEdges (int order, bool dual = true);
  
  template <int ELEMENT_TYPE>
  Table<array<int,4>> GetMicroCellFaces (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  Table<array<int,2>> GetMicroCellFaceBoundaries (int order, bool dual = true);
  


  //new convention
  template <ELEMENT_TYPE ET>
  tuple<Array<Vec<Dim(ET)>>,Array<array<int,Dim(ET)+1>>,int> GetNanoPoints (int order, bool dual = true);

  template<>
  tuple<Array<Vec<2>>,Array<array<int,3>>,int> GetMicroCellPoints<ET_TRIG> (int order, bool dual);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,2>>,int> GetNanoEdges (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,4>>,int> GetNanoFaces (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Array<array<int,8>>,int> GetNanoCells (int order, bool dual = true);




  template <int ELEMENT_TYPE>
  tuple<Table<array<int,2>>,Array<array<int,2>>> 
        GetNanoLines (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Table<array<int,4>>,Table<array<int,2>>>
        GetNanoSurfaces (int order, bool dual = true);

  template <int ELEMENT_TYPE>
  tuple<Table<array<int,8>>,Table<array<int,4>>>
      GetNanoVolumes (int order, bool dual = true);
}
#endif
