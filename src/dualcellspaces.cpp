#include <python_comp.hpp>
#include "dcs_common.hpp"
#include "h1dualcells.hpp"
#include "h1dualcells3d.hpp"
#include "h1primalcells.hpp"

#include "hcurldualcells.hpp"
#include "hcurlprimalcells.hpp"

#include "hdivprimalcells.hpp"
#include "hdivdualcells.hpp"
#include "l2cells.hpp"

namespace std {
template <typename T, size_t S>
ostream & operator<< (ostream & ost, std::array<T,S> a)
{
  for (auto v : a)
    ost << v << ", ";
  return ost;
}
}

PYBIND11_MODULE(dualcellspaces,m) {
  // import ngsolve such that python base classes are defined
  auto ngs = py::module::import("ngsolve");

  using namespace ngcomp;




  auto h1primalpy = ExportFESpace<H1PrimalCells>(m, "H1PrimalCells"/*, true*/)
      .def("GetIntegrationRules", &H1PrimalCells::GetIntegrationRules, py::arg("fix_lo")=true)
    // .def("Rot", &H1PrimalCells::GetRotOperator, py::arg("dual")=true)
    ;
  m.attr("H1PrimalCells3D") = h1primalpy;

  ExportFESpace<H1DualCells>(m, "H1DualCells"/*, true*/)
      .def("GetIntegrationRules", &H1DualCells::GetIntegrationRules, py::arg("fix_lo")=true)
      .def("Gradient", &H1DualCells::GetGradientOperator2D, py::arg("dual")=true)
    ;

  auto h1dualpy = 
    ExportFESpace<H1DualCells3D>(m, "H1DualCells3D"/*, true*/)
      .def("GetIntegrationRules", &H1DualCells3D::GetIntegrationRules, py::arg("fix_lo")=true)
      .def("Gradient", &H1DualCells3D::GetGradientOperator3D,py::arg("dual")=true)
    ;
  m.attr("H1DualCells3D") = h1dualpy;
  


  auto hcurlprimalpy = 
    ExportFESpace<HCurlPrimalCells>(m, "HCurlPrimalCells"/*, true*/)
    .def("GetIntegrationRules", &HCurlPrimalCells::GetIntegrationRules, py::arg("fix_lo")=true)
    ;
  m.attr("HCurlPrimalCells3D") = hcurlprimalpy;  

  auto hcurldualpy =
    ExportFESpace<HCurlDualCells>(m, "HCurlDualCells"/*, true*/)
    .def("GetIntegrationRules", &HCurlDualCells::GetIntegrationRules, py::arg("fix_lo")=true)
    .def("GetPotentialSpace", &HCurlDualCells::GetPotentialSpace,
         py::arg("include_central")=false, py::arg("distributional")=false)
    .def("Curl", &HCurlDualCells::GetCurlOperator, py::arg("dual")=true, py::arg("altshapes")=true, py::arg("lumping")=true,
         py::arg("Kronecker")=false)
    // .def("Rot", &HCurlDualCells::GetRotOperator, py::arg("zerobd")=true)
    .def("Gradient", &HCurlDualCells::GetGradientOperator, py::arg("nanocells")=true)
    ;
  m.attr("HCurlDualCells3D") = hcurldualpy;

  

  
  ExportFESpace<HDivPrimalCells>(m, "HDivPrimalCells"/*, true*/)
    .def("GetIntegrationRules", &HDivPrimalCells::GetIntegrationRules, py::arg("fix_lo")=true)
    //.def("ConvertGR2GOperator", &HDivPrimalCells::ConvertGR2GOperator)
    ; 
  ExportFESpace<HDivDualCells>(m, "HDivDualCells")
    .def("GetIntegrationRules", &HDivDualCells::GetIntegrationRules, py::arg("fix_lo")=true)
    ; 

  py::class_<PotentialFESpace, FESpace, shared_ptr<PotentialFESpace>>(m, "PotentialFESpace")
    .def("GetInnerDofs", &PotentialFESpace::GetInnerDofs)
    ;  

  py::class_<HCurlDualCellsPotential3D, PotentialFESpace, shared_ptr<HCurlDualCellsPotential3D>> (m, "HCurlDualCellsPotential3D")
    ;

  ExportFESpace<L2Cells> (m, "L2Cells")
    ;

  
  
  m.def("GetIntegrationRules", [](int order, bool innerfacets)
  {
    if (!innerfacets)
      return GetIntegrationRules(order);
    else
      return GetIntegrationRulesInnerFacets(order);
  }, py::arg("order"), py::arg("innerfacets")=false);
  
  m.def("GetWebGuiPoints", &GetWebGuiPoints);

  m.def("MicroJacobiDet", [](int dim) -> shared_ptr<CoefficientFunction>
        {
          switch (dim)
            {
            case 2:
              return make_shared<MicroJacobiDet<2>> ();
            case 3:
              return make_shared<MicroJacobiDet<3>> ();              
            }
          return nullptr;
        });


  

  
  ExportArray<std::array<int,2>>(m);
  ExportArray<std::array<int,4>>(m);
  ExportArray<std::array<int,8>>(m);
  ExportTable<std::array<int,2>>(m);
  ExportTable<std::array<int,4>>(m);
  ExportTable<std::array<int,8>>(m);
  //Nanocells legacy
  m.def("GetMicroCellPoints", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                //issue with Array<array<int,len>>
                //return GetMicroCellPoints<ET_TRIG> (order, dual);
                auto [pnts, inds, num] = GetMicroCellPoints<ET_TRIG> (order, dual);
                std::vector<array<int,3>> i2;
                for (auto ind : inds)
                  i2.push_back(ind);
                std::vector<Vec<2>> p2;
                for (auto p : pnts)
                  p2.push_back(p);
                // return tuple(p2, i2, num);
                return tuple(py::cast(p2), py::cast(i2), py::cast(num));
              }
            case ET_TET:
              {
                //issue with Array<array<int,len>>
                //return GetMicroCellPoints<ET_TRIG> (order, dual);
                auto [pnts, inds, num] = GetMicroCellPoints<ET_TET> (order, dual);
                std::vector<array<int,4>> i2;
                for (auto ind : inds)
                  i2.push_back(ind);
                std::vector<Vec<3>> p2;
                for (auto p : pnts)
                  p2.push_back(p);
                return tuple(py::cast(p2), py::cast(i2), py::cast(num));                
              }
            default:
              throw Exception("MicroCellPoints not available for "+ToString(type));
            }
        });

  
  
  m.def("GetMicroCellEdges", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                auto [edges, num] = GetMicroCellEdges<ET_TRIG> (order, dual);
                return make_tuple(std::move(edges), num);
              }
            case ET_TET:
              {
                auto [edges, num] = GetMicroCellEdges<ET_TET> (order, dual);
                return make_tuple(std::move(edges), num);
              }
            default:
              throw Exception("MicroCellEdges not available for "+ToString(type));
            }
        });

  m.def("GetMicroCellFaces", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                return GetMicroCellFaces<ET_TRIG> (order, dual);
              }
            case ET_TET:
              {
                return GetMicroCellFaces<ET_TET> (order, dual);
              }
            default:
              throw Exception("MicroCellFaces not available for "+ToString(type));
            }
        });
  
  m.def("GetMicroCellFaceBoundaries", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                return GetMicroCellFaceBoundaries<ET_TRIG> (order, dual);
              }
            case ET_TET:
              {
                return GetMicroCellFaceBoundaries<ET_TET> (order, dual);
              }
            default:
              throw Exception("MicroCellFaceBoundaries not available for "+ToString(type));
            }
        });
  ngs.attr("_add_flags_doc")(m);



  //Nanocell new
  m.def("GetNanoPoints", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                //issue with Array<array<int,len>>
                //return GetMicroCellPoints<ET_TRIG> (order, dual);
                auto [pnts, inds, num] = GetNanoPoints<ET_TRIG> (order, dual);
                std::vector<array<int,3>> i2;
                for (auto ind : inds)
                  i2.push_back(ind);
                std::vector<Vec<2>> p2;
                for (auto p : pnts)
                  p2.push_back(p);
                // return tuple(p2, i2, num);
                return tuple(py::cast(p2), py::cast(i2), py::cast(num));
              }
            case ET_TET:
              {
                //issue with Array<array<int,len>>
                //return GetMicroCellPoints<ET_TRIG> (order, dual);
                auto [pnts, inds, num] = GetNanoPoints<ET_TET> (order, dual);
                std::vector<array<int,4>> i2;
                for (auto ind : inds)
                  i2.push_back(ind);
                std::vector<Vec<3>> p2;
                for (auto p : pnts)
                  p2.push_back(p);
                return tuple(py::cast(p2), py::cast(i2), py::cast(num));                
              }
            default:
              throw Exception("NanoPoints not available for "+ToString(type));
            }
        });

  m.def("GetNanoEdges", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                auto [edges, num] = GetNanoEdges<ET_TRIG> (order, dual);
                return make_tuple(std::move(edges), num);
              }
            case ET_TET:
              {
                auto [edges, num] = GetNanoEdges<ET_TET> (order, dual);
                return make_tuple(std::move(edges), num);
              }
            default:
              throw Exception("NanoEdges not available for "+ToString(type));
            }
        });
  m.def("GetNanoFaces", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            case ET_TRIG:
              {
                auto [faces, num] = GetNanoFaces<ET_TRIG> (order, dual);
                return make_tuple(std::move(faces), num);
              }
            /*case ET_TET:
              {
                auto [edges, num] = GetNanoFaces<ET_TET> (order, dual);
                return make_tuple(std::move(edges), num);
              }*/
            default:
              throw Exception("NanoFaces not available for "+ToString(type));
            }
        });
  m.def("GetNanoCells", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            /*
            case ET_TRIG:
              {
                auto [edges, num] = GetNanoCells<ET_TRIG> (order, dual);
                return make_tuple(std::move(edges), num);
              }
            case ET_TET:
              {
                auto [edges, num] = GetNanoCells<ET_TET> (order, dual);
                return make_tuple(std::move(edges), num);
              }
              */
            default:
              throw Exception("NanoCells not available for "+ToString(type));
            }
        });

  m.def("GetNanoLines", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            
            case ET_TRIG:
              {
                auto [lines, line_ends] = GetNanoLines<ET_TRIG> (order, dual);
                std::vector<array<int,2>> i2;
                for (auto ind : line_ends)
                  i2.push_back(ind);
                return tuple(std::move(lines), py::cast(i2));
                //return tuple(py::cast(i2));
              }
            case ET_TET:
              {
                auto [lines, line_ends] = GetNanoLines<ET_TET> (order, dual);
                std::vector<array<int,2>> i2;
                for (auto ind : line_ends)
                  i2.push_back(ind);
                return tuple(std::move(lines), py::cast(i2));
                //return tuple(py::cast(i2));
              }
            default:
              throw Exception("NanoLines not available for "+ToString(type));
            }
        });
  m.def("GetNanoSurfaces", [](ELEMENT_TYPE type, int order, bool dual)
        {
          switch (type)
            {
            
            case ET_TRIG:
              {
                auto [lines, line_ends] = GetNanoSurfaces<ET_TRIG> (order, dual);
                return make_tuple(std::move(lines), std::move(line_ends));
                //return tuple(py::cast(i2));
              }
              /*
            case ET_TET:
              {
                auto [lines, line_ends] = GetNanoSurfaces<ET_TET> (order, dual);
                return make_tuple(std::move(lines), std::move(line_ends));
                //return tuple(py::cast(i2));
              }
              */
            default:
              throw Exception("NanoLines not available for "+ToString(type));
            }
        });
}

