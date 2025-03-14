cmake_minimum_required(VERSION 3.16)

# change this to your desired addon name
# make sure, that the name matches exactly the argument of PYBIND11_MODULE()
set(addon_name dualcellspaces)

project(${addon_name})
include(ngsolve_addon.cmake)

# change the source file arguments to your source files
add_ngsolve_addon(${addon_name} 
   src/dualcellspaces.cpp
   src/h1dualcells.cpp 
   #   src/h1dualcells3d.cpp 
   src/h1primalcells.cpp 
   src/hcurldualcells.cpp    
   src/hcurlprimalcells.cpp    
   src/hdivprimalcells.cpp
   src/hdivdualcells.cpp
   src/l2cells.cpp    
   src/intrules.cpp
   src/supersparse.cpp
  )

# for opimization
# set(CMAKE_CXX_FLAGS "-O3")
# install the compiled python module and __init__.py ( don't change this )
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/${addon_name})
install(TARGETS ${addon_name} DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})

# install additional python files/demos/examples
install(FILES demos/dc_intrules.py DESTINATION ${install_dir}/demos)
install(FILES demos/test_spaces.py DESTINATION ${install_dir}/demos)

# generate stub files for autocomplete in IDEs
# this must be done at the very end (such that the stubgen generation happens after the python modules are installed)
ngsolve_generate_stub_files(${addon_name})
