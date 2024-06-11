(installation)=
# Installation

The `dualcellspaces` package is an add-on to the high-order finite element library [NGSolve](https://ngsolve.org). Thus the main premise is that a sufficiently recent version of [NGSolve](https://ngsolve.org) is installed beforehand.


## Get NGSolve

The most simple way to install the `NGSolve` is to use `pip` (the [python package installer](https://pypi.org/project/pip/)):
```
python -m pip install numpy scipy matplotlib jupyter ipyparallel scikit-build
python -m pip install --upgrade --pre ngsolve webgui_jupyter_widgets
```
For troubleshooting we refer to the various `NGSolve` installation tutorials [NGS24](https://docu.ngsolve.org/ngs24/intro.html) or the [NGSolve Documentation](https://docu.ngsolve.org/latest/).

If you do not want to update or modify your `NGSolve` installation you might want to consider installing `NGSolve` and the `dualcellspaces` add on in a [virtual environment](https://docu.ngsolve.org/ngs24/intro.html#installing-parallel-ngsolve).

## Install the Dual Cell Method add-on

As the Dual Cell Method is not available in the core code of `NGSolve` the `dualcellspaces` add-on has to be installed. 

### Using pip

The most simple way to install the add on is again to use `pip` via

```` {tab-set}
```{tab-item} NGSolve pip installation

  if you have installed NGSolve via pip or binaries

    pip install git+https://github.com/NGSolve/dcm.git
```
```{tab-item} NGSolve  built from sources

  if your NGSolve was built from sources

    python -m pip install scikit-build-core pybind11_stubgen 
    python -m pip install --no-build-isolation git+https://github.com/NGSolve/dcm.git
```
````

### Working with the code

If you also want to have the sources of the `dualcellspaces` module available you may clone them using

```
git clone https://github.com/NGSolve/dcm.git
```

Installation can be done either again using pip or building using `CMake`

```` {tab-set}
``` {tab-item} pip
    cd dcm
    python -m pip install --no-build-isolation .
```
``` {tab-item} CMake
    cd dcm
    mkdir build
    cd build
    cmake ..
    make -j4 install
```
````


### Test the `dualcellspaces` installation

Test whether the installation worked using
```
python -m dualcellspaces.demos.dc_intrules
```

