(installation)=
# Installation


## Get NGSolve

The main premise is that a sufficiently recent version of [NGSolve](https://ngsolve.org) is already installed. If this is not the case the most simple way to install the packe is done using `pip` (the [python package installer](https://pypi.org/project/pip/)):
```
python -m pip install numpy scipy matplotlib jupyter ipyparallel scikit-build
python -m pip install --upgrade --pre ngsolve webgui_jupyter_widgets
```

## Install the Dual Cell Method add-on

As the Dual Cell Method is not available in the core code of `NGSolve` the `dcm` add-on has to be installed. 

### Using pip

The most simple way to install the add on is again to use `pip` via
```
pip install git+https://github.com/NGSolve/dcm.git
```

### Compiling the code

Compiling the code may be done using

```
git clone https://github.com/NGSolve/dcm.git
cd dcm
mkdir build
cd build
cmake ..
make -j4 install
```
