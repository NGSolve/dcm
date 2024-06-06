import glob
import os
import sys
import site
import ngsolve, netgen

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib
import pkg_resources

version = os.environ["CI_COMMIT_TAG"]
name = "dualcellspaces"

py_install_dir = get_python_lib(1, 0, "").replace("\\", "/")

ngsolve_version = ngsolve.config.NGSOLVE_VERSION_PYTHON
root_dir = os.path.abspath(
    os.path.join(netgen.__file__, "../" * (len(py_install_dir.split("/")) + 2))
)

if netgen.config.NG_INSTALL_DIR_CMAKE.startswith("netgen"):
    netgen_dir = os.path.abspath(
        os.path.join(os.path.dirname(netgen.__file__), "cmake")
    )
    ngsolve_dir = os.path.abspath(
        os.path.join(os.path.dirname(ngsolve.__file__), "cmake")
    )
else:
    netgen_dir = root_dir
    ngsolve_dir = root_dir

install_requires = [
    f"ngsolve == {ngsolve_version}",
]

_cmake_args = [
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_CXX_COMPILER=clang++",
    f"-DNGSolve_DIR={ngsolve_dir}",
    f"-DNetgen_DIR={netgen_dir}",
]

if "PYDIR" in os.environ:
    _cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}']

setup(
    name=name,
    version=version,
    description="Dualcellwaves",
    author="",
    license="custom",
    install_requires=install_requires,
    py_modules=["duallcellspaces"],
    # tests_require=['pytest','scipy','numpy'],
    cmake_args=_cmake_args
)
