from ngsolve import *
import dualcellspaces as dcs


print("### 2d spaces ###")
mesh = Mesh(unit_square.GenerateMesh())

print("number of elements = {}".format(mesh.ne))
h1primal = dcs.H1PrimalCells(mesh,order=0)
print("number of H1PrimalCells dofs = {}".format(h1primal.ndof))
h1dual = dcs.H1DualCells(mesh,order=0)
print("number of H1DualCells dofs = {}".format(h1dual.ndof))
hDivprimal = dcs.HDivPrimalCells(mesh,order=0)
print("number of HDivPrimalCells dofs = {}".format(hDivprimal.ndof))
#hDivdual = dcs.HDivDualCells(mesh,order=0)
#print("number of HDivDualCells dofs = {}".format(hDivdual.ndof))
hCurlprimal = dcs.HCurlPrimalCells(mesh,order=0)
print("number of HCurlPrimalCells dofs = {}".format(hCurlprimal.ndof))
hCurldual = dcs.HCurlDualCells(mesh,order=0)
print("number of HCurlDualCells dofs = {}".format(hCurldual.ndof))
