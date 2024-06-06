from ngsolve import *
import dualcellspaces as dcs


mesh = Mesh(unit_square.GenerateMesh(maxh=1))

h1 = dcs.H1PrimalCells(mesh)

print(h1.GetIntegrationRules()[TRIG])
