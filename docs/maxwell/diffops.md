---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


(diffops)=
# The discrete differential operators


Missing for the whole spacial discretization are the discretization matrices for the differential operators. 
Assembling the discrete $\mathrm{curl}$ matrices requires integration over element boundaries. This is done the following way


```{code-cell} ipython
from ngsolve import *
import dualcellspaces as dcs
from ngsolve.webgui import Draw
from time import time



mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))
fesE = dcs.HCurlDualCells(mesh,order=2)
fesH = dcs.HCurlPrimalCells(mesh,order=2)

H = fesH.TrialFunction()
dE = fesE.TestFunction()

irs = dcs.GetIntegrationRules(5)
dx_vol = dx(intrules=irs)
dx_edge = dx(element_boundary=True,intrules = irs)

normal = specialcf.normal(3)
Curl = BilinearForm(curl(H)*dE*dx_vol-H*Cross(dE,normal)*dx_edge)
```

Due to the covariant transformation, the geometric contributions (i.e., the Jacobian matrices of the transformation from reference to physical cell) cancel out. Thus the contribution of each element to the global matrix is the same (up to possible permutations of the basis functions). Thus we merely need to assemble and store one element matrix for each equivalent class of possible element matrices which results in a reduction in computational and memory costs. This is realized by the flag `geom_free = True`:

```{code-cell} ipython

Curl_gf = BilinearForm(curl(H)*dE*dx_vol-H*Cross(dE,normal)*dx_edge, geom_free=True)

gfE = GridFunction(fesE)
gfH = GridFunction(fesH)

gfH.vec.SetRandom()

with TaskManager():
  now = time()
  Curl.Assemble()
  curlt = time()-now

  now = time()
  Curl_gf.Assemble()
  gft = time()-now

  print('#### assembling ####')
  print('full: {}s'.format(curlt))
  print('geometry_free: {}s'.format(gft))

  n = 100
  now = time()
  for i in range(n):
    gfE.vec.data = Curl.mat * gfH.vec    
  curlt = time()-now

  now = time()
  for i in range(n):
    gfE.vec.data = Curl_gf.mat * gfH.vec    
  gft = time()-now

  print('#### application ####')
  print('full: {}s'.format(curlt))
  print('geometry_free: {}s'.format(gft))
```


## Étude: computing the gradient of a Gaussian peak

Before we launch into the full time-domain wave problem we take some time to verify our implementation by projecting a Gaussian peak 
```{math}
f(\mathbf x) = \frac{1}{2}\exp(-100 \|\mathbf x - (\tfrac{1}{2},\tfrac{1}{2})^\top\|^2),
```
 in $\mathbb R^2$ into the discrete space $\tilde X^{\mathrm{grad}}_P(\tilde{\mathcal T})$ and computing the discrete gradient in $X^{\mathrm{div}}_P({\mathcal T})$.

We define the mesh and spaces

```{code-cell} ipython

mesh = Mesh(unit_square.GenerateMesh(maxh = 0.03))

order = 4
h1 = dcs.H1DualCells(mesh, order = order)
hdiv = dcs.HDivPrimalCells(mesh, order = order)

print("DoFs H1Primal: {}".format(h1.ndof))
print("DoFs HDivDual: {}".format(hdiv.ndof))
```

As in {numref}`mass_lumping` we obtain the mass matrix. We assemble the right hand side for the projection using the lumped integration rule and solve the projection problem to find $p\in \tilde X^{\mathrm{grad}}_P(\tilde {\mathcal T})$

```{math}
(p,q)_h = (f,q)_h
```
for all $q\in  X^{\mathrm{grad}}_P(\mathcal T)$q.

```{code-cell} ipython
mass_h1_inv = h1.Mass(1).Inverse()

dx_h1 = dx(intrules = h1.GetIntegrationRules())
peak = CF( 0.5 * exp(-100*( (x-0.5)**2 + (y-0.5)**2 ))  )

p,q = h1.TnT()
rhs = LinearForm(peak*q*dx_h1).Assemble().vec

gfp = GridFunction(h1)
gfp.vec.data = mass_h1_inv * rhs

Draw(gfp, order = 2, points = dcs.GetWebGuiPoints(2), deformation = True, euler_angles = [-40,-4,-150]);
```

Similar to above we assemble the discrete gradient including the distributional (boundary) terms

```{math}
b(p,v) = \sum_{T\in\mathcal T} -\int_T p \mathrm{div} v dx +\int_{\partial T} p v\cdot n ds.
```

```{code-cell} ipython
n = specialcf.normal(2)

dSw = dx(element_boundary = True, intrules = dcs.GetIntegrationRules(2*order - 1))
dxw = dx(intrules = dcs.GetIntegrationRules(2*order -1))

v = hdiv.TestFunction()
grad = BilinearForm(-p*div(v)*dxw + p*(v*n)*dSw, geom_free = True).Assemble().mat

```
Lastly we solve the weak problem to find $u\in X^{\mathrm{div}}_P(\mathcal T)$ such that
```{math}
(u,v)_h = b(p,v)
```
for all $v$.


```{code-cell} ipython
gfu = GridFunction(hdiv)

mass_hdiv_inv = hdiv.Mass().Inverse()

gfu.vec.data = mass_hdiv_inv @ grad * gfp.vec


Draw(gfu, order = 2, points = dcs.GetWebGuiPoints(2), vectors = True, euler_angles = [-40,-4,-150]);
```
