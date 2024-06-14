# The Dual Cell Method in NGSolve

by *M. Wess*, *J. Schöberl*

TU Wien, Institute of Analysis and Scientific Computing, 

*based on joint work with B. Kapidani and L. Codecasa*


---


This book is designed to provide an introduction and examples to the implementation of the Dual Cell Method in the high-order finite element library [NGSolve](https://ngsolve.org).


The *Dual Cell Method* (DCM) is a Galerkin Method for the simulation of time-domain waves (e.g., electromagnetic or acoustiv waves) in mixed formulation. It is a Disconitinuous Galerkin variant where the two wave-fields are approximated by conforming functions on meshes dual to each other. Thus the respective ansatz functions feature discontinuities on different element boundaries.

For a full mathematical introduction to the method we refer to {cite}`WKCS2023,KCS2021`.

---

## Table of Contents
```{tableofcontents}
```


---

## References
```{bibliography}
:filter: docname in docnames
```
