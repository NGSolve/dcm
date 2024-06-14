(maxwell)=
# The dual cell method for the time-domain Maxwell system

In this section we describe in detail how the dual cell method is applied to the time-domain Maxwell equations.

## Problem setting

The time-domain Maxwell system (we neglect possible sources here) is the problem to find fields $\mathbf D,\mathbf B,\mathbf H,\mathbf E$ such that
```{card}
\begin{align*}
\partial_t \mathbf{D}(t,\mathbf{x})-\mathrm{curl}\mathbf{H}(t,\mathbf{x})&=0,\\
\partial_t \mathbf{B}(t,\mathbf{x})+\mathrm{curl}\mathbf{E}(t,\mathbf{x})&=0.\\
&+\text{b.c., i.c.,}
\end{align*}
```
for $t\in(0,T),\mathbf x\in\Omega$ and some $T>0$ and a suitable domain $\Omega\subset\mathbb R^3$.
To close the system one also needs the constitutive relations
```{card}
\begin{align*}
  \mathbf{D} &= {\color{emph1}\varepsilon}\mathbf{E},&
  \mathbf{B} &= {\color{emph2}\mu}\mathbf{H},
\end{align*}
```
where $\varepsilon,\mu$ are the **permittivity** and **permeability** of the medium in question.

In weak (EH-)formulation, assuming homogeneous boundary conditions $\mathbf E \times \mathbf n = 0$, the above problem may be rewritten as the problem to find $\mathbf E,\mathbf H:[0,T]\to H(\mathrm{curl})(\Omega)$  such that
```{card}
\begin{align*}
  \int_\Omega\varepsilon \partial_t \mathbf{E}(t,\mathbf{x})\mathbf E'(\mathbf x)d\mathbf x-\int_\Omega\mathrm{curl}\mathbf{H}(t,\mathbf{x})\mathbf E'(\mathbf x)d\mathbf x&=0,\\
\int_\Omega\mu\partial_t \mathbf{H}(t,\mathbf{x})\mathbf H'(\mathbf x)d\mathbf x+\int_\Omega\mathbf{E}(t,\mathbf{x})\mathrm{curl}\mathbf H'(\mathbf x)d\mathbf x&=0.\\
&+\text{i.c.,}
\end{align*}
```
for all $\mathbf E',\mathbf H'\in H(\mathrm{curl})(\Omega)$.
