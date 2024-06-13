(first_steps)= 
# A quick start into the dual cell method

````{warning}
This is not a full introduction to the method. For a more detailed mathematical exposition of the method we refer to {cite}`WKCS2023` and {numref}`
````

## Problem setting

The dual cell method is designed for the spacial discretization of time-domain wave-type equations in first order form, i.e.,
for $\Omega\in \mathbb R^{d}, d=1,2,3$ and $T>0$

````{card}
```{math}
\begin{aligned}
\partial_t p(t,x) - B^\top v(t,x) &= 0,&t&\in (0,T),x\in\Omega\\
\partial_t v(t,x) + B p(t,x) &= 0,&t&\in (0,T),x\in\Omega\\
\text{+ i.c., + b.c.}
\end{aligned}
```
```` 
where $p:[0,T]\to V_p, v:[0,T]\to V_v$ with suitable spaces $V_p,V_v$ and a (differential in space) operator $B:V_v\to V_p$.

Examples for such equations include the time-domain Maxwell system or the acoustic wave equation. where the Operator $B$ is the $\mathrm{curl}$ or $\nabla$ Operator.

We aim for a [method of lines](https://en.wikipedia.org/wiki/Method_of_lines) approach, i.e, we discretize the problem in space to obtain a system of ODEs in time. This resulting system is treated using standard time-stepping schemes for stiff ODEs (e.g., [Leap Frog](https://en.wikipedia.org/wiki/Leapfrog_integration) or variants).

The ansatz spaces $V_v^h$ and $V_p^h$ for the spacial discretization are chosen locally conforming on a primal and dual mesh respectively.


## Spacial discretization

As an example we choose the acoustic wave equation in a two-dimensional domain $\Omega$ in weak first order form
````{card}
```{math}
\begin{aligned}
\partial_t \int_\Omega p(t,x)q(x)dx -\int_\Omega v(t,x)\cdot \nabla q(x)dx &= 0,&t&\in (0,T)\\
\partial_t \int_\Omega v(t,x)\cdot w(x)dx + \int_\Omega  \nabla p(t,x)\cdot w(x)dx &= 0,&t&\in (0,T)\\
p(0,x) &= p_0,&x&\in\Omega,\\
v(0,x) &= v_0,&x&\in\Omega,
\end{aligned}
```
```` 

We choose 
