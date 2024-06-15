(time_stepping)=
# Time stepping

For the semi-discrete system of differential equations of the form

```{math}
:label: semi_disc_1o
\begin{aligned}
\mathbf M_p\partial_t{\mathbf p}&=\mathbf B\mathbf v,\\
\mathbf M_v\partial_t{\mathbf v}&=-\mathbf B^\top\mathbf p,\\
\mathbf p(0)&=\mathbf p_0,\quad {\mathbf v}(0)=\mathbf v_0,
\end{aligned}
```
where $\mathbf M_p,\mathbf M_v$ are the discrete mass matrices and $\mathbf B$ is the matrix of the discrete differential operator,
We use the Leap-Frog scheme given by
```{math}
:label: leap_frog
\mathbf v_{1/2}&=\mathbf v_0 -\frac{\tau}{2} \mathbf M_v^{-1}\mathbf B^\top\mathbf p_0,\\
\mathbf p_{j+1}&=\mathbf p_{j} +\tau \mathbf M_p^{-1}\mathbf B\mathbf v_{j+1/2},\\
\mathbf v_{j+1/2}&=\mathbf v_{j-1/2} -\tau \mathbf M_v^{-1}\mathbf B^\top\mathbf p_j.
```

It is well-known that this scheme is stable for timesteps $\tau>0$ fulfilling

```{math}
\tau^2 < \frac{4}{\sigma \left( \mathbf M_p^{-1}\mathbf B\mathbf M_v^{-1}\mathbf B^\top \right)}
```
where $\sigma$ is the spectral radius.
