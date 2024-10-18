# [Algorithm](@id algorithm-theory)

The main motivation of this project is to develop an efficient algorithm for vector-valued
_p_-Laplace problems that occur in shape optimization using a _p_-harmonic approach
\[[2](@ref references-theory)\] such as
```math
\underset{v \in W^{1,p}_g(\Omega,\R^d)}
    {\mathrm{arg\;min}}\; \frac{1}{p} \int_{\Omega} \Vert \nabla v \Vert_2^p \;\mathrm{d} x
    - \int_{\Omega} f v \;\mathrm{d} x
    - \int_{\Gamma} h v  \;\mathrm{d}\Gamma
```
Nevertheless, various kinds of problems can be solved using this package.
Assuming sufficient regularity, the given minimization problem is associated
with the general Euler-Lagrange-equation
```math
\left.\begin{array}{rl}
    -\Delta_p v = f &\text{in } \Omega,\\
    \Vert\nabla v \Vert_2^{p-2} \partial_{\eta} v = h & \text{on } \Gamma,\\
    v = g &\text{on } \partial\Omega\setminus\Gamma
    \end{array}\right\}
```
for the _p_-Laplacian ``\Delta_p v = \nabla \cdot (\Vert \nabla v \Vert_2^{p-2} \nabla v)``.

The main intention of this package in julia is not high-performance computing, 
but rather creating understandable code and reproducible results.
This enables improved opportunities for analysis and quick prototyping.

The implementation is based on an approach proposed by Loisel
\[[1](@ref references-theory)\] using a reformulation of the problem as a convex
optimization problem and then solve it using an interior point method.
While originally proposed for pure Dirichlet conditions and just scalar-valued problems,
it was extended to general vector-valued problems as well as mixed Dirichlet and Neumann
boundary conditions \[[4](@ref references-theory)\].

In the following, we give a short introduction into the implementation of the algorithm.
We start with some basics and notation for the used finite element library and extension.
Then we show the reformulation as a standard convex problem in finite dimension.
Finally, we outline the numerical technique for actually solving the problem.

## [FEM](@id fem-theory)

As the finite element base layer, the package
[`MinFEM.jl`](https://minfem.github.io/MinFEM.jl/stable/) is used.
In particular, linear triangular Lagrange elements with non-negative weights.

Let a quasi-uniform triangulation with ``n`` nodes and  ``m`` elements of
a computational domain ``\Omega \subset \mathbb{R}^d`` be given.
Each element ``i`` is associated with its weight ``\omega_i``.
In the following, for a function ``u: \Omega \rightarrow \mathbb{R}^{d^\prime}``,
we also denote by ``u \in \mathbb{R}^{nd^\prime}`` the finite element coefficient vector
corresponding to the global basis functions ``\Phi_k``.
This vector contains ``n`` blocks of length ``d^\prime`` denoting
all the components at each node ``k``.

To keep the spirit of the algorithm proposal and use julia's matrix operations efficiently,
we expand the FEM base by introducing
discrete derivative matrices for vector-valued functions.
Note that this might not be suitable for high performance computing,
but results in readable and understandable code, 
which is one of the primary intentions of this project.
Let
```math
D^{(j,r)}_{i,d^\prime (k-1)+r} = \frac{\partial}{\partial x_j} \Phi_k(x^{(i)})
```
where ``x^{(i)}`` denotes the midpoint of the element ``i``.
Multiplication to a coefficient vector results in a vector containing the partial
derivatives of the r-th component of a function in direction ``x_j``
at the midpoints of all elements.
Note that there are ``d \times d^\prime`` partial derivatives for all components
and directions at the midpoint of each element.
By construction, the vector has full length and all of the other entries are zero.
This requires slightly more memory, but allows for easy summation of all the derivatives.

## [Convex Problem](@id convex-problem-theory)
The first term of the objective functional can be understood as a norm and
can be discretized using the previously introduced derivative matrices by
```math
\Vert u+g \Vert_{X^p(\Omega)}^p 
    = \int_{\Omega} \Vert \nabla u+g \Vert_2^p \;\mathrm{d} x
    \approx \sum_{i=1}^m \omega_i
        \left(\sum_{j=1}^{d} \sum_{r=1}^{d^\prime} [D^{(j,r)}u + b]_i^2\right)^{\frac{p}{2}}
```
Note that for functions ``v \in W^{1,p}_g(\Omega,\R^d)`` it is possible to find
a decomposition ``v = u + g `` such that ``u \in W^{1,p}_0(\Omega,\R^d)`` and
``g \in W^{1,p}_g(\Omega,\R^d)`` fixed.
This means, we can consider a fixed prolongation of the Dirichlet condition to the
full domain, pre-compute ``b^{(j,r)} = D^{(j,r)}g`` and then just minimize over
``u \in W^{1,p}_0(\Omega,\R^d)`` with homogeneous Dirichlet conditions
in the relevant parts.

The second and the third term are standard integrals and can be discretized by
mass matrices ``M`` and ``\bar{M}`` on the domain and the boundary, respectively.
Now we can simply obtain a finite version of the problem by reformulating
the summation as a scalar product and introduce a slack variable for the non-linear term:  

*Minimizing the _p_-Laplace objective functional over the relevant finite element space*
*is equivalent to the standard convex problem*
```math
\min_{x \in \mathcal{Q}_p} \langle c,x \rangle \text{ with } c = 
    \begin{bmatrix}
        -Mf-\bar{M}h\\ 
        \frac{\omega}{p}
    \end{bmatrix}
```
*with a bounded and convex search set given by*
```math
\mathcal{Q}_p = \left\{ 
    (u,s) \in \mathbb{R}^n \times \mathbb{R}^m \, : \,
    s_i \geq \left(\sum_{j=1}^{d} \sum_{r=1}^{d^\prime} 
        [D^{(j,r)}(u+g)]_i^2\right)^{\frac{p}{2}}
    \land \omega_i s_i \leq R \right
\}.
```

In order to obtain the boundedness of the set with the second condition, an additional
constraint
```math
\omega_i \Vert \nabla (u+g)\vert_{K_i}\Vert_2^p \leq R
``` on the local gradients was added.
Here, the constant ``R`` can be chosen large enough to not influence the solution. 

## [Path-Following](@id path-following-theory)

The main part of the implementation is to efficiently solve this convex problem efficiently.
Thus, in \[[1](@ref references-theory)\] is proposed to apply an interior-point method
\[[4](@ref references-theory)\] which is well known in convex optimization.
In this method, the problem is solved by approximately following the central path
```math
x^*(t) = \underset{x \in \mathcal{Q}_p}{\mathrm{arg\;min}}\;
    t \langle c,x \rangle + F(x)
```
for ``t \rightarrow \infty`` and during that staying sufficiently close by satisfying
```math
\Vert tc + F^\prime \Vert^*_x \leq \beta.
```
Here ``\beta`` is a given parameter and ``\Vert \cdot \Vert_x^*`` denotes the norm induced
by ``F`` or in particular by its symmetric positive definite Hessian ``F^{\prime\prime}``
at ``x``.

Thus, in order to apply this method, we need a self-concordant barrier ``F`` of the 
constrained search set as well as its first and second derivatives.
We obtain this self-concordant barrier for ``\mathcal{Q}_p`` by 
```math
F_p(u,s) = -\sum_{i=1}^m \log{z_i} -\sum_{i=1}^m \log{s_i} -\sum_{i=1}^m \log{\tau_i}
```
where
```math
z_i(u,s) = s^{\frac{2}{p}} - \sum_{j=1}^{d} \sum_{r=1}^{d^\prime} [D^{(j,r)}u + b]_i^2 
    \qquad\text{and}\qquad \tau_i(s) = R - \omega_i s_i.
```
The derivatives can be found in the literature \[[4](@ref references-theory)\].

The general solution procedure then consists of two phases.
An auxiliary phase to compute a suitable starting point, and after that a main phase to
actually solve the constrained problem up to the given accuracy.
Details can be found in \[[1,3](@ref references-theory)\].
Here, we will only show the key steps for the simplest variant,
the short-step path-following.

For the **auxiliary path-following**, we start with ``t_0 = 1``, an
arbitrary ``x_0 \in \mathcal{Q}_p`` and set ``G =  -F^\prime(x_0)``.
Then we iterate over

1. ``t_{k+1} = t_k - \frac{\gamma}{\Vert G \Vert^*_{x_k}}``
2. ``x_{k+1} = x_k - [F^{\prime\prime}(x_k)]^{-1} (t_{k+1} G + F^\prime(x_k))``

until ``\Vert F^\prime(x_k) \Vert^*_{x_k} \leq \frac{\sqrt{\beta}}{1+\sqrt{\beta}}``.
Doing an additional update step, we then get an ``x \in \mathcal{Q}_p`` such that 
``\Vert F^\prime(x) \Vert^*_{x_k} \leq \beta``,
which is a necessary condition for the starting guess in the next phase.

In the **main path-following**, the iteration is quite similar. However, we start with
``t_0 = 0``, pick the final iterate from the previous phase as starting guess and use the
system vector ``c`` instead of the negative initial gradient.
The iteration then reads

1. ``t_{k+1} = t_k + \frac{\gamma}{\Vert c \Vert^*_{x_k}}``
2. ``x_{k+1} = x_k - [F^{\prime\prime}(x_k)]^{-1} (t_{k+1} c + F^\prime(x_k))``

until t exceeds some given (inverse) tolerance. Note that 2. now corresponds to a
classical Newton step for the central path.

In the implementation, by default, an adaptive stepping with the proposed parameters from
\[[1](@ref references-theory)\] is used.
For comparison or numerically challenging problems, it can be beneficial to use
the also implemented short- and long-step variants.
But be aware that they usually require many more iterations,
even if their theoretical worst-case estimate is better.
Further, you can adjust all the parameters like step length and its updates
by keyword arguments, if necessary.

## [References](@id references-theory)

\[1\] S. Loisel. “[*Efficient algorithms for solving the p-Laplacian in polynomial time*]
    (https://link.springer.com/article/10.1007/s00211-020-01141-z)”.
    Numerische Mathematik. 2020.\
\[2\] P.M. Müller, N. Kühl, M. Siebenborn, K. Deckelnick, M. Hinze and T. Rung.
    “[*A Novel p-Harmonic Descent Approach Applied to Fluid Dynamic Shape Optimization*]
    (https://link.springer.com/article/10.1007/s00158-021-03030-x)”.
    Structural and Multidisciplinary Optimization. 2021.\
\[3\] Y. Nesterov. “[*Introductory lectures on convex optimization*]
    (https://link.springer.com/book/10.1007/978-1-4419-8853-9)”.
    Applied Optimization. 2004.\
\[4\] H. Wyschka and M. Siebenborn. “[*Towards computing high-order p-harmonic descent
    directions and their limits in shape optimization*]
    (https://link.springer.com/chapter/10.1007/978-3-031-45158-4_9)”.
    In: "Modeling, Simulation and Optimization of Fluid Dynamic Applications". 2023.
