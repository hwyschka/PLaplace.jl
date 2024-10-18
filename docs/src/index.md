![PLaplace Logo](assets/logo.png)

# PLaplace.jl
*A Solver for Problems of the p-Laplacian PDE Operator.*

### Introduction

This package provides a solver for problems of the scalar or vector-valued
[_p_-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) with finite _p_ 
including source terms and mixed Dirichlet-Neumann boundary conditions.
The solver works iteratively based on a piece-wise linear finite element discretization
and interior-point methods.
Theory can be found in the publications \[[1](@ref References)\] and 
\[[3](@ref References)\] as well as a short description on the
[algorithm](@ref algorithm-theory) page. 

The implementation is based on the finite element library
[MinFEM.jl](https://minfem.github.io/MinFEM.jl/stable/).
Thus, it is able to import meshes in GMSH v1, v2 and v4 format
and outputs VTK format for Paraview as well as statistics in plain txt.

### Getting Started

Start by adding the PLaplace package to our julia installation and test it.
Note that the test might take some time depending on your machine, since it solves a full
validation problem provided by the method of
manufactured solutions \[[2](@ref References)\].
Thus, open the julia REPL, hit the **]** key and type

```
add PLaplace
test PLaplace
```

### Examples

The **Examples** section features some demonstrations for different types of problems
that can be solved using this package.
This only showcases parts of the functionalities, but gives an idea of the usage.
Original code files are in addition available in the 
[GitHub repository](https://github.com/hwyschka/PLapalce.jl/tree/master/examples).

```@contents
Pages = [
    "examples/dirichlet-square.md",
    "examples/neumann-square.md",
    "examples/vector-neumann-cube.md"
]
Depth = 1
```

### Library

For more experienced users, we further offer the **Library** with a full list of types,
functions and methods provided by PLaplace featuring explicit documentation.
This also includes extensions of the algorithm not used for the tutorials
showed in this documentation.

```@contents
Pages = ["lib/public.md", "lib/internal.md"]
Depth = 1
```

### References
\[1\] S. Loisel. “[*Efficient algorithms for solving the p-Laplacian in polynomial time*]
    (https://link.springer.com/article/10.1007/s00211-020-01141-z)”.
    Numerische Mathematik. 2020.\
\[2\] K. Salari and P. Knupp. “[*Code Verification by the Method of Manufactured Solutions*]
    (https://www.osti.gov/biblio/759450-wLI4Ux/native/)”. Sandia Report. 2000.\
\[3\] H. Wyschka and M. Siebenborn. “[*Towards computing high-order p-harmonic descent
    directions and their limits in shape optimization*]
    (https://link.springer.com/chapter/10.1007/978-3-031-45158-4_9)”.
    In: "Modeling, Simulation and Optimization of Fluid Dynamic Applications". 2023.