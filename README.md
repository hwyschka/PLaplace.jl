<center><img src="https://user-images.githubusercontent.com/66529462/184855093-33aca1eb-96a3-447f-adb6-42d2ac0a5a10.png" alt="fem" width="400"/></center>

## A Solver for Problems of the _p_-Laplacian PDE Operator.

This package provides a solver for problems of the scalar or vector-valued [_p_-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) with finite _p_ including source terms and mixed Dirichlet-Neumann boundary conditions.
The solver works iteratively based on a piece-wise linear finite element discretization and interior-point methods.
Theory can be found in the publications \[[1](#references)\] and \[[2](#references)\]. 

The implementation is based on the finite element library [MinFEM](https://github.com/msiebenborn/MinFEM.jl).
Thus, it is able to import meshes in GMSH v1, v2 and v4 format and outputs VTK format for Paraview as well as statistics in plain txt.

Start by adding the PLaplace package to our julia installation and test it.
Note that the test might take some time depending on your machine since it solves a full validation problem provided by the method of manufactured solutions \[[3](#references)\].
Thus, open the julia REPL, hit the **]** key and type

```
add PLaplace
test PLaplace
```

## An Example Problem

Lets go through a code for the _p_-Laplace equation on a unit square once with inhomogeneous Dirichlet boundary conditions and once with mixed homogenous Dirichlet and inhomogeneous Neumann boundary conditions.

First we have to load the package PLaplace.
Here we also load the package MinFEM to gain access to it's mesh generation functions.
With these we then generate a uniform, triangular 30x30 mesh for the unit square.

```julia
using PLaplace, MinFEM

mesh = unit_square(30)
```

As an alternative: Download the package from github to obtain the examples and meshes and navigate, within the julia console, to the **examples** folder.
Then import one of the mesh files generated with GMSH.
Here the additional include of MinFEM is not necessary.

```julia
using PLaplace

mesh = import_mesh("../meshes/square.msh")
```

In the next step we specify the boundary conditions.
As the default procedure we simply specify julia functions as shown here for a non-homogenous dirichlet boundary condition. 

```julia
g(x) = x[1]^2
```
The same has to be done for mixed boundary conditions denoting one function for the Dirichlet part and one for the Neumann part. Note that also homogenous conditions need to be specified. 

```julia
g(x) = 0
h(x) = x[1]^1 - x[2]^2
```

Alternatively, the package also support providing all boundary conditions and the source term as discrete values.
Either as FEM coefficients on the nodes or values on quadrature points of the boundary elements.
Note that this procedure is only recommended for experienced users or if an analytical description of the function is not available. 

```julia
g = evaluate_mesh_function(mesh, x -> x[1]^2)
```

Next, we need to specify the boundary sets for the respective conditions.
For pure Dirichlet problems it is easiest to select all physical boundaries of the mesh.

```julia
dirichletBoundary = select_boundaries(mesh)
```

For mixed problems it is required to select the boundaries manually.
Note that the mesh in this example is designed to have four physical boundaries identified by the indices 1001-1004.

```julia
dirichletBoundary = select_boundaries(mesh, 1001, 1004)
neumannBoundary = select_boundaries(mesh, 1002, 1003)
```

Finally, we solve the problem and write the solution in a file for visualization with Paraview.
The object data of the custom type ```PLaplaceData``` contains further relevant data generated during the iteration.
The most important parameters, such as iteration counts and obtained accuracy, can be appended to a log file in .txt format. 

```julia
p = 3.0
data = solve_plaplace(p, mesh, g, dirichletBoundary)

mkpath("results")
write_result_to_vtk("results/example_$p", data)

mkpath("logs")
write_log("logs/example", data)
```

While the solution routine always needs the parameter _p_, the mesh and the Dirichlet boundary information, it can be provided with a number of optional keyword arguments.
The most prominent ones for Neumann boundary handling are shown here.
For a full list including solver parameters see the documentation.

```julia
data = solve_plaplace(p, mesh, g, dirichletBoundary, h=h, neumannBoundary=neumannBoundary)
```

## References

\[1\] S. Loisel. “[Efficient algorithms for solving the p-Laplacian in polynomial time](https://link.springer.com/article/10.1007/s00211-020-01141-z)”. Numerische
Mathematik 146.2, pp. 369–400. 2020.<br/>
\[2\] H. Wyschka and M. Siebenborn. “[Towards computing high-order p-harmonic descent directions and their limits in shape optimization](https://arxiv.org/abs/2208.06897)”. arXiv (Preprint). 2022<br/>
\[3\] K. Salari and P. Knupp. “[Code Verification by the Method of Manufactured Solutions](https://www.osti.gov/biblio/759450-wLI4Ux/native/)”. Sandia Report. 2000.
