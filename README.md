<center>
    <img src="https://user-images.githubusercontent.com/66529462/184855093-33aca1eb-96a3-447f-adb6-42d2ac0a5a10.png" alt="fem" width="400"/>
</center>

## A Solver for Problems of the _p_-Laplacian PDE Operator.

[![][license-badge]][license-url]
[![][docs-badge]][docs-url]
[![][test-badge]][test-url]

This package provides a solver for problems of the scalar or vector-valued
[_p_-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) with finite _p_ including
source terms and mixed Dirichlet-Neumann boundary conditions.
The solver works iteratively based on a piece-wise linear finite element discretization
and interior-point methods.
Theory can be found in the publications \[[1](#references)\] and \[[3](#references)\]
or summarized in the documentation. 

The implementation is based on the finite element library
[MinFEM](https://github.com/msiebenborn/MinFEM.jl).
Thus, it is able to import meshes in GMSH v1, v2 and v4 format and outputs VTK format
for Paraview as well as statistics in plain txt.

Start by adding the PLaplace package to our julia installation and test it.
Note that the test might take some time depending on your machine since it solves a full
validation problem provided by the method of manufactured solutions \[[3](#references)\].
Thus, open the julia REPL, hit the **]** key and type

```
add PLaplace
test PLaplace
```

## An Example Problem

Lets go through a code for the _p_-Laplace equation on a unit square once with 
inhomogeneous Dirichlet boundary conditions and once with mixed homogenous Dirichlet
and inhomogeneous Neumann boundary conditions.

First we have to load the package PLaplace.
Here we also load the package MinFEM to gain access to it's mesh generation functions.
With these we then generate a uniform, triangular 30x30 mesh for the unit square.

```julia
using PLaplace, MinFEM

mesh = unit_square(30)
```

As an alternative: Download the package from github to obtain the examples and meshes
and navigate, within the julia console, to the **examples** folder.
Then import one of the mesh files generated with GMSH.
Here, the additional include of MinFEM is not necessary.

```julia
using PLaplace

mesh = import_mesh("meshes/square.msh")
```

As preparation, it is often useful to create or specify a directory for 
all the output files.

```julia
output_path::String = "results/readme-example/"
mkpath(output_path)
```

In particular, we can also prepare a file where we can later write statistics.
In this command, the `guarded` flag means that the file is not overwritten if it already
exist.
This way, we can keep this command in the script even if we already used it before and
just want to append new results, e.g., for different parameters to the file.

```julia
statistics_file = output_path * "statistics.txt"
write_statistics_header(statistics_file, guarded = true)
```

In the next step we specify the boundary conditions.
As the default procedure we simply specify julia functions as shown here for a
non-homogenous Dirichlet boundary condition. 

```julia
g(x) = x[1]^2
```
The same has to be done for mixed boundary conditions denoting one function for the
Dirichlet part and one for the Neumann part.
Note that also homogenous conditions need to be specified. 

```julia
g(x) = 0
h(x) = x[2]^2 - x[1]^2
```

Alternatively, the package also support providing all boundary conditions
and the source term as discrete values.
Either as FEM coefficients on the nodes or
values on quadrature points of the boundary elements.
Note that this procedure is only recommended for experienced users or
if an analytical description of the function is not available. 

```julia
g = evaluate_mesh_function(mesh, x -> x[1]^2)
```

Next, we need to specify the boundary sets for the respective conditions.
For pure Dirichlet problems it is easiest to select all physical boundaries of the mesh.

```julia
dirichlet_boundary = select_boundaries(mesh)
```

For mixed problems it is required to select the boundaries manually.
Note that the mesh in this example is designed to have four physical boundaries
identified by the indices 1001-1004.

```julia
dirichlet_boundary = select_boundaries(mesh, 1001, 1004)
neumann_boundary = select_boundaries(mesh, 1002, 1003)
```

Finally, we solve the problem and then write the solution in a file for visualization with
Paraview.
The object data of the custom type ```PLaplaceData``` contains further relevant data
generated during the iteration.
The most important parameters, such as iteration counts and obtained accuracy,
can now be stored to the previously created statistics file. 

```julia
p = 3.0
data = solve_plaplace(p, mesh, g, dirichlet_boundary)

write_statistics(statistics_file, data)

write_result_to_vtk(output_path * "result_p=$p", data)
```

While the solution routine always needs the parameter _p_, the mesh and the Dirichlet
boundary information, it can be provided with a number of optional keyword arguments.
The most prominent ones for Neumann boundary handling are shown here.
For a full list including solver parameters see the documentation.

```julia
data = solve_plaplace(
    p,
    mesh,
    g,
    dirichlet_boundary,
    h = h,
    neumann_boundary = neumann_boundary
)
```

Both of the presented example and many more are available in the examples directory for
experimenting.

## References

\[1\] S. Loisel. “[Efficient algorithms for solving the p-Laplacian in polynomial time](https://link.springer.com/article/10.1007/s00211-020-01141-z)”. Numerische Mathematik 146(2). 2020.\
\[2\] K. Salari and P. Knupp. “[Code Verification by the Method of Manufactured Solutions](https://www.osti.gov/biblio/759450-wLI4Ux/native/)”. Sandia Report. 2000.\
\[3\] H. Wyschka and M. Siebenborn.“[Towards computing high-order p-harmonic descent directions and their limits in shape optimization](https://link.springer.com/chapter/10.1007/978-3-031-45158-4_9)”. In: "Modeling, Simulation and Optimization of Fluid Dynamic Applications". 2023.

[license-url]: https://github.com/hwyschka/PLaplace.jl/blob/master/LICENSE
[license-badge]: https://img.shields.io/badge/License-MIT-brightgreen.svg
[docs-url]: https://hwyschka.github.io/PLaplace.jl/stable/
[docs-badge]: https://img.shields.io/badge/docs-stable-blue.svg
[test-url]: https://github.com/hwyschka/PLaplace.jl/actions/workflows/test.yml
[test-badge]: https://github.com/hwyschka/PLaplace.jl/actions/workflows/test.yml/badge.svg
[cov-url]: https://codecov.io/gh/hwyschka/PLaplace.jl
[cov-badge]: https://codecov.io/gh/hwyschka/PLaplace.jl/branch/master/graph/badge.svg
