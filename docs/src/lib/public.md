# Public Documentation

---

```@meta
CurrentModule = PLaplace
```

Documentation for `PLaplace.jl`'s public interface.

See the [Internal](internal.md) page of the library for the documentation 
of internal types and functions.

## Contents

```@contents
Pages = ["public.md"]
Depth = 3
```

## Module
```@docs
PLaplace
```

## Algorithm
```@docs
PLaplaceData
```

```@docs
solve_plaplace
```

### Keyword Types
```@docs
Stepsize
LinearSolver
Preconditioner
```

### Initial prolongations
```@docs
compute_prolongation_harmonic
compute_prolongation_zero
```

## FEM
```@docs
MinFEM.Mesh
MinFEM.Boundary
```

```@docs
MinFEM.import_mesh
MinFEM.select_boundaries
```

### Modified Derivatives
```@docs
assemble_derivativetensor
assemble_derivativetensor_boundary
assemble_derviativetensor_modified
compute_derivative
compute_normalderivative
```

## Statistics

### PLaplace Data
```@docs
hasresult
write_result_to_vtk
write_result_to_txt
print_defaultdata
print_statistics
write_statistics_header
check_statistics_header
write_statistics
```

### Log Data
```@docs
AlgorithmLogData
read_algorithmlog
```

## Error Analysis

```@docs
ErrorData
```

```@docs
objective_functional
compute_errors
```

```@docs
write_error_header
check_error_header
write_error
read_error
append!(A::ErrorData, B::ErrorData) 
```

## Other Utilities

```@docs
xpnorm
compute_lipschitzconstant_boundary
```
