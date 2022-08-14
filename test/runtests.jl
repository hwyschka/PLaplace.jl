using PLaplace
using Test

tests = ["utility_test", "scalar_test", "vector_test"]

for t in tests
    include("$(t).jl")
end
