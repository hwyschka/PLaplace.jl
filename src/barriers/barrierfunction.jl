"""
    @enum Barrier

Variable used to specify the barrier function in the solution algorithm.
"""
@enum Barrier begin
    DEFAULT
end

struct BarrierFunction
    isadmissible::Function
    value::Function
    gradient::Function
    hessian::Function
    initialguess::Function
end

function FiniteBarrier()
    return BarrierFunction(
        isadmissible_finite,
        compute_value_finite,
        compute_gradient_finite,
        compute_hessian_finite,
        compute_initialguess_finite
    )
end

function BarrierFunction(barrier::Barrier, p::Float64)
    return FiniteBarrier()
end
