"""
    BarrierFunction

Object holdigng the relevant functions for evaluating a barrier function.
Includes gradient and hessian as well as a specific initial guess.

Is not intended to be initialized directly but via sub-constructors corresponding
to specific barrier files. See their documentations belov to identify the behaviour
of specific barrier implemtations. 

# Fields
$(TYPEDFIELDS)
"""
struct BarrierFunction
    "Returns boolean value wether ``x`` is in the corresponding admissible set."
    isadmissible::Function
    
    "Returns float value of the barrier function at ``x``."
    value::Function
    
    "Returns vector-valued gradient of the barrier function at ``x``."
    gradient::Function
    
    "Returns matrix-valued Hessian of the barrier function at ``x``."
    hessian::Function
    
    "Returns inital vector guess for starting an auxilliary path-following."
    initialguess::Function
end

"""
$(TYPEDSIGNATURES)

Returns BarrierFunction object corresponding to default barrier function for finite p.

In this setting the admissible set is given by

```math 
    \\mathcal{Q}_p 
        = \\left\\{ (u,s) \\in \\mathbb{R}^n \\times \\mathbb{R}^m \\, : \\,
        s_i \\geq \\left(\\sum_{j=1}^{d} \\sum_{r=1}^{d'}
        [D^{(j,r)}(u+g)]_i^2 \\right)^{\\frac{p}{2}}
        \\land \\omega_i s_i \\leq R \\right\\}
```

Defining the variables
```math 
    z_i(u,s) = s^{\\frac{2}{p}} - \\sum\\limits_{j=1}^d \\sum\\limits_{r=1}^d
        [\\underbrace{D^{(j,r)}(u+g)}_{=: y^{(j,r)}(u)}]_i^2
        \\quad \\text{and} \\quad \\tau_i(s) = R -\\omega_i s_i
```
it can be equipped with the self-concordant barrier
```math 
    F(x) = F(u,s) = -\\sum\\limits_{i=1}^{m} \\log(z_i)
        -\\sum\\limits_{i=1}^{m} \\log(s_i)
        -\\sum\\limits_{i=1}^{m} \\log(\\tau_i).
```

The derivatives
```math 
    F'(x) = \\begin{bmatrix} F_u \\\\ F_s \\end{bmatrix}
        \\text{ and }
        F''(x) = \\begin{bmatrix} 
            F_{uu} & F_{us} \\\\
            F_{us}^{\\intercal} & F_{ss}
        \\end{bmatrix}
```
are then given by
```math 
    \\begin{aligned}
        F_{u} &=2 \\sum_{j=1}^{d} \\sum_{r=1}^{d'} [D^{(j,r)}]^{\\intercal}
            \\frac{y^{(j,r)}}{z} \\\\
        F_{s} &= -\\frac{2}{p}\\frac{1}{z}s^{\\frac{2}{p}-1} + \\frac{\\omega}{\\tau}\\\\
        F_{uu} &= 2 \\sum_{j=1}^{d} \\sum_{r=1}^{d'} [D^{(j,r)}]^{\\intercal} Z^{-1}
            D^{(j,r)} \\\\
        &\\quad+ 4 \\sum_{j_1=1}^{d} \\sum_{r_1=1}^{d'}\\sum_{j_2=1}^{d}
            \\sum_{r_2=1}^{d'} (Y^{(j_1,r_1)}D^{(j_1,r_1)})^{\\intercal} Z^{-1} D^{(j)}
            (Y^{(j_2,r_2)}D^{(j_2,r_2)}),\\\\
        F_{us} &= -\\frac{4}{p}\\sum_{j=1}^{d} \\sum_{r=1}^{d'}
            (Y^{(j,r)}D^{(j,r)})^{\\intercal} Z^{-2} S^{\\frac{2}{p}-1},\\\\
        F_{ss} &= -\\frac{2}{p} \\left(\\frac{2}{p} - 1 \\right) Z^{-1} S^{\\frac{2}{p}-2}
            + \\frac{4}{p^2} Z^{-2} S^{\\frac{4}{p}-2} + W^2 T^{-2},\\\\
        \\text{where } S &= \\mathrm{diag}(s),\\, W = \\mathrm{diag}(\\omega),\\,
            Y = \\mathrm{diag}(y),\\, Z = \\mathrm{diag}(z)
            \\text{ and } T = \\mathrm{diag}(\\tau).
    \\end{aligned}
```

An initial guess inside ``\\mathcal{Q}_{\\infty}`` is given by
```math
    \\hat{x} = \\begin{bmatrix} 
            0 \\\\
            \\hat{s}
        \\end{bmatrix}
        \\in \\mathbb{R}^{n \\times m}
```
where
```math
    \\hat{s} = 1 + \\left(\\sum_{j=1}^{d} \\sum_{r=1}^{d'}
        [D^{(j,r)} g]_i^2 \\right)^{\\frac{p}{2}}
        \\in \\mathbb{R}^m.
```
"""
function FiniteBarrier()
    return BarrierFunction(
        isadmissible_finite,
        compute_value_finite,
        compute_gradient_finite,
        compute_hessian_finite,
        compute_initialguess_finite
    )
end

"""
$(TYPEDSIGNATURES)

Wrapping constructor for the selection of a barrier function.

For the specific behaviour of a barrier see
the documentation of the sub-constructors below.
"""
function BarrierFunction() :: BarrierFunction
    return FiniteBarrier()
end
