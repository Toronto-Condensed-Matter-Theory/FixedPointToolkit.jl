# FixedPointToolkit.Updates

```@autodocs
Modules = [FixedPointToolkit, FixedPointToolkit.Updates]
Private = false
Pages = ["Updates.jl"]
```

# Hot to write your own Update
If you want to design your own update function, you need to follow certain guidelines
* Define your function as follows (`F` is the function whose fixed point is being calculated)
```julia
function CustomUpdate(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict(), kwargs..., _extra...) :: Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}}
    ##### Get VInNext by doing something on VIn, VOut, and stuff in kwargs.
    VOutNext    =   F(VInNext, F_args... ; F_kwargs...)

    ##### If you want to include scheduling, you can now include a function which updates the kwargs themselves --> call then new_kwargs (MUST contain all the keys of kwargs, even if they are unchanged)

    return Dict("VInNext"   => VInNext, 
                "VOutNext"  => VOutNext,
                "Delta"     => norm(VOutNext - VInNext) / sqrt(length(VInNext)), 
                "kwargs"    => new_kwargs)
end
```
* Every other argument you pass to your `CustomUpdate` must be a keyword argument. 
* Unlike `CustomScheduler`, you will be able to pass whatever `kwargs` you want to `CustomUpdate` when you begin the Fixed Point iterations.. 


# Choosing an update method
In a real self-consistency loop the cost is dominated by calls to `F`, so what matters is how many iterations a method needs.
The table below (from `Sample/MixingComparison.jl`) gives the number of `F` calls needed to reach `max|F(x) - x| < 1e-8`
on standard test problems; "–" means not converged within 2000 iterations.

| Problem | Simple α=0.5 | Broyden α=0.5 | Anderson m=5 | Anderson m=10 |
|---|---|---|---|---|
| `cos(x)` (scalar) | 12 | 6 | 6 | 6 |
| 2D polynomial system | 49 | 10 | 8 | 8 |
| Chandrasekhar H-equation, c = 0.99 | 155 | 18 | 14 | 13 |
| Chandrasekhar H-equation, c = 0.9999 | 1078 | 34 | 23 | 24 |
| Jacobi iteration, 1D Poisson, N = 30 | – | 36 | 232 | 96 |
| Nonlinear map with spectrum in [-0.9, 0.97], N = 60 | 132 | 45 | 38 | 34 |
| Complex linear map, N = 40 | 112 | 45 | 71 | 62 |

* [`AndersonMixing`](@ref) is the recommended default when `F` is expensive. It stores only `m` past steps (O(N m) memory).
  Stiff, nearly linear problems benefit from a longer history `m`.
* [`BroydenMixing`](@ref) is robust and remembers all past steps, but it stores a dense N×N matrix, so it suits N up to a few thousand.
* [`SimpleMixing`](@ref) is enough when `F` is cheap and the iteration is already a good contraction.
