module Updates
    export SimpleMixing, ScheduledMixing, BroydenMixing, AndersonMixing

    using LinearAlgebra
    using ..FixedPointToolkit.Scheduler: QuadraticScheduling, ExponentialScheduling


    ##  Largest per-component |VOut - VIn| (works for scalars and vectors without allocating a temporary).
    max_residual(VIn, VOut) = maximum(((a, b),) -> abs(b - a), zip(VIn, VOut))


@doc """
```julia
SimpleMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, _extra...) --> Dict   where {T<:Function, N<:Union{Number, Vector{<:Number}}}
```

Returns a dictionary containing the following
- "VInNext" : Simple mixing of VIn and VOut with weight (1-`alpha`) and `alpha` respectively.
- "VOutNext" : The function `F` evaluated at VInNExt, with fixed args `F_args` and kwargs `F_kwargs`.
- "Delta" : The maximum absolute per-component difference b/w the next input and output, so convergence requires every component (not just the aggregate) to be within `tol`.
- "kwargs" : the kwargs of this function itself, to be used by a subsequent function call in the next fixed point iteration. Since this is simple mixing, `alpha` is held constant.

"""
    function SimpleMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, _extra...) :: Dict   where {T<:Function, N<:Union{Number, Vector{<:Number}}}

        VInNext     =    alpha * VOut + (1 - alpha) * VIn
        VOutNext    =    F(VInNext, F_args... ; F_kwargs...)

        return Dict("VInNext"   => VInNext,
                    "VOutNext"  => VOutNext,
                    "Delta"     => max_residual(VInNext, VOutNext),
                    "kwargs"    => Dict(:alpha => alpha))
    end


@doc """
```julia
ScheduledMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, iter::Int64 = 1, max_iter::Int64 = 1, Scheduler :: R = QuadraticScheduling, _extra...) --> Dict   where {T<:Function, N<:Union{Number, Vector{<:Number}}, R<:Function}
```

Returns a dictionary containing the following
- "VInNext" : Simple mixing of VIn and VOut with weight (1-beta) and beta respectively, where beta = `Scheduler(alpha)`, which may change with iteration.
- "VOutNext" : The function `F` evaluated at VInNExt, with fixed args `F_args` and kwargs `F_kwargs`.
- "Delta" : The maximum absolute per-component difference b/w the next input and output, so convergence requires every component (not just the aggregate) to be within `tol`.
- "kwargs" : the kwargs of this function itself, to be used by a subsequent function call in the next fixed point iteration. Here, the alpha value is updated every iteration using a Scheduler.

"""
    function ScheduledMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, iter::Int64 = 1, max_iter::Int64 = 1, Scheduler :: R = QuadraticScheduling, _extra...) :: Dict   where {T<:Function, N<:Union{Number, Vector{<:Number}}, R<:Function}

        ScheduledAlpha    =   Scheduler(alpha ; iter = iter, max_iter = max_iter)

        VInNext     =    ScheduledAlpha * VOut + (1 - ScheduledAlpha) * VIn
        VOutNext    =    F(VInNext, F_args... ; F_kwargs...)

        return Dict("VInNext"   => VInNext,
                    "VOutNext"  => VOutNext,
                    "Delta"     => max_residual(VInNext, VOutNext),
                    "kwargs"    => Dict(:alpha => ScheduledAlpha, :Scheduler => Scheduler))
    end


@doc """
```julia
BroydenMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, B::Z = alpha*Matrix(1.0I, length(VIn), length(VIn)), _extra...) --> Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}, Z<:Union{Number,Matrix{<:Number}}}
```

Quasi-Newton ("good" Broyden) update. With residual `g = VOut - VIn`, the next input is `VIn + B * g`,
where `B` approximates `-J⁻¹` for the Jacobian `J` of `g`. The initial `B = alpha * I` is exactly simple mixing with weight `alpha`;
afterwards `B` is improved by a rank-1 secant update every iteration (O(N²) work and memory).

Returns a dictionary containing the following
- "VInNext" : `VIn + B * (VOut - VIn)`.
- "VOutNext" : The function `F` evaluated at VInNExt, with fixed args `F_args` and kwargs `F_kwargs`.
- "Delta" : The maximum absolute per-component difference b/w the next input and output, so convergence requires every component (not just the aggregate) to be within `tol`.
- "kwargs" : the kwargs of this function itself, to be used by a subsequent function call in the next fixed point iteration. Here, the B-matrix is updated every iteration according to Broyden mixing rules.

!!! note
    Before v1.2 the step was `VIn - B * (VOut - VIn)`, so the default `B = alpha * I` initially stepped *away* from the fixed point.
    A `B` matrix saved in a checkpoint by an older version has the opposite sign convention.

"""
    function BroydenMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, B::Z = alpha*Matrix(1.0I, length(VIn), length(VIn)), _extra...) :: Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}, Z<:Union{Number,Matrix{<:Number}}}

        if N<:Number
            B       =   B[begin] ##This is to handle when the VIn is a float. We need to turn B to a float too
        end

        g           =   VOut - VIn
        VInNext     =   VIn + B * g
        VOutNext    =   F(VInNext, F_args... ; F_kwargs...)

        dVIn        =   VInNext - VIn
        dg          =   (VOutNext - VInNext) - g

        ##  Rank-1 update enforcing the secant condition BNext * dg = -dVIn, in O(N²):
        ##  BNext = B + (-dVIn - B dg) (B' dVIn)' / (dVIn' B dg)
        v           =   B' * dVIn
        denom       =   dot(v, dg)
        if isfinite(denom) && abs(denom) > sqrt(eps(Float64)) * norm(v) * norm(dg)
            u       =   -dVIn - B * dg
            BNext   =   B .+ u .* (v' ./ denom)
        else
            BNext   =   B       ##  degenerate step (e.g. already converged): keep the old B rather than producing NaNs
        end

        return Dict("VInNext"   => VInNext,
                    "VOutNext"  => VOutNext,
                    "Delta"     => max_residual(VInNext, VOutNext),
                    "kwargs"    => Dict(:alpha => alpha, :B => BNext))

    end


@doc """
```julia
AndersonMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict(), alpha::Float64 = 0.5, m::Int64 = 5, X::Vector = [], G::Vector = [], max_cond::Float64 = 1e10, restart::Float64 = 1.0, _extra...) --> Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}}
```

Anderson (Pulay / DIIS) mixing with history depth `m`. With residual `g = VOut - VIn`, the next input is

    VInNext = VIn + alpha * g - (ΔX + alpha * ΔG) * γ,     γ = argmin ‖g - ΔG γ‖,

where the columns of `ΔX` and `ΔG` are differences of the last `m + 1` inputs and residuals.
Costs O(N m²) per iteration and O(N m) memory, so it suits large inputs. `m = 0` is exactly [`SimpleMixing`](@ref).

Two safeguards:
- If the least-squares problem becomes ill-conditioned (condition number above `max_cond`), the oldest history entries are dropped.
- If the residual grows, `‖g‖ > restart * ‖g_previous‖`, the history is cleared and the step falls back to simple mixing.
  Stale history from a different nonlinear regime otherwise steers the extrapolation wrong, and a longer history makes it worse.
  Pass `restart = Inf` to disable.

Returns a dictionary containing the following
- "VInNext" : the Anderson-extrapolated next input.
- "VOutNext" : The function `F` evaluated at VInNext, with fixed args `F_args` and kwargs `F_kwargs`.
- "Delta" : The maximum absolute per-component difference b/w the next input and output.
- "kwargs" : the kwargs of this function itself, including the input history `X` and residual history `G` for the next iteration.

"""
    function AndersonMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict(), alpha::Float64 = 0.5, m::Int64 = 5, X::Vector = [], G::Vector = [], max_cond::Float64 = 1e10, restart::Float64 = 1.0, _extra...) :: Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}}

        @assert m >= 0 "History depth `m` must be non-negative."

        x       =   N<:Number ? [VIn] : VIn     ##  scalars are handled as 1-vectors internally
        g       =   N<:Number ? [VOut - VIn] : VOut - VIn

        ##  Residual grew: the history no longer describes the local behaviour of F, so start afresh.
        if !isempty(G) && norm(g) > restart * norm(G[end])
            X, G    =   empty(X), empty(G)
        end

        ##  Keep the last mEff+1 (input, residual) pairs; more than length(x) differences can never be independent.
        mEff    =   min(m, length(x))
        keep    =   max(0, length(X) - mEff)
        Xh      =   [X[keep+1:end]..., x]
        Gh      =   [G[keep+1:end]..., g]

        step    =   alpha * g
        while length(Xh) > 1
            ΔX  =   reduce(hcat, diff(Xh))
            ΔG  =   reduce(hcat, diff(Gh))
            QR  =   qr(ΔG)
            if cond(UpperTriangular(QR.R), 1) < max_cond
                γ       =   QR \ g
                step    =   step - (ΔX + alpha * ΔG) * γ
                break
            end
            popfirst!(Xh)       ##  ill-conditioned: forget the oldest pair and retry
            popfirst!(Gh)
        end

        xNext       =   x + step
        VInNext     =   N<:Number ? xNext[1] : xNext
        VOutNext    =   F(VInNext, F_args... ; F_kwargs...)

        return Dict("VInNext"   => VInNext,
                    "VOutNext"  => VOutNext,
                    "Delta"     => max_residual(VInNext, VOutNext),
                    "kwargs"    => Dict(:alpha => alpha, :m => m, :max_cond => max_cond, :restart => restart, :X => Xh, :G => Gh))
    end


end
