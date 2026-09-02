module Updates
    export SimpleMixing, ScheduledMixing, BroydenMixing

    using LinearAlgebra
    using ..FixedPointToolkit.Scheduler: QuadraticScheduling, ExponentialScheduling


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
                    "Delta"     => maximum(abs.(VOutNext - VInNext)),
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
                    "Delta"     => maximum(abs.(VOutNext - VInNext)),
                    "kwargs"    => Dict(:alpha => ScheduledAlpha, :Scheduler => Scheduler))
    end


@doc """
```julia
BroydenMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, B::Z = alpha*Matrix(1.0I, length(VIn), length(VIn)), _extra...) --> Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}, Z<:Union{Number,Matrix{<:Number}}}
```

Returns a dictionary containing the following
- "VInNext" : simple mixing of VIn and VOut with weight (1-beta) and beta respectively, where beta = `Scheduler(alpha)`, which may change with iteration. 
- "VOutNext" : The function `F` evaluated at VInNExt, with fixed args `F_args` and kwargs `F_kwargs`.
- "Delta" : The maximum absolute per-component difference b/w the next input and output, so convergence requires every component (not just the aggregate) to be within `tol`.
- "kwargs" : the kwargs of this function itself, to be used by a subsequent function call in the next fixed point iteration. Here, the B-matrix is updated every iteration according to Broyden mixing rules.

"""
    function BroydenMixing(VIn::N, VOut::N, F::T ; F_args::Tuple = (), F_kwargs::Dict = Dict() , alpha::Float64 = 0.5, B::Z = alpha*Matrix(1.0I, length(VIn), length(VIn)), _extra...) :: Dict where {T<:Function, N<:Union{Number, Vector{<:Number}}, Z<:Union{Number,Matrix{<:Number}}}
        
        if N<:Number
            B       =   B[begin] ##This is to handle when the VIn is a float. We need to turn B to a float too
        end

        VInNext     =   VIn - (B * (VOut - VIn))
        VOutNext    =   F(VInNext, F_args..., F_kwargs...)

        dVIn        =   VInNext - VIn
        dF          =   (VOutNext - VInNext) - (VOut - VIn)

        denom       =   transpose(dVIn) * B * dF
        BNext       =   B + ((((dVIn - B * dF) * transpose(dVIn)) * B) / (denom))

        return Dict("VInNext"   => VInNext,
                    "VOutNext"  => VOutNext,
                    "Delta"     => maximum(abs.(VOutNext - VInNext)),
                    "kwargs"    => Dict(:alpha => alpha, :B => BNext))

    end


end
