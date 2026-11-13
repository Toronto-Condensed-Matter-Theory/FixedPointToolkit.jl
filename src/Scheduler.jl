module Scheduler
    export QuadraticScheduling, ExponentialScheduling

@doc """
```julia
QuadraticScheduling(alpha::Float64; iter::Int64, max_iter::Int64) --> Float64
```
Returns an alpha value descreased quadratically with the iteration `iter`.

"""
    function QuadraticScheduling(alpha::Float64; iter::Int64, max_iter::Int64, rate::Float64 = 3.0) :: Float64

        return  alpha * ((rate^2 - 1) * (iter^2) + (max_iter^2 - rate^2)) / ((rate^2 - 1) * ((iter + 1)^2) + (max_iter^2 - rate^2))
    end


@doc """
```julia
ExponentialScheduling(alpha::Float64; iter::Int64, max_iter::Int64, rate::Float64 = 1.0) --> Float64
```
Returns an alpha value descreased exponentially with the iteration `iter` with some rate `rate`.

"""
    function ExponentialScheduling(alpha::Float64; iter::Int64, max_iter::Int64, rate::Float64 = 1.0) :: Float64

        return alpha * exp(- (rate / max_iter))
    end

end
