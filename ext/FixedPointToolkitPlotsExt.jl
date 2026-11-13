module FixedPointToolkitPlotsExt

    using Plots, LinearAlgebra
    using FixedPointToolkit: SelfCons
    import FixedPointToolkit.PlotSelfCons: Plot_History


    function Plot_History(sc::SelfCons{T, R, S}, arg::String ; indices::Vector{Int64} = collect(1:length(sc.Initial)), plot_legend::Bool = true, log_plot::Bool = false) where {T<:Function, R<:Function, S<:Union{Number, Vector{<:Number}}}

        @assert arg in ["inputs", "outputs", "convergence"] "Unsupported arguement `arg` passed."

        plt     	=    plot(grid=false, legend = plot_legend, bg_legend = :transparent)

        if arg == "inputs"
            
            for dim in indices
                data    =   getindex.(sc.VIns, dim)
                plot!(data, marker = :circle, label = "input $(dim)", lw = 2.0)
            end

        elseif arg == "outputs"

            for dim in indices
                data    =   getindex.(sc.VOuts, dim)
                plot!(data, marker = :circle, label = "output $(dim)", lw = 2.0)
            end

        elseif arg == "convergence"

            data    =   norm.(sc.VOuts .- sc.VIns) / sqrt(length(sc.Initial))
            if log_plot
                data = log.(data)
            end

            plot!(data, marker = :circle, label = "convergence", lw = 2.0)
        end

        N       =   length(sc.VIns)
        ticks   =   collect(1 : max(N ÷ 10, 2) : N)

        if ticks[end] != N
            push!(ticks, N)
        end

        xticks!(ticks)
        xlabel!("Iteration", guidefontsize = 9)
        title!("Fixed point function $(arg)", titlefontsize = 12)

    end

end
