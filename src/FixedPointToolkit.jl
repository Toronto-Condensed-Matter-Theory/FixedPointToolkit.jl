module FixedPointToolkit

    include("Scheduler.jl")
    using .Scheduler
    export QuadraticScheduling, ExponentialScheduling

    include("Updates.jl")
    using .Updates
    export SimpleMixing, ScheduledMixing, BroydenMixing, AndersonMixing

    include("SelfCons.jl")
    using .SelfConsistency
    export SelfCons

    include("Checkpointer.jl")
    using .Checkpointer
    export save_checkpoint, read_checkpoint, ReCreateSelfCons

    include("FixedPoint.jl")
    using .FPoint
    export FixedPoint!, ContinueFixedPoint!

    include("PlotSelfCons.jl")
    using .PlotSelfCons
    export Plot_History

    function __init__()
        ##  Plot_History is implemented in the Plots extension; point users at it if Plots is not loaded.
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f === Plot_History && Base.get_extension(@__MODULE__, :FixedPointToolkitPlotsExt) === nothing
                print(io, "\nPlot_History needs the Plots package: run `using Plots` first.")
            end
        end
    end

end # module FixedPointToolkit
