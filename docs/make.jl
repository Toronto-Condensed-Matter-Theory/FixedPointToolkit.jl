using Documenter
using FixedPointToolkit

makedocs(
    build       =   "build" ,
    sitename    =   "FixedPointToolkit.jl"    ,
    modules     =   [FixedPointToolkit, FixedPointToolkit.Scheduler, FixedPointToolkit.Updates, FixedPointToolkit.SelfConsistency, FixedPointToolkit.Checkpointer, FixedPointToolkit.FPoint, 
                            FixedPointToolkit.PlotSelfCons]   ,
    pages = [
        "Introduction"              =>  "index.md",
        "Schedulers"                =>  "Scheduler.md",
        "Updates"                   =>  "Updates.md",
        "Self-Consistency"          =>  "SelfCons.md",
        "Checkpointing"             =>  "Checkpointer.md",
        "Fixed Point calculation"   =>  "FixedPoint.md",
        "Plotting Results"          =>  "PlotSelfCons.md",
    ]
)

deploydocs(
    repo = "github.com/Toronto-Condensed-Matter-Theory/FixedPointToolkit.jl.git",
    devbranch = "main"
)
