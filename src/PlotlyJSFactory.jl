module PlotlyJSFactory

using PlotlyJS, LinearAlgebra, SparseArrays
using LightGraphs, NetworkLayout
# package code goes here
export create_chord_plot
include("chord_plot.jl")

end # module
