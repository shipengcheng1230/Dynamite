module Dynamite

using Reexport

@reexport using GmshTools
using Combinatorics
using LinearAlgebra
using ReverseDiff
using TupleTools

export gmsh

include("element.jl")
include("basis.jl")
include("matrix.jl")
include("kernel.jl")
include("mesh.jl")

function __init__()
    try
        gmsh.initialize()
    catch ex
        error("GmshTools not installed properly.")
    end
end

end # module
