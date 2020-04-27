module NiBundleAdjustment

using NiLang, NiLang.AD
using LinearAlgebra
using ForwardDiff
using StaticArrays
using Requires

export compute_ba_J, Camera, read_ba_instance

struct Camera{T,V3,V2}
    r::V3
    c::V3
    f::T
    x0::V2
    κ::V2
end

include("forwarddiff.jl")
include("reversible.jl")
include("load.jl")

function __init__()
    @require KernelAbstractions="63c18a36-062a-441e-b654-da1e3ab1ce7c" include("cuda.jl")
end

end # module
