module NiBundleAdjustment

using NiLang, NiLang.AD
using LinearAlgebra
using ForwardDiff
using StaticArrays

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

end # module
