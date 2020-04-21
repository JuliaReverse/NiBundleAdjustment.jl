using BenchmarkTools
using NiBundleAdjustment
using NiBundleAdjustment: vec2cam, vec2scam, compute_reproj_err
using StaticArrays

function load(b, n, m, p)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "ba")
    fn = "ba$(b)_n$(n)_m$(m)_p$(p)"
    fn_in = joinpath(dir_in, fn)
    NiBundleAdjustment.read_ba_instance(string(fn_in,".txt"))
end

cams, X, w, obs, feats = load(4, 372, 47423, 204472)
CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
XX = [P3(X[:,i]...) for i=1:size(X, 2)]
FEATS = [P2(feats[:,i]...) for i=1:size(feats, 2)]
println("Normal Objective")
SCAMS = [NiBundleAdjustment.vec2scam(cams[:,i]) for i = 1:size(cams,2)]
display(@benchmark NiBundleAdjustment.ba_objective($SCAMS, $X, $w, $obs, $feats))

println()
println("Reversible Objective")
reproj_err! = zeros(P2{Float64}, size(feats, 2))
w_err! = zero(w)
reproj_err_cache! = zeros(P2{Float64}, size(feats, 2))
display(@benchmark NiBundleAdjustment.ba_objective!($reproj_err!, $w_err!, $reproj_err_cache!, $CAMS, $XX, $w, $obs, $FEATS))

println("ForwardDiff Gradient")
display(@benchmark compute_ba_J(Val(:ForwardDiff), $cams, $X, $w, $obs, $feats))

println()
println("NiLang Gradient")
display(@benchmark compute_ba_J(Val(:NiLang), $CAMS, $XX, $w, $obs, $FEATS))
