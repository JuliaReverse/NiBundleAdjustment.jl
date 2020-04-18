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
println("Normal Objective")
display(@benchmark compute_reproj_err($(vec2scam(cams[:,1])), $(X[:,1]), $(w[1]), $(feats[:,1])))

println()
println("Reversible Objective")
display(@benchmark compute_reproj_err($(P2(0.0, 0.0)), $(P2(0.0, 0.0)),
    $(vec2cam(cams[:,1])), $(P3(X[:,1]...)), $(w[1]), $(P2(feats[:,1]...))))
println()
println("NiLang Gradient")
display(@benchmark compute_ba_J(Val(:NiLang), $CAMS, $X, $w, $obs, $feats))

println("ForwardDiff Gradient")
display(@benchmark compute_ba_J(Val(:ForwardDiff), $cams, $X, $w, $obs, $feats))
