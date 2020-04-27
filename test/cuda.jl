using KernelAbstractions
using CuArrays
using NiBundleAdjustment
using NiBundleAdjustment: vec2cam, vec2scam, compute_reproj_err
using Test

function load(b, n, m, p)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "ba")
    fn = "ba$(b)_n$(n)_m$(m)_p$(p)"
    fn_in = joinpath(dir_in, fn)
    NiBundleAdjustment.read_ba_instance(string(fn_in,".txt"))
end

@testset "multi-threading and cuda" begin
    cams, X, w, obs, feats = load(4, 372, 47423, 204472)
    CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
    XX = [P3(X[:,i]...) for i=1:size(X, 2)]
    FEATS = [P2(feats[:,i]...) for i=1:size(feats, 2)]

    res1 = compute_ba_J(Val(:NiLang), copy(CAMS), copy(XX), copy(w), copy(obs), copy(FEATS))
    res2 = compute_ba_J_multithread(Val(:NiLang), copy(CAMS), copy(XX), copy(w), copy(obs), copy(FEATS); nthread=4)
    @test res1[1] ≈ res2[1]
    @test res1[2] ≈ res2[2]
end

@testset "multi-threading and cuda" begin
    cams, X, w, obs, feats = load(4, 372, 47423, 204472)
    CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
    XX = [P3(X[:,i]...) for i=1:size(X, 2)]
    FEATS = [P2(feats[:,i]...) for i=1:size(feats, 2)]

    res1 = compute_ba_J(Val(:NiLang), copy(CAMS), copy(XX), copy(w), copy(obs), copy(FEATS))
    res3 = compute_ba_J_cuda(Val(:NiLang), (CuArray(CAMS)), (CuArray(XX)), (CuArray(w)), (CuArray(obs)), (CuArray(FEATS)))
    @test res1[1] ≈ Array(res3[1])
    @test res1[2] ≈ Array(res3[2])
end
