using NiBundleAdjustment
using NiBundleAdjustment: icross!, vec2cam, compute_reproj_err,
    project, rand_camera, mzeros
using NiLang
using LinearAlgebra, StaticArrays
using Test

function load(b, n, m, p)
    dir_in = joinpath(dirname(dirname(@__FILE__)), "data", "ba")
    fn = "ba$(b)_n$(n)_m$(m)_p$(p)"
    fn_in = joinpath(dir_in, fn)
    NiBundleAdjustment.read_ba_instance(string(fn_in,".txt"))
end

@testset "cross" begin
    x = randn(3)
    y = randn(3)
    c = zeros(3)
    @instr icross!(c, x, y)
    @test c ≈ cross(x, y)
end

@testset "check grad" begin
    cams, X, w, obs, feats = load(4, 372, 47423, 204472)
    CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
    j1, j2 = compute_ba_J(Val(:NiLang), CAMS, X, w, obs, feats)
    @test j1[1] ≈ -461.44632100159936
    j1_, j2_ = compute_ba_J(Val(:ForwardDiff), cams, X, w, obs, feats)
    @test j1 ≈ j1_
    @test j2 ≈ j2_
end

@testset "project" begin
    cam = rand_camera()
    X = MVector{3}(randn(3))
    @test project(mzeros(2), cam, X)[1] ≈ project(cam, X)
    @test check_inv(project,(mzeros(2), cam, X))

    w = randn()
    feat = MVector{2}(randn(2))
    cam = rand_camera()
    X = MVector{3}(randn(3))
    compute_reproj_err(mzeros(2), mzeros(2), cam, X, w, feat)[1]
    res = compute_reproj_err(zero(cam.x0), zero(cam.x0), cam, X, w, feat)[1]
    @test compute_reproj_err(cam, X, w, feat) ≈ res
    @test check_inv(compute_reproj_err, (zero(cam.x0), zero(cam.x0), cam, X, w, feat))
end
