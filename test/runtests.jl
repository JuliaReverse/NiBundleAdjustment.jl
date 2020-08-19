using NiBundleAdjustment
using NiBundleAdjustment: i_cross!, vec2cam, compute_reproj_err,
    project, rand_camera, P2, P3
using NiLang, NiLang.AD
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
    @instr i_cross!(c, x, y)
    @test c ≈ cross(x, y)
end

@testset "check grad" begin
    cams, X, w, obs, feats = load(4, 372, 47423, 204472)
    XX = [P3(X[:,i]...) for i=1:size(X, 2)]
    FEATS = [P2(feats[:,i]...) for i=1:size(feats, 2)]
    CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
    j1, j2 = compute_ba_J(Val(:NiLang), CAMS, XX, w, obs, FEATS)
    @test j1[1] ≈ -461.44632100159936
    j1_, j2_ = compute_ba_J(Val(:ForwardDiff), cams, X, w, obs, feats)
    @test j1 ≈ j1_
    @test j2 ≈ j2_
end

@testset "project" begin
    cam = rand_camera()
    scam = Camera(SVector(cam.r), SVector(cam.c), cam.f, SVector(cam.x0), SVector(cam.κ))
    X = P3(randn(3)...)
    @test SVector(project(zero(P2{Float64}), cam, X)[1]) ≈ project(scam, SVector(X))
    @test check_inv(project,(zero(P2{Float64}), cam, X))

    w = randn()
    feat = MVector{2}(randn(2))
    cam = rand_camera()
    scam = Camera(SVector(cam.r), SVector(cam.c), cam.f, SVector(cam.x0), SVector(cam.κ))
    X = MVector{3}(randn(3))
    res = compute_reproj_err(zero(cam.x0), zero(cam.x0), cam, P3(X...), w, P2(feat...))[1]

    @test compute_reproj_err(scam, X, w, feat) ≈ SVector(res)
    @test check_inv(compute_reproj_err, (zero(cam.x0), zero(cam.x0), cam, P3(X...), w, P2(feat...)))
end

@testset "P2, P3" begin
    p = P2(2.0, 4.0)
    @test (~GVar)(GVar(p)) == p
    @test GVar(P2(2.0, 3.0), P2(1.0, 0.5)) == P2(GVar(2.0, 1.0), GVar(3.0, 0.5))
    p = P3(2.0, 4.0, 8.0)
    @test (~GVar)(GVar(p)) == p
    @test GVar(P3(2.0, 3.0, 4.0), P3(1.0, 0.5, 0.25)) == P3(GVar(2.0, 1.0), GVar(3.0, 0.5), GVar(4.0, 0.25))
end

@testset "objectives" begin
    cams, X, w, obs, feats = load(4, 372, 47423, 204472)
    XX = [P3(X[:,i]...) for i=1:size(X, 2)]
    FEATS = [P2(feats[:,i]...) for i=1:size(feats, 2)]
    SCAMS = [NiBundleAdjustment.vec2scam(cams[:,i]) for i = 1:size(cams,2)]
    CAMS = [vec2cam(cams[:,i]) for i = 1:size(cams,2)]
    reproj_err! = zeros(P2{Float64}, size(feats, 2))
    w_err! = zero(w)
    reproj_err_cache! = zeros(P2{Float64}, size(feats, 2))
    NiBundleAdjustment.ba_objective!(reproj_err!, w_err!, reproj_err_cache!, CAMS, XX, w, obs, FEATS)
    reproj_err, w_err = NiBundleAdjustment.ba_objective(SCAMS, X, w, obs, feats)
    correct = true
    for j=1:size(feats, 2)
        correct = correct && reproj_err[1,j] ≈ reproj_err![j].x
        correct = correct && reproj_err[2,j] ≈ reproj_err![j].y
    end
    @test correct
    @test w_err ≈ w_err!
end
