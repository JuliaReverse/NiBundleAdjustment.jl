using CuArrays
using KernelAbstractions
CuArrays.allowscalar(false)

export compute_ba_J_cuda, compute_ba_J_multithread

println("NiBundleAdjustment: You find CUDA!")

@kernel function ba_kernel(reproj_err_d, obs, cams, X, w, feats)
    i = @index(Global)
    idx = (2*(i-1))+1
    @inbounds j = obs[1,i]
    @inbounds l = obs[2,i]
    @inbounds cam = cams[j]
    out! = zero(cam.x0)
    diff = zero(cam.x0)
    @inbounds res = compute_reproj_err(out!, diff, cam, X[l], w[i], feats[i])
    _, _, gcam, gX, gw, _ = (~compute_reproj_err)(P2(GVar(res[1].x, 1.0), GVar(res[1].y, 0.0)),
                    GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
    NiBundleAdjustment.packgrad!(reproj_err_d, idx, gcam, gX, gw)
    _, _, gcam, gX, gw, _ = (~compute_reproj_err)(P2(GVar(res[1].x, 0.0), GVar(res[1].y, 1.0)),
                    GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
    NiBundleAdjustment.packgrad!(reproj_err_d, idx+1, gcam, gX, gw)
end

@kernel function ba_kernel_w(w_err_d, w)
    i = @index(Global)
    @inbounds w_err_d[i] = grad(NiBundleAdjustment.compute_w_err'(Val(1), 0.0, w[i])[3])
end

function compute_ba_J_multithread(::Val{:NiLang}, cams::AbstractArray{<:Camera{T}}, X::AbstractArray{<:P3}, w, obs, feats::AbstractArray{<:P2}; nthread=4) where T
    p = size(obs,2)
    reproj_err_d = zeros(2*p, 15)
    event1 = ba_kernel(CPU(), nthread)(reproj_err_d, obs, cams, X, w, feats; ndrange=p)
    wait(event1)
    w_err_d = zeros(1,p)
    event2 = ba_kernel_w(CPU(), nthread)(w_err_d, w; ndrange=p)
    wait(event2)
    (reproj_err_d, w_err_d)
end

function compute_ba_J_cuda(::Val{:NiLang}, cams::CuArray{<:Camera{T}}, X::CuArray{<:P3}, w::CuArray, obs::CuArray, feats::CuArray{<:P2}; blockdim=256) where T
    p = size(obs,2)
    reproj_err_d = CuArrays.zeros(T, 2*p, 15)
    event1 = ba_kernel(CUDA(), blockdim)(reproj_err_d, obs, cams, X, w, feats; ndrange=p)
    wait(event1)
    w_err_d = CuArrays.zeros(T, 1, p)
    event2 = ba_kernel_w(CUDA(), blockdim)(w_err_d, w; ndrange=p)
    wait(event2)
    (reproj_err_d, w_err_d)
end
