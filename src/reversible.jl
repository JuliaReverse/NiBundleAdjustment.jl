function vec2cam(v::AbstractVector)
    Camera(MVector{3}(v[1:3]), MVector{3}(v[4:6]), v[7], MVector{2}(v[8:9]), MVector{2}(v[10:11]))
end

rand_camera() = vec2cam(randn(11))
mzeros(n) = MVector{n}(zeros(n))
NiLang.AD.GVar(cam::Camera) = Camera(GVar(cam.r), GVar(cam.c), GVar(cam.f), GVar(cam.x0), GVar(cam.κ))

##################### objective #############################

@i @inline function idot(out!, x, y)
    @safe @assert length(x) == length(y)
    @invcheckoff @inbounds for i=1:length(x)
        out! += x[i] * y[i]
    end
end

@i @inline function inorm2(out!, x)
    @invcheckoff @inbounds for i=1:length(x)
        out! += abs2(x[i])
    end
end

@i @inline function icross!(out!, x, y)
    @inbounds begin
        out![1] += x[2] * y[3]
        out![1] -= x[3] * y[2]
        out![2] += x[3] * y[1]
        out![2] -= x[1] * y[3]
        out![3] += x[1] * y[2]
        out![3] -= x[2] * y[1]
    end
end

@i @inline function iaxpy!(a, X, Y)
    @safe @assert length(X) == length(Y)
    @invcheckoff @inbounds for i=1:length(Y)
        Y[i] += a * X[i]
    end
end

@i @inline function rodrigues_rotate_point(out!::AbstractVector{T}, rot::AbstractVector{T}, X::AbstractVector{T}) where T
    @invcheckoff sqtheta ← zero(T)
    inorm2(sqtheta, rot)
    @invcheckoff if (sqtheta > 1e-10, ~)
        @routine begin
            tmp1 ← zero(T)
            tmp2 ← zero(T)
            wX ← zero(T)
            w ← zero(out!)
            w_cross_X ← zero(out!)
            invtheta ← zero(sqtheta)
            s ← zero(T)
            c ← zero(T)
            theta ← zero(T)
            theta += sqtheta ^ 0.5
            (s, c) += sincos(theta)
            tmp1 += 1 - c
            invtheta += 1 / theta
            iaxpy!(invtheta, rot, w)
            icross!(w_cross_X, w, X)
            idot(wX, w, X)
            tmp2 += wX * tmp1
        end
        iaxpy!(c, X, out!)
        iaxpy!(s, w_cross_X, out!)
        iaxpy!(tmp2, w, out!)
        ~@routine
    else
        icross!(w_cross_X, rot, X)
        out! .+= X .+ w_cross_X
    end
    (~inorm2)(sqtheta, rot)
    @invcheckoff sqtheta → zero(T)
end

@i @inline function radial_distort(out!::AbstractVector{T}, rad_params, proj) where T
    @routine @invcheckoff begin
        rsq ← zero(T)
        L ← zero(T)
        inorm2(rsq, proj)
        rsqsq ← zero(T)
        rsqsq += rsq ^ 2
        @inbounds L += rad_params[2] * rsqsq
        @inbounds L += rad_params[1] * rsq
        L += identity(1)
    end
    iaxpy!(L, proj, out!)
    ~@routine
end

@i function loss1(l, out!::AbstractVector{T}, rad_params, proj) where T
    radial_distort(out!, rad_params, proj)
    l += identity(out![1])
end

@i function loss2(l, out!::AbstractVector{T}, rad_params, proj) where T
    rodrigues_rotate_point(out!, rad_params, proj)
    l += identity(out![1])
end

@i function project(out!, cam::Camera{T}, X) where T
    @routine @invcheckoff begin
        Xcam ← zero(cam.r)
        Xcam2 ← zero(cam.x0)
        distorted ← zero(cam.x0)
        X .-= identity.(cam.c)
        rodrigues_rotate_point(Xcam, cam.r, X)
        @inbounds Xcam2[1] += Xcam[1] / Xcam[3]
        @inbounds Xcam2[2] += Xcam[2] / Xcam[3]
        radial_distort(distorted, cam.κ, Xcam2)
    end
    iaxpy!(cam.f, distorted, out!)
    out! .+= identity.(cam.x0)
    ~@routine
end

@i function compute_reproj_err(out!, diff, cam, X, w, feat)
    project(diff, cam, X)
    diff .-= identity.(feat)
    iaxpy!(w, diff, out!)
end

@i @inline function compute_w_err(out!, w)
    out! += identity(1.0)
    out! -= w ^ 2
end

function compute_ba_J(::Val{:NiLang}, cams::AbstractArray{<:Camera{T}}, X, w, obs, feats) where T
    p = size(obs,2)
    reproj_err_d = zeros(2*p, 15)
    for i in 1:p
        idx = (2*(i-1))+1
        @inbounds j = obs[1,i]
        @inbounds l = obs[2,i]
        cam = cams[j]
        out! = zero(cam.x0)
        diff = zero(cam.x0)
        res = compute_reproj_err(out!, diff, cam, MVector(X[1,l], X[2,l], X[3,l]), w[i], MVector(feats[1,i], feats[2,i]))
        _, _, gcam, gX, gw, _ = (~compute_reproj_err)(MVector(GVar(res[1][1], 1.0), GVar(res[1][2], 0.0)),
                        GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
        packgrad!(reproj_err_d, idx, gcam, gX, gw)
        _, _, gcam, gX, gw, _ = (~compute_reproj_err)(MVector(GVar(res[1][1], 0.0), GVar(res[1][2], 1.0)),
                        GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
        packgrad!(reproj_err_d, idx+1, gcam, gX, gw)
    end
    w_err_d = zeros(1,p)
    for i in 1:p
        @inbounds w_err_d[i] = grad(compute_w_err'(Val(1), 0.0, w[i])[3])
    end
    (reproj_err_d, w_err_d)
end

@inline function packgrad!(out, idx, gcam, gX, gw)
    @inbounds out[idx,1:3] = grad(gcam.r)
    @inbounds out[idx,4:6] = grad(gcam.c)
    @inbounds out[idx,7] = grad(gcam.f)
    @inbounds out[idx,8:9] = grad(gcam.x0) # ok
    @inbounds out[idx,10:11] = grad(gcam.κ)
    @inbounds out[idx,12:14] = grad(gX)
    @inbounds out[idx,15] = grad(gw)  # ok
end

#=
@i @inline function rodrigues_rotate_point(out!::AbstractVector{T}, rot::AbstractVector{T}, X::AbstractVector{T}, wX::T, sqtheta::T, w_cross_X::AbstractVector{T}, w::AbstractVector{T}) where T
    idot(sqtheta, rot, rot)
    @invcheckoff if (sqtheta > 1e-10, ~)
        @routine begin
            tmp1 ← zero(T)
            tmp2 ← zero(T)
            invtheta ← zero(sqtheta)
            s ← zero(T)
            c ← zero(T)
            theta ← zero(T)
            theta += sqtheta ^ 0.5
            (s, c) += sincos(theta)
            tmp1 += 1 - c
            invtheta += 1 / theta
        end
        iaxpy!(invtheta, rot, w)
        icross!(w_cross_X, w, X)
        idot(wX, w, X)
        iaxpy!(c, X, out!)
        iaxpy!(s, w_cross_X, out!)
        tmp2 += wX * tmp1
        iaxpy!(tmp2, w, out!)
        tmp2 -= wX * tmp1
        ~@routine
    else
        icross!(w_cross_X, rot, X)
        out! .+= X .+ w_cross_X
    end
end

@i @inline function radial_distort(out!, L::T, rsq::T, rad_params, proj) where T
    idot(rsq, proj, proj)
    @routine @invcheckoff begin
        rsqsq ← zero(T)
        rsqsq += rsq ^ 2
    end
    @inbounds L += rad_params[2] * rsqsq
    ~@routine
    @inbounds L += rad_params[1] * rsq
    L += identity(1)
    iaxpy!(L, proj, out!)
end

@i function project(out!, cam::Camera{T}, X) where T
    @routine @invcheckoff begin
        Xcam ← zero(cam.r)
        w_cross_X ← zero(cam.r)
        w ← zero(cam.r)
        Xcam2 ← zero(cam.x0)
        distorted ← zero(cam.x0)
        rsq ← zero(T)
        L ← zero(T)
        wX ← zero(T)
        sqtheta ← zero(T)
        X .-= identity.(cam.c)
        rodrigues_rotate_point(Xcam, cam.r, X, wX, sqtheta, w_cross_X, w)
        @inbounds Xcam2[1] += Xcam[1] / Xcam[3]
        @inbounds Xcam2[2] += Xcam[2] / Xcam[3]
        radial_distort(distorted, L, rsq, cam.κ, Xcam2)
    end
    iaxpy!(cam.f, distorted, out!)
    out! .+= identity.(cam.x0)
    ~@routine
end
=#
