import NiLang: i_dot, i_norm2, i_axpy!
abstract type Point{T,N} end
export P2, P3

struct P2{T} <: Point{T,2}
    x::T
    y::T
end

struct P3{T} <: Point{T,3}
    x::T
    y::T
    z::T
end

Base.zero(::Type{P2{T}}) where T = P2(zero(T), zero(T))
Base.zero(::Type{P3{T}}) where T = P3(zero(T), zero(T), zero(T))
Base.zero(::P2{T}) where T = P2(zero(T), zero(T))
Base.zero(::P3{T}) where T = P3(zero(T), zero(T), zero(T))
NiLang.AD.GVar(x::P2) = P2(GVar(x.x), GVar(x.y))
(_::Type{Inv{GVar}})(x::P2) = P2((~GVar)(x.x), (~GVar)(x.y))
NiLang.AD.GVar(x::P2, g::P2) = P2(GVar(x.x, g.x), GVar(x.y, g.y))

NiLang.AD.GVar(x::P3) = P3(GVar(x.x), GVar(x.y), GVar(x.z))
(_::Type{Inv{GVar}})(x::P3) = P3((~GVar)(x.x), (~GVar)(x.y), (~GVar)(x.z))
NiLang.AD.GVar(x::P3, g::P3) = P3(GVar(x.x, g.x), GVar(x.y, g.y), GVar(x.z, g.z))

StaticArrays.SVector(p::P3) = SVector(p.x, p.y, p.z)
StaticArrays.SVector(p::P2) = SVector(p.x, p.y)

function vec2cam(v::AbstractVector)
    Camera(P3(v[1:3]...), P3(v[4:6]...), v[7], P2(v[8:9]...), P2(v[10:11]...))
end

rand_camera() = vec2cam(randn(11))
#mzeros(n) = MVector{n}(zeros(n))
NiLang.AD.GVar(cam::Camera) = Camera(GVar(cam.r), GVar(cam.c), GVar(cam.f), GVar(cam.x0), GVar(cam.κ))

##################### objective #############################

@i @inline function i_dot(out!, x::P2, y::P2)
    out! += x.x * y.x
    out! += x.y * y.y
end

@i @inline function i_dot(out!, x::P3, y::P3)
    out! += x.x * y.x
    out! += x.y * y.y
    out! += x.z * y.z
end

@i @inline function i_norm2(out!, x::P2)
    out! += abs2(x.x)
    out! += abs2(x.y)
end

@i @inline function i_norm2(out!, x::P3)
    out! += abs2(x.x)
    out! += abs2(x.y)
    out! += abs2(x.z)
end

@i @inline function i_cross!(out!, x, y)
    @inbounds begin
        out![1] += x[2] * y[3]
        out![1] -= x[3] * y[2]
        out![2] += x[3] * y[1]
        out![2] -= x[1] * y[3]
        out![3] += x[1] * y[2]
        out![3] -= x[2] * y[1]
    end
end

@i @inline function i_cross!(out!, x::P3, y::P3)
    @inbounds begin
        out!.x += x.y * y.z
        out!.x -= x.z * y.y
        out!.y += x.z * y.x
        out!.y -= x.x * y.z
        out!.z += x.x * y.y
        out!.z -= x.y * y.x
    end
end

@i @inline function i_axpy!(a, X::P2, Y::P2)
    Y.x += a * X.x
    Y.y += a * X.y
end

@i @inline function i_axpy!(a, X::P3, Y::P3)
    Y.x += a * X.x
    Y.y += a * X.y
    Y.z += a * X.z
end

@i @inline function rodrigues_rotate_point(out!::P3{T}, rot::P3{T}, X::P3{T}) where T
    @invcheckoff sqtheta ← zero(T)
    i_norm2(sqtheta, rot)
    @invcheckoff if (sqtheta > 1e-10, ~)
        @routine begin
            @zeros T tmp1 tmp2 wX s c theta invtheta
            @zeros P3{T} w w_cross_X
            theta += sqtheta ^ 0.5
            (s, c) += sincos(theta)
            tmp1 += 1 - c
            invtheta += 1 / theta
            i_axpy!(invtheta, rot, w)
            i_cross!(w_cross_X, w, X)
            i_dot(wX, w, X)
            tmp2 += wX * tmp1
        end
        i_axpy!(c, X, out!)
        i_axpy!(s, w_cross_X, out!)
        i_axpy!(tmp2, w, out!)
        ~@routine
    else
        i_cross!(w_cross_X, rot, X)
        out! += X + w_cross_X
    end
    (~i_norm2)(sqtheta, rot)
    @invcheckoff sqtheta → zero(T)
end

@i @inline function radial_distort(out!::P2{T}, rad_params::P2{T}, proj::P2{T}) where T
    @routine @invcheckoff begin
        @zeros T rsq L rsqsq
        i_norm2(rsq, proj)
        rsqsq += rsq ^ 2
        @inbounds L += rad_params.y * rsqsq
        @inbounds L += rad_params.x * rsq
        L += 1
    end
    i_axpy!(L, proj, out!)
    ~@routine
end

@i @inline function ⊕(identity)(x::P2, y::P2)
    x.x += y.x
    x.y += y.y
end

@i @inline function ⊕(+)(x::P2, y::P2, z::P2)
    x.x += y.x + z.x
    x.y += y.y + z.y
end

@i @inline function ⊕(-)(x::P2, y::P2, z::P2)
    x.x += y.x - z.x
    x.y += y.y - z.y
end

@i @inline function ⊕(identity)(x::P3, y::P3)
    x.x += y.x
    x.y += y.y
    x.z += y.z
end

@i @inline function ⊕(+)(x::P3, y::P3, z::P3)
    x.x += y.x + z.x
    x.y += y.y + z.y
    x.z += y.z + z.z
end

@i @inline function ⊕(-)(x::P3, y::P3, z::P3)
    x.x += y.x - z.x
    x.y += y.y - z.y
    x.z += y.z - z.z
end

@i function project(out!, cam::Camera{T}, X) where T
    @routine @invcheckoff begin
        Xcam ← zero(cam.r)
        Xcam2 ← zero(cam.x0)
        distorted ← zero(cam.x0)
        X -= cam.c
        rodrigues_rotate_point(Xcam, cam.r, X)
        Xcam2.x += Xcam.x / Xcam.z
        Xcam2.y += Xcam.y / Xcam.z
        radial_distort(distorted, cam.κ, Xcam2)
    end
    i_axpy!(cam.f, distorted, out!)
    out! += cam.x0
    ~@routine
end

@i function compute_reproj_err(out!, diff, cam, X, w, feat)
    project(diff, cam, X)
    diff -= feat
    i_axpy!(w, diff, out!)
end

@i @inline function compute_w_err(out!, w)
    out! += 1.0
    out! -= w ^ 2
end

function compute_ba_J(::Val{:NiLang}, cams::AbstractArray{<:Camera{T}}, X::AbstractArray{<:P3}, w, obs, feats::AbstractArray{<:P2}) where T
    p = size(obs,2)
    reproj_err_d = zeros(2*p, 15)
    for i in 1:p
        idx = (2*(i-1))+1
        @inbounds j = obs[1,i]
        @inbounds l = obs[2,i]
        cam = cams[j]
        out! = zero(cam.x0)
        diff = zero(cam.x0)
        res = compute_reproj_err(out!, diff, cam, X[l], w[i], feats[i])
        _, _, gcam, gX, gw, _ = (~compute_reproj_err)(P2(GVar(res[1].x, 1.0), GVar(res[1].y, 0.0)),
                        GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
        packgrad!(reproj_err_d, idx, gcam, gX, gw)
        _, _, gcam, gX, gw, _ = (~compute_reproj_err)(P2(GVar(res[1].x, 0.0), GVar(res[1].y, 1.0)),
                        GVar(res[2]), GVar(res[3]), GVar(res[4]), GVar(res[5]), GVar(res[6]))
        packgrad!(reproj_err_d, idx+1, gcam, gX, gw)
    end
    w_err_d = zeros(1,p)
    for i in 1:p
        @inbounds w_err_d[i] = grad(Grad(compute_w_err)(Val(1), 0.0, w[i])[3])
    end
    (reproj_err_d, w_err_d)
end

@inline function packgrad!(out, idx, gcam, gX, gw)
    @inbounds begin
        out[idx,1] = grad(gcam.r.x)
        out[idx,2] = grad(gcam.r.y)
        out[idx,3] = grad(gcam.r.z)
        out[idx,4] = grad(gcam.c.x)
        out[idx,5] = grad(gcam.c.y)
        out[idx,6] = grad(gcam.c.z)
        out[idx,7] = grad(gcam.f)
        out[idx,8] = grad(gcam.x0.x)
        out[idx,9] = grad(gcam.x0.y)
        out[idx,10] = grad(gcam.κ.x)
        out[idx,11] = grad(gcam.κ.y)
        out[idx,12] = grad(gX.x)
        out[idx,13] = grad(gX.y)
        out[idx,14] = grad(gX.z)
        out[idx,15] = grad(gw)
    end
end

@i function ba_objective!(reproj_err!, w_err!, reproj_err_cache!, cams, X, w, obs, feats)
    for i in 1:length(feats)
        compute_reproj_err(reproj_err![i], reproj_err_cache![i], cams[obs[1,i]], X[obs[2,i]], w[i], feats[i])
    end
    @inbounds for i=1:length(w_err!)
        w_err![i] += 1.0
        w_err![i] -= w[i] ^ 2
    end
end

