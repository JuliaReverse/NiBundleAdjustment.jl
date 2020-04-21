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

@i @inline function idot(out!, x, y)
    @safe @assert length(x) == length(y)
    @invcheckoff @inbounds for i=1:length(x)
        out! += x[i] * y[i]
    end
end

@i @inline function idot(out!, x::P2, y::P2)
    out! += x.x * y.x
    out! += x.y * y.y
end

@i @inline function idot(out!, x::P3, y::P3)
    out! += x.x * y.x
    out! += x.y * y.y
    out! += x.z * y.z
end

@i @inline function inorm2(out!, x)
    @invcheckoff @inbounds for i=1:length(x)
        out! += abs2(x[i])
    end
end

@i @inline function inorm2(out!, x::P2)
    out! += abs2(x.x)
    out! += abs2(x.y)
end

@i @inline function inorm2(out!, x::P3)
    out! += abs2(x.x)
    out! += abs2(x.y)
    out! += abs2(x.z)
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

@i @inline function icross!(out!, x::P3, y::P3)
    @inbounds begin
        out!.x += x.y * y.z
        out!.x -= x.z * y.y
        out!.y += x.z * y.x
        out!.y -= x.x * y.z
        out!.z += x.x * y.y
        out!.z -= x.y * y.x
    end
end

@i @inline function iaxpy!(a, X, Y)
    @safe @assert length(X) == length(Y)
    @invcheckoff @inbounds for i=1:length(Y)
        Y[i] += a * X[i]
    end
end

@i @inline function iaxpy!(a, X::P2, Y::P2)
    Y.x += a * X.x
    Y.y += a * X.y
end

@i @inline function iaxpy!(a, X::P3, Y::P3)
    Y.x += a * X.x
    Y.y += a * X.y
    Y.z += a * X.z
end

@i @inline function rodrigues_rotate_point(out!::P3{T}, rot::P3{T}, X::P3{T}) where T
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
        out! += X + w_cross_X
    end
    (~inorm2)(sqtheta, rot)
    @invcheckoff sqtheta → zero(T)
end

@i @inline function radial_distort(out!::P2{T}, rad_params::P2{T}, proj::P2{T}) where T
    @routine @invcheckoff begin
        rsq ← zero(T)
        L ← zero(T)
        inorm2(rsq, proj)
        rsqsq ← zero(T)
        rsqsq += rsq ^ 2
        @inbounds L += rad_params.y * rsqsq
        @inbounds L += rad_params.x * rsq
        L += identity(1)
    end
    iaxpy!(L, proj, out!)
    ~@routine
end

@i @inline function ⊕(identity)(x::P2, y::P2)
    x.x += identity(y.x)
    x.y += identity(y.y)
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
    x.x += identity(y.x)
    x.y += identity(y.y)
    x.z += identity(y.z)
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
        X -= identity(cam.c)
        rodrigues_rotate_point(Xcam, cam.r, X)
        Xcam2.x += Xcam.x / Xcam.z
        Xcam2.y += Xcam.y / Xcam.z
        radial_distort(distorted, cam.κ, Xcam2)
    end
    iaxpy!(cam.f, distorted, out!)
    out! += identity(cam.x0)
    ~@routine
end

@i function compute_reproj_err(out!, diff, cam, X, w, feat)
    project(diff, cam, X)
    diff -= identity(feat)
    iaxpy!(w, diff, out!)
end

@i @inline function compute_w_err(out!, w)
    out! += identity(1.0)
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
        @inbounds w_err_d[i] = grad(compute_w_err'(Val(1), 0.0, w[i])[3])
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
