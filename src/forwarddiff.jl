function vec2scam(v::AbstractVector)
    Camera(SVector{3}(v[1:3]), SVector{3}(v[4:6]), v[7], SVector{2}(v[8:9]), SVector{2}(v[10:11]))
end

##################### objective #############################

function rodrigues_rotate_point(rot::AbstractVector{T}, X::AbstractVector{T}) where T
    sqtheta = sum(rot.*rot)
    if sqtheta > 1e-10
        theta = sqrt(sqtheta)
        costheta = cos(theta)
        sintheta = sin(theta)
        theta_inverse = 1. / theta

        w = theta_inverse * rot
        w_cross_X = cross(w,X)
        tmp = dot(w,X) * (1. - costheta)

        X*costheta + w_cross_X * sintheta + w * tmp
    else
        X + cross(rot,X)
    end
end

@inline function radial_distort(rad_params,proj)
    rsq = sum(proj.*proj)
    L = 1. + rad_params[1]*rsq + rad_params[2]*rsq*rsq
    proj*L
end

function project(cam, X)
    Xcam = rodrigues_rotate_point(cam.r, X - cam.c)
    distorted = radial_distort(cam.κ, SVector(Xcam[1], Xcam[2])/Xcam[3])
    distorted*cam.f + cam.x0
end

function compute_reproj_err(cam, X, w, feat)
    return w*(project(cam,X) - feat)
end

function pack(cam,X,w)
    vcat(cam.r, cam.c, cam.f, cam.x0, cam.κ, X, w)
end

function unpack(packed)
    N = length(packed)
    vec2scam(view(packed,1:N-4)), view(packed,N-3:N-1),packed[N]
end

function compute_w_err(w)
    1.0 - w*w
end

function compute_reproj_err_d(params, feat)
    cam, X, w = unpack(params)
    compute_reproj_err(cam,X,w,feat)
end

function compute_ba_J(::Val{:ForwardDiff}, cams::AbstractArray{T},X,w,obs,feats) where T
    p = size(obs,2)
    reproj_err_d = zeros(2*p, 11 + 3 + 1)
    for i in 1:p
        idx = (2*(i-1))+1
        @inbounds j = obs[1,i]
        @inbounds l = obs[2,i]
        @inbounds xi = SVector{15,T}(cams[:,j]...,X[:,l]...,w[i])
        @inbounds reproj_err_d[idx:idx+1,:] = ForwardDiff.jacobian(
            x -> compute_reproj_err_d(x, SVector{2,T}(feats[1,i], feats[2,i])),
              xi)
    end
    w_err_d = zeros(1,p)
    for i in 1:p
        @inbounds w_err_d[i] = ForwardDiff.derivative(compute_w_err, w[i])
    end
    (reproj_err_d, w_err_d)
end
