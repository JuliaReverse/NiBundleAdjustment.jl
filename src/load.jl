function read_ba_instance(fn)
    fid = open(fn)
    lines = readlines(fid)
    close(fid)
    line=split(lines[1]," ")
    n = parse(Int,line[1])
    m = parse(Int,line[2])
    p = parse(Int,line[3])
    off = 2

    one_cam = zeros(Float64,11,1)
    line=split(lines[off]," ")
    for i in 1:11
        one_cam[i] = parse(Float64,line[i])
    end
    cams = repeat(one_cam,1,n)
    off += 1

    one_X = zeros(Float64,3,1)
    line=split(lines[off]," ")
    for i in 1:3
        one_X[i] = parse(Float64,line[i])
    end
    X = repeat(one_X,1,m)
    off += 1

    one_w = parse(Float64,lines[off])
    w = repeat([one_w],1,p)
    off += 1

    one_feat = zeros(Float64,2,1)
    line=split(lines[off]," ")
    for i in 1:2
        one_feat[i] = parse(Float64,line[i])
    end
    feats = repeat(one_feat,1,p)

    camIdx = 1
    ptIdx = 1
    obs = zeros(Int,2,p)
    for i in 1:p
        obs[1,i] = camIdx
        obs[2,i] = ptIdx
        camIdx = (camIdx%n) + 1
        ptIdx = (ptIdx%m) + 1
    end

    (cams,X,w,obs,feats)
end
