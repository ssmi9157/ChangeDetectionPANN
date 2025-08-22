"""
    getS2params()

Returns a dictionary with the values for the Sentinel-2 decision tree.
"""
function getS2params()
    s2 = Dict(
              "firesBands"=> [5,6,7,8,9,12,13],
              "firesIdx"=> [13,9],
              "firesThres"=> 0.1,
              "floodsBands"=> [2,6,9,12],
              "floodsIdx"=> [13,2],
              "floodsThres"=> -0.2,
              "hurricanesBands"=> [2,3,4,13],
              "landslidesBands"=> [2,3,4,6,7,8,9],
              "landslidesIdx"=> [9,4],
              "landslidesThres"=> -0.35,
             )

    return s2
end

"""
    getL8params()

Returns a dictionary with the values for the LandSat-8 decision tree.
"""
function getL8params()
    l8 = Dict(
              "firesBands"=> [4,5,6,7],
              "firesIdx"=> [6,4],
              "firesThres"=> 0.1,
              "floodsBands"=> [1,2,3,4,5,8,9],
              "landslidesBands"=> [2,3,4,5,6,9],
              "landslidesIdx"=> [4,3],
              "landslidesThres"=> -0.55,
             )

    return l8
end

"""
    loadTIF(handle, bands::Vector{Int})

Returns the data of the .tif given by the handle.

# Arguments
- 'handle': The ArchGDAL handle of the .tif file.
- 'bands': The index of the bands in the .tif file to be loaded.
"""
function loadTIF(handle, bands::Vector{Int})
    dat = AG.read(handle, Int32.(bands))
    return permutedims(dat, [2,1,3])
end

"""
    getIndexScore(mode::String, handles)

Returns the index score for the given mode.

# Arguments
- 'mode': String used to determine index score to calculate.
- 'handles': Vector of ArchGDAL handles for the .tif files in the natural 
  disaster scene.
"""
function getIndexScore(mode::String, handles, params)
    bands = params[mode*"Idx"]
    threshold = params[mode*"Thres"]

    imgs = [loadTIF(h, bands) for h in handles]

    idxs = []
    for img in imgs
        idx = (img[:,:,1] .- img[:,:,2]) ./ (img[:,:,1] .+ img[:,:,2])
        pro = zeros(size(idx))
        pro[idx .> threshold] .= 1
        push!(idxs, pro)
    end

    final = abs.(idxs[end] .- idxs[end-1])
    replace!(final, NaN=>0)
    score = sum(final) / length(idxs[1])

    return score
end

"""
    s2Tree(handles)

Returns a set of bands as determined by the Sentinel-2 decison tree.

# Arguments
- 'handles': Vector of ArchGDAL handles for the .tif files in the natural 
  disaster scene.
"""
function s2Tree(handles)
    params = getS2params()
    bai = getIndexScore("fires", handles, params)
    if bai >= 0.05
        return params["firesBands"]
    end

    landScore = getIndexScore("landslides", handles, params)
    if landScore < 0.005
        return params["landslidesBands"]
    end

    floodScore = getIndexScore("floods", handles, params)
    if floodScore > 0.15
        return params["floodsBands"]
    else
        return params["hurricanesBands"]
    end
end

"""
    l8Tree(handles)

Returns a set of bands as determined by the LandSat-8 decison tree.

# Arguments
- 'handles': Vector of ArchGDAL handles for the .tif files in the natural 
  disaster scene.
"""
function l8Tree(handles)
    params = getL8params()
    bai = getIndexScore("fires", handles, params)
    if bai >= 0.05
        return params["firesBands"]
    end

    landScore = getIndexScore("landslides", handles, params)
    if landScore < 0.0001
        return params["landslidesBands"]
    end

    return params["floodsBands"]
end

"""
    bandSelection(path::String)

Returns the bands to be used for the scene based on the index scores.

# Arguments
- 'path': File path to the directory with the .tif files for the natural 
  disaster scene.
- 'sensor': The directory name that contains the data from that sensor.
"""
function bandSelection(path::String, sensor::String)
    files = readdir(path)[end-1:end]
    handles = [AG.read(joinpath(path, fi)) for fi in files]

    if sensor == "S2"
        return s2Tree(handles)
    elseif sensor == "L8"
        return l8Tree(handles)
    end

    return nothing
end
