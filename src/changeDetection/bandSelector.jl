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
function getIndexScore(mode::String, handles)
    if mode == "fires"
        bands = [13,9]
        threshold = 0.1
    elseif mode == "landslides"
        bands = [9,4]
        threshold = -0.35
    elseif mode == "floods"
        bands = [13,2]
        threshold = -0.2
    else
        return nothing
    end

    imgs = [loadTIF(h, bands) for h in handles]

    idxs = []
    for img in imgs
        idx = (img[:,:,1] .- img[:,:,2]) ./ (img[:,:,1] .+ img[:,:,2])
        pro = zeros(size(idx))
	if mode == "floods"
            pro[idx .< threshold] .= 1
        else
            pro[idx .> threshold] .= 1
        end
        push!(idxs, pro)
    end

    final = abs.(idxs[end] .- idxs[end-1])
    replace!(final, NaN=>0)
    score = sum(final) / length(idxs[1])

    return score
end

"""
    bandSelection(path::String)

Returns the bands to be used for the scene based on the index scores.

# Arguments
- 'path': File path to the directory with the .tif files for the natural 
  disaster scene.
"""
function bandSelection(path::String)
    files = readdir(path)[end-1:end]
    handles = [AG.read(joinpath(path, fi)) for fi in files]

    bai = getIndexScore("fires", handles)
    if bai >= 0.05
        return([5,6,7,8,9,12,13])
    end

    landScore = getIndexScore("landslides", handles)
    if landScore < 0.005
        return([2,3,4,6,7,8,9])
    end

    floodScore = getIndexScore("floods", handles)
    if floodScore > 0.15
        return([2,6,9,12])
    else
        return([2,3,4,13])
    end

end
