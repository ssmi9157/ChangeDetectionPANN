"""
    loadChangeMap(fileName::String)

Returns the data of the .csv file.

# Arguments
- 'fileName': The path to the .csv file to be loaded.
"""
function loadChangeMap(fileName::String)
    return readdlm(fileName, ',')
end

"""
    loadChangeMask(path::String)

Returns the change mask of target values.

# Arguments
- 'path': Path to the change mask to be loaded.
"""
function loadChangeMask(path::String)
    handle = AG.read(path)
    dat = AG.read(handle, [Int32(1)])
    return permutedims(dat, [2,1,3])[:,:,1]
end

"""
# Arguments
- 'window': Only load data within the window. The first to points define the
  left corner of the window and the last two points define how far beyond to
  load in the respective direction.
"""
function loadChangeMask(path::String, window::NTuple{4, Int32})
    handle = AG.read(path)
    buffer = Array{Float32}(undef, (window[4], window[3],1))
    AG.read!(handle, buffer, [Int32(1)], 
             window[1]:(window[1]+window[3]-Int32(1)),
             window[2]:(window[2]+window[4]-Int32(1)))
    return permutedims(buffer, [2,1,3])[:,:,1]
end

"""
    evaluateChangeMaps(expName::String, eventPath::String; 
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))

Returns the AUPRC for the natural disaster category given by 'eventPath'.

# Arguments
- 'expName': The name of the experiment that produced the change maps.
- 'eventPath': The path to the data for the natural disaster category.
- 'verbose': Determines if the auprc should be printed for the individual 
  natural disasters.
- 'tileSz': The size of the tile used to create the change maps.
"""
function evaluateChangeMaps(expName::String, eventPath::String; 
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))

    locations = readdir(eventPath)
    locationPaths = [joinpath(eventPath, loc, "changes") for loc in locations]

    totalChanges = []
    totalTargets = []
    for (i,loc) in enumerate(locations)

        targetPath = readdir(locationPaths[i])[1]
        targetPath = joinpath(locationPaths[i],targetPath)

        handle = AG.read(targetPath)
        imgSize = (AG.height(handle), AG.width(handle))
        outputSize = imgSize .- (imgSize .% tileSz)
	target = loadChangeMask(targetPath, 
                                Int32.((1,1,outputSize[1],outputSize[2])))
        push!(totalTargets, target)

        changeMap = expName*"_"*loc*"_changeMap.csv"
        changes = loadChangeMap(changeMap)
        auprc = calculateAUPRC(changes, target, tileSz)
        if verbose
            println("AUPRC for $loc: $auprc\n")
        end
        push!(totalChanges, changes)
    end

    auprc = calculateTotalAUPRC(totalChanges, totalTargets, tileSz)
    if verbose
        println("Overall AUPRC: $auprc.")
    end
    return auprc
end

"""
    evaluateALLChangeMaps(expName::String, basePath::String; 
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))

Returns the AUPRC for each natural disaster category.

# Arguments
- 'expName': The name of the experiment that produced the change maps.
- 'basePath': The path to the root directory that contains the data for each 
  natural disaster category.
- 'verbose': Determines if the auprc should be printed for the individual 
  natural disasters.
- 'tileSz': The size of the tile used to create the change maps.
"""
function evaluateAllChangeMaps(expName::String, basePath::String;
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))
    events = readdir(basePath)
    eventPaths = [joinpath(basePath, e) for e in events]
    for (e,eventPath) in enumerate(eventPaths)
        auprc = evaluateChangeMaps(expName, eventPath, verbose=verbose,
                                   tileSz=tileSz)
	println("AUPRC for event $(events[e]): $auprc\n")
    end
end

"""
function evaluateErrors(expNames::Vector{String}, eventPath::String;
                       verbose::Bool=false, tileSz::Tuple{Int, Int}= (32,32))

Returns the mean AUPRC for each experiment and the SEM for a single natural 
disaster category.

# Arguments
- 'expNames': The list of experiment names.
- 'eventPath': The path to the data for the natural disaster category.
- 'verbose': Determines if the auprc should be printed for the individual 
  natural disasters.
- 'tileSz': The size of the tile used to create the change maps.
"""
function evaluateErrors(expNames::Vector{String}, eventPath::String;
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))

    locations = readdir(eventPath)
    locationPaths = [joinpath(eventPath, loc, "changes") for loc in locations]

    targets = Dict{String, Matrix{Float64}}()
    for (i,loc) in enumerate(locations)

        targetPath = readdir(locationPaths[i])[1]
        targetPath = joinpath(locationPaths[i],targetPath)

        handle = AG.read(targetPath)
        imgSize = (AG.height(handle), AG.width(handle))
        outputSize = imgSize .- (imgSize .% tileSz)
        targets[loc] = loadChangeMask(targetPath,
                                Int32.((1,1,outputSize[1],outputSize[2])))
    end

    overallAuprc = []
    for expName in expNames

        totalChanges = []
        totalTargets = []
        for loc in locations

            push!(totalTargets, targets[loc])

            changeMap = expName*"_"*loc*"_changeMap.csv"
            changes = loadChangeMap(changeMap)
	    auprc = calculateAUPRC(changes, targets[loc], tileSz)
            push!(totalChanges, changes)
        end
        auprc = calculateTotalAUPRC(totalChanges, totalTargets, tileSz)
        push!(overallAuprc, auprc)
        if verbose
            println("AUPRC for $name: $auprc\n")
        end
    end

    result = mean(overallAuprc)
    sem = std(overallAuprc) / sqrt(length(expNames))
    if verbose
        println("Final AUPRC: $result +- $sem")
    end

    return result, sem
end

"""
function evaluateAllErrors(expNames::Vector{String}, basePath::String;
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))

Returns the mean AUPRC for each experiment and the SEM for natural disaster 
category.

# Arguments
- 'expNames': The list of experiment names.
- 'basePath': The path to the root directory that contains the data for each 
  natural disaster category.
- 'verbose': Determines if the auprc should be printed for the individual 
  natural disasters.
- 'tileSz': The size of the tile used to create the change maps.
"""
function evaluateAllErrors(expNames::Vector{String}, basePath::String;
                       verbose::Bool=false, tileSz::Tuple{Int, Int}=(32,32))
    events = readdir(basePath)
    eventPaths = [joinpath(basePath, e) for e in events]
    for (e,eventPath) in enumerate(eventPaths)
        result, sem = evaluateChangeMaps(expNames, eventPath, verbose=verbose, 
                                   tileSz=tileSz)
	println("Final AUPRC for event $(events[e]): $result +- $sem\n")
    end
end
