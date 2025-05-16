"""
    getPathsAndLocations(dataPath::String, seed::Int)

Returns the paths to .tif files for all the natural disasters in a random order.

# Arguments
- 'dataPath': Path to the root directory containing all the data.
- 'seed': Seed used to randomise the order of the paths.
"""
function getPathsAndLocations(dataPath::String, seed::Int)
    events = readdir(dataPath)
    events = [joinpath(dataPath, e) for e in events]
    paths = []
    locations = []
    for e in events
        location = readdir(e)
        append!(paths, joinpath(e, loc, "S2") for loc in location)
        append!(locations, location)
    end
    Random.seed!(seed)
    paths = shuffle(paths)
    Random.seed!(seed)
    locations = shuffle(locations)

    return paths, locations
end

"""
    createDataLoaders(bands::Vector{Int}, minVal::Real, maxVal::Real,
                      history::Int, tileSz::Tuple{Int, Int}, 
                      poolSz::Int, stride::Int)

Returns a list of DataLoaders for each band.

# Arguments
- 'bands': Bands to create a DataLoader for.
"""
function createDataLoaders(bands::Vector{Int}, minVal::Real, maxVal::Real,
                           history::Int, tileSz::Tuple{Int, Int}, 
                           poolSz::Int, stride::Int)
    dataLoaders = Dict{Int32, DataLoader}()
    for b in bands
        dataLoaders[b] = createDL(minVal, maxVal, b, tileSz, history, 
                                  poolSz, stride)
    end
    return dataLoaders
end

"""
    createPannModels(bands::Vector{Int}, connectivityFile::String, 
              inputNodes::Vector{Int}, numReadoutNodes::Int, T::Real, dt::Real,
              seed::Int)

Returns a list of PannModels for each band.

# Arguments
- 'bands': Bands to create a PannModel for.
"""
function createPannModels(bands::Vector{Int}, connectivityFile::String, 
              inputNodes::Vector{Int}, numReadoutNodes::Int, T::Real, dt::Real,
              seed::Int)
    panns = Dict{Int32, PannModel}()
    edgeParams = createEdgeParams()
    for b in bands
        panns[b] = createPannModel(connectivityFile,inputNodes,numReadoutNodes, 
                              T, dt, edgeParams, seed=seed)
    end
    return panns
end

"""
    getInputNodes(gridSz::Int, numberOfNodes::Int)

Returns a list of nodes to be used as inputs based on the grid size.

# Arguments
- 'gridSz': Number of input nodes per row. The grid is square.
- 'numberOfNodes': Number of non-input nodes in the networks.
"""
function getInputNodes(gridSz::Int, numberOfNodes::Int)
    numInputNodes = gridSz^2
    inputNodes = reshape(collect(1:numInputNodes), (gridSz,gridSz))
    inputNodes = reverse(inputNodes, dims=1)
    inputNodes = vec(permutedims(inputNodes, (2,1)))
    inputNodes .+= numberOfNodes
    return inputNodes
end

"""
     runScenes(connectivityFile::String, gridSz::Int, numberOfNodes::Int, 
               numReadoutNodes::Int, name::String, dataPath::String; 
               tileSz::Tuple{Int, Int}=(32,32), 
               history::Int=4, saveFiles::Bool=false,
               enableProgBar::Bool=true, distanceMetric=CorrDist(),
               bands::Vector{Int}=[2,3,4,5,6,7,8,9,12,13], seed::Int=1, 
               minVal::Real=-1, maxVal::Real=1,
               poolSz::Int=2, stride::Int=2)

Runs the desired experiment with the given variables.

# Arguments
- 'connectivityFile': Name of the connectivity file that contains the network 
  structure information.
- 'gridSz': Number of input nodes per row. The grid is square.
- 'numberOfNodes': Number of non-input nodes in the networks.
- 'numReadoutNodes': Number of nodes to be used as readouts for each network.
- 'name': Name of the experiment being run.
- 'dataPath': Path to the root directory containing all the data.
"""
function runScenes(connectivityFile::String, gridSz::Int, numberOfNodes::Int, 
                   numReadoutNodes::Int, name::String, dataPath::String;
                   tileSz::Tuple{Int, Int}=(32,32), 
                   history::Int=4, saveFiles::Bool=false,
                   enableProgBar::Bool=true, distanceMetric=CorrDist(),
                   bands::Vector{Int}=[2,3,4,5,6,7,8,9,12,13], seed::Int=1, 
                   minVal::Real=-1, maxVal::Real=1,
                   poolSz::Int=2, stride::Int=2)

    dt = 0.005
    T = (history+1) * dt
    inputNodes = getInputNodes(gridSz, numberOfNodes)
    panns = createPannModels(bands, connectivityFile, inputNodes, 
                             numReadoutNodes, T, dt, seed)
    dataLoaders = createDataLoaders(bands, minVal, maxVal, history, tileSz, 
                                    poolSz, stride)

    paths, locations = getPathsAndLocations(dataPath, seed)

    so = SceneOptions("", saveFiles, tileSz, history, seed,-1,-1,-1,
                      ["" for i in 1:history+1], (-1,-1))

    distances = Vector{Float64}(undef, history)
    outputValues = Dict{Int32, Matrix{Float64}}()

    for (path, location) in zip(paths, locations)
	imgBands = bandSelection(path)
        updateLocation!(so, name, location, path)
        changeMap = zeros(so.changeMapSz)

        if saveFiles
	    metaData = Dict("experimentName"=>so.name,
                       "numTiles"=>so.numTiles,
                       "dataType"=>typeof(panns[imgBands[1]].nw.nodeVoltage[1]),
                       "dataShape"=>size(getReadoutValues(panns[imgBands[1]])),
                       "history"=>so.history,
                       "seed"=>so.seed,
                       "bands"=>imgBands,
                       "files"=>so.files)
            open(metaData["experimentName"]*"_metaData.json", "w") do f
                JSON3.write(f, metaData)
                println(f)
            end
        end

        p = Progress(so.numTiles, desc="Running experiment ", enabled=enableProgBar)
        i, j = 1, 1
        x, y = 1, 1
        while (i <= so.imgWidth - so.tileSz[1]+1) &&
              (j <= so.imgHeight - so.tileSz[2]+1)

            window = Int32.((j,i,so.tileSz[2],so.tileSz[1]))
            for b in imgBands
                dataLoaders[b].files .= so.files
                loadSignals!(dataLoaders[b], window)
                simulatePANN!(panns[b], dataLoaders[b].signals, T, dt)
                outputValues[b] = getReadoutValues(panns[b])
            end

            allReadouts = nothing
            for b in imgBands
                if b == imgBands[1]
                    allReadouts = outputValues[b]
                else
                    allReadouts = vcat(allReadouts, outputValues[b])
                end
            end

            for idx in 1:history
                distances[idx] = distanceMetric(allReadouts[:,idx], 
                                                allReadouts[:,end])
            end
            changeMap[x,y] = minimum(distances)


            if so.saveFiles
                fiName = name*"_"*location*"_output.dat"
		io = open(fiName, "a")
                write(io, allReadouts)
                close(io)
            end

            i += so.tileSz[1]
            y += 1
            if i > so.imgWidth - so.tileSz[1]+1
                i = 1
                y = 1
                j += so.tileSz[2]
                x += 1
            end
            next!(p)
        end

        replace!(changeMap, NaN=>0)
        mapName = name*"_"*location*"_changeMap.csv"
        open(mapName, "w") do f
            writedlm(f, changeMap, ",")
        end
    end
end
