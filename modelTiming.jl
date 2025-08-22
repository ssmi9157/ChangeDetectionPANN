"Struct for storing all the required variables to run a mock simulation"
struct setupVariables
    history::Int
    imgWidth::Int
    imgHeight::Int
    T::Float64
    dt::Float64
    signalsSize::Tuple{Int, Int}
    tileSz::Tuple{Int, Int}
    imgBands::Vector{Int}
    distances::Vector{Float64}
    outputValues::Dict{Int32, Matrix{Float64}}
    changeMap::Matrix{Float64}
    distanceMetric
    panns
end

"""
    setupV()

Returns a setupVariables struct that contains all the required variables to run
a mock simulation.
"""
function setupV()
    seed = 1
    history = 4
    imgWidth = 574
    imgHeight = 509
    tileSz = (32,32)
    connectivityFile = "803nw12279junctions256electrodes.mat"
    dt = 0.005
    T = (history+1) * dt
    bands = [2,3,4,5,6,7,8,9,12,13]
    inputNodes = getInputNodes(16, 256)
    numReadoutNodes = 400
    panns = createPannModels(bands, connectivityFile, inputNodes, 
                             numReadoutNodes, T, dt, seed)

    distanceMetric = CorrDist()
    distances = Vector{Float64}(undef, history)
    outputValues = Dict{Int32, Matrix{Float64}}()

    imgBands = [2]

    numSignals = length(1:2:tileSz[1])* length(1:2:tileSz[2])
    signalsSize = (history+1, numSignals)
    changeMapSz = (trunc(Int, imgHeight/tileSz[1]),
                      trunc(Int, imgWidth/tileSz[2]))
    changeMap = zeros(changeMapSz)

    v = setupVariables(history, imgWidth, imgHeight,
                       T, dt, signalsSize, tileSz,
                       imgBands, distances,
                       outputValues, changeMap,
                       distanceMetric, panns
                      )

    return v
end

"""
    processImage(v::setupVariables)

Runs a mock simulation as defined by the setup variables.

# Arguments
- 'v': Struct containing all the required variables to run a mock simulation.
"""
function processImage(v::setupVariables)
    i, j = 1, 1
    x, y = 1, 1
    while (i <= v.imgWidth - v.tileSz[1]+1) &&
          (j <= v.imgHeight - v.tileSz[2]+1)

        for b in v.imgBands
            simulatePANN!(v.panns[b], rand(v.signalsSize...), v.T, v.dt)
            v.outputValues[b] = getReadoutValues(v.panns[b])
        end

        allReadouts = nothing
        for b in v.imgBands
            if b == v.imgBands[1]
                allReadouts = v.outputValues[b]
            else
                allReadouts = vcat(allReadouts, v.outputValues[b])
            end
        end

        for idx in 1:v.history
            v.distances[idx] = v.distanceMetric(allReadouts[:,idx], 
                                              allReadouts[:,end])
        end
        v.changeMap[x,y] = minimum(v.distances)


        i += v.tileSz[1]
        y += 1
        if i > v.imgWidth - v.tileSz[1]+1
            i = 1
            y = 1
            j += v.tileSz[2]
            x += 1
        end
    end

    replace!(v.changeMap, NaN=>0)
end

"""
    timeModel()

Times the time to run a mock simulation of the PANN model.
"""
function timeModel()
    b = @benchmark processImage(v) setup = (v=setupV())
    display(b)
end
