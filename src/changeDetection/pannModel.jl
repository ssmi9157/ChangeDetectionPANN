"Struct for containing all the information for running an individual PANN 
network"
struct PannModel
    "Contains all the activity of the network."
    nw::Network
    "Contains the current state of each edge in the network."
    es::EdgeState
    "The left-hand side of the main matrix solved in the network."
    lhs::Matrix{Float64}
    "The right-hand side of the main matrix solved in the network."
    rhs::Vector{Float64}
    "List of the input nodes of the network."
    inputNodes::Vector{Int32}
    "List of the readout nodes of the network."
    readoutNodes::Vector{Int32}
end


"""
    createEdgeParams(criticalFlux=1e-2, maxFlux=1.5e-2, boost=1,
                        Ron=1.287e4, Roff=1.287e7, Vset=1e-2, Vreset=5e-3)

Returns a dictionary of the hyperparameters for the network.
"""
function createEdgeParams(criticalFlux=1e-2, maxFlux=1.5e-2, boost=1,
                        Ron=1.287e4, Roff=1.287e7, Vset=1e-2, Vreset=5e-3)

    edgeParams = Dict("criticalFlux" => criticalFlux,
                      "maxFlux" => maxFlux,
                      "Ron" => Ron, "Roff" => Roff,
                      "Vset" => Vset, "Vreset" => Vreset,
                      "boost" => boost)
    return edgeParams
end

"""
    createPannModel(connectivityFile::String, inputNodes::Vector{Int},
                        numReadoutNodes::Int, T::Real, dt::Real, 
                        edgeParams::Dict{String, Real}; seed::Int=1)

Returns the PannModel Struct with all the specified options.

# Arguments
- 'connectivityFile': Name of the connectivity file that contains the network 
  structure information.
- 'inputNodes': List of nodes that will be the inputs.
- 'numReadoutNodes': Number of nodes to be used as readouts.
- 'T': Simulated time for each tile.
- 'dt': Simulated sample rate.
- 'edgeParams': Hyperparameters for the network.
- 'seed': Seed used to randomly select the readout nodes.
"""
function createPannModel(connectivityFile::String, inputNodes::Vector{Int},
                        numReadoutNodes::Int, T::Real, dt::Real, 
                        edgeParams::Dict{String, Real}; seed::Int=1)

    connectivity = createConnectivity(connectivityFile)
    N = size(connectivity.adjMat)[1]
    E = connectivity.numOfEdges
    timeVector = dt:dt:T
    niterations = length(timeVector)
    numOfInputNodes = length(inputNodes)

    Random.seed!(seed)
    readoutNodes = shuffle(1:connectivity.numOfNodes)[1:numReadoutNodes]

    nw = Network(connectivity, zeros(numOfInputNodes, niterations),
                      zeros(E, niterations), zeros(E, niterations),
                      zeros(E, niterations), edgeParams, zeros(N, niterations),
                      zeros(niterations))

    es = createEdgeState(nw.edgeVoltage[:,end], nw.edgeConductance[:,end],
                         nw.filamentState[:,end])

    lhs = zeros((N+numOfInputNodes, N+numOfInputNodes))
    rhs = zeros(N+numOfInputNodes)

    return PannModel(nw, es, lhs, rhs, inputNodes, readoutNodes)
end

"""
    simulatePANN!(pann::PannModel, signals::Matrix, T::Real, dt::Real)

Updates the values of the PannModel for simulating the input 'signals'.

# Arguments
- 'pann': PannModel used to simulate the input signals being fed into.
- 'signals': Input value signal streams to be fed into the network.
- 'T': Simulated time for the signals.
- 'dt': Simulated sample rate.
"""
function simulatePANN!(pann::PannModel, signals::Matrix, T::Real, dt::Real)

    stimulus = [createStimulus(biasType="Custom", T=T, dt=dt, 
           customSignal=signals[:,i]) for i in 1:size(signals)[2]]

    stimulusOptions = createSimulationOptions(stimulus, pann.inputNodes,
                                            dt:dt:T, dt=dt, T=T,
                                            connectivity=pann.nw.connectivity,
                                            contactMode="preSet")

    simulateNetworkBatch!(pann.nw, pann.es, pann.lhs, pann.rhs, 
                          stimulusOptions, enableProgBar=false)

    pann.es.voltage .= pann.nw.edgeVoltage[:,end]
    pann.es.conductance .= pann.nw.edgeConductance[:,end]
    pann.es.filamentState .= pann.nw.filamentState[:,end]
end

"""
    getReadoutValues(pann::PannModel)

Returns the output values from the readout nodes.
"""
function getReadoutValues(pann::PannModel)
    return pann.nw.nodeVoltage[pann.readoutNodes,:]
end
