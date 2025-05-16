" The connectivity information of the network."
struct Connectivity
    "Adjacency matrix of the network."
    adjMat::Matrix{Int}
    "A EX2 matrix of vertex indices, where each column represents an edge. 
    This field follows two conventions:
    1. edgeList(:,1) < edgeList(:,2)
    2. The index of an edge is defined as the index of the column containing 
        the indices of the two vertices it connects."
    edgeList::Matrix{Int}
    "Number of nodes in the network."
    numOfNodes::Int
    "Number of edges in the network."
    numOfEdges::Int
    "Number of electrodes in the network."
    numOfElectrodes::Int
    "x coordinates of edge positions."
    xi::Matrix{Float64}
    "x coordinates of one end of the edges."
    xa::Matrix{Float64}
    "x coordinates of the other end of the edges."
    xb::Matrix{Float64}
    "x coordinates of the nodes."
    xc::Matrix{Float64}
    "y coordinates of edge positions."
    yi::Matrix{Float64}
    "y coordinates of one end of the edges."
    ya::Matrix{Float64}
    "y coordinates of the other end of the edges."
    yb::Matrix{Float64}
    "y coordinates of the nodes."
    yc::Matrix{Float64}
    "x coordinates of the electrodes."
    xe::Matrix{Float64}
    "y coordinates of the electrodes."
    ye::Matrix{Float64}
    "Diameter of the electrodes (μm)."
    elecDiameter::Float64
    "x length of the area the network is located on."
    lengthX::Int
    "y length of the area the network is located on."
    lengthY::Int
end

"""
    createConnectivityElectrodes(filename::String)

Returns the connectivity struct from the given file name.

# Arguments
- 'filename': The name of a pregenerated .mat file with all the fields for 
  the connectivity struct.
"""
function createConnectivity(filename::String) 

    try
        global matfile = matread(filename)
    catch
        rootPath = dirname(@__FILE__)
        rootPath = normpath(rootPath, "..")
        fullPath = joinpath(rootPath, "connectivityData", filename)
        global matfile = matread(fullPath)
    end

    con = Connectivity(matfile["adj_matrix"],
                       matfile["edge_list"],
                       matfile["number_of_wires"],
                       matfile["number_of_junctions"],
                       matfile["number_of_electrodes"],
                       matfile["xi"],
                       matfile["xa"],
                       matfile["xb"],
                       matfile["xc"],
                       matfile["yi"],
                       matfile["ya"],
                       matfile["yb"],
                       matfile["yc"],
                       matfile["xe"],
                       matfile["ye"],
                       matfile["electrode_diameter"],
                       matfile["length_x"],
                       matfile["length_y"])

    # Correct for index starting at 1 in julia not 0 like python
    con.edgeList .+= 1

    return con
end

" Contains the current electrical state of each edge in the network."
struct EdgeState
    "The voltages across each edge."
    voltage::Vector{Float64}
    "The conductance of each edge."
    conductance::Vector{Float64}
    "The on conductance value of the edges."
    onConductance::Float64
    "The off conductance value of the edges."
    offConductance::Float64
    "The filament state (lambda) of each edge."
    filamentState::Vector{Float64}
    "Value of the set voltage."
    setVoltage::Float64
    "Value of the reset voltage."
    resetVoltage::Float64
    "Value of the critical flux."
    criticalFlux::Float64
    "Value of the maximum flux."
    maxFlux::Float64
    "The memory decay rate of the filaments."
    boost::Float64
end


"""
    createEdgeState(numOfEdges::Int; 
                [setVoltage::Float64=1e-2, resetVoltage::Float64=5e-3,
                onResistance::Float64=1.287e4, offResistance::Float64=1.287e7,
                criticalFlux::Float64=1e-2, maxFlux::Float64=1.5e-2,
                boost::Float64=1.0])

Return an initial edge state, where each edge is in the off state.

# Arguments:
- 'numOfEdges': The number of edges in the network.
"""
function createEdgeState(numOfEdges::Int;
                setVoltage::Float64=1e-2, resetVoltage::Float64=5e-3,
                onResistance::Float64=1.287e4, offResistance::Float64=1.287e7,
                criticalFlux::Float64=1e-2, maxFlux::Float64=1.5e-2,
                boost::Float64=1.0)

    voltage = zeros(numOfEdges)
    conductance = zeros(numOfEdges)
    filamentState = zeros(numOfEdges)

    onConductance = 1/onResistance
    offConductance = 1/offResistance

    return EdgeState(voltage, conductance, onConductance, 
                     offConductance, filamentState, setVoltage, 
                     resetVoltage, criticalFlux, maxFlux, boost)
end

"""
    createEdgeState(voltage::Vector{Float64},
                conductance::Vector{Float64}, filamentState::Vector{Float64};
                [setVoltage::Float64=1e-2, resetVoltage::Float64=5e-3,
                onResistance::Float64=1.287e4, offResistance::Float64=1.287e7,
                criticalFlux::Float64=1e-2, maxFlux::Float64=1.5e-2,
                boost::Float64=1.0])

Return an edge state with the given voltage, conductance and filamentState.

# Arguments:
- 'voltage': The voltage across each edge in the network.
- 'conductance': The conductance for each edge in the network.
- 'filamentState': The filamentState (lambda) for each edge in the network.
"""
function createEdgeState(voltage::Vector{Float64}, 
		conductance::Vector{Float64}, filamentState::Vector{Float64};
                setVoltage::Float64=1e-2, resetVoltage::Float64=5e-3,
                onResistance::Float64=1.287e4, offResistance::Float64=1.287e7,
                criticalFlux::Float64=1e-2, maxFlux::Float64=1.5e-2,
                boost::Float64=1.0)

    onConductance = 1/onResistance
    offConductance = 1/offResistance

    return EdgeState(voltage, conductance, onConductance,
                     offConductance, filamentState, setVoltage,
                     resetVoltage, criticalFlux, maxFlux, boost)
end

"""
    updateConductance!(es::EdgeState)

Updates the conductance value of each edge for the next timestep."

# Arguments:
- 'es': The edge state being updated.
"""
function updateConductance!(es::EdgeState)
    phi = 0.81
    C0 = 10.19
    J1 = 0.0000471307
    A = 0.17
    d = (es.criticalFlux .- abs.(es.filamentState)).*5 ./es.criticalFlux
    d[d.<0] .= 0
    tun = 1/A .* d ./ phi^0.5 .* exp.(C0*(phi^0.5) .* d)./J1
    es.conductance .= 1 ./(tun .+ 1/es.onConductance) .+ es.offConductance
end


"""
    updateEdgeState!(es::EdgeState, dt::Float64)

Updates the filament state of each edge for the next timestep."

# Arguments:
- 'es': The edge state being updated.
- 'dt': The amount of time between timesteps (s).
"""
function updateEdgeState!(es::EdgeState, dt::Real)
    es.filamentState .= es.filamentState .+ 
                            (abs.(es.voltage) .> es.setVoltage) .*
                            (abs.(es.voltage) .- es.setVoltage) .*
                            sign.(es.voltage) .* dt

    es.filamentState .= es.filamentState .-
                            (es.resetVoltage .> abs.(es.voltage)) .*
                            (es.resetVoltage .- abs.(es.voltage)) .*
                            sign.(es.filamentState) .* dt .* es.boost
    es.filamentState[abs.(es.filamentState) .> es.maxFlux] =
            sign.(es.filamentState[abs.(es.filamentState) .> es.maxFlux]) .* es.maxFlux

end

" Contains details of the external stimulus."
struct Stimulus
    "Avaliable types of signals:
    Drain | DC | AC | Square | Triangular | Sawtooth | Pulse | MKG | Custom"
    biasType::String
    "Time betweeen timesteps (s)."
    dt::Real
    "Duration of simulation (s)."
    T::Real
    "The start time for the signal."
    onTime::Float64
    "The end time for the signal."
    offTime::Float64
    "Amplitude of the signal when on (V)."
    onAmp::Float64
    "Amplitude of the signal when off (V). Normally non-zero."
    offAmp::Float64
    "Value of the frequency."
    f::Float64
    "The phase shift of the signal."
    shift::Float64
    "The input stimulus signal."
    signal::Vector{Float64}
end


"""
    createStimulus(; biasType::String="DC", customSignal=nothing, 
                dt::Float64 = 1e-3, T::Float64=10.0,
                onTime::Float64=0.0, offTime::Float64=5e7,
                onAmp::Float64=1.0, offAmp::Float64=5e-3,
                f::Float64=1.0, shift::Float64=0.0)

Return a stimulus struct that contains the details of the external stimulus.
"""
function createStimulus(; biasType::String="DC", customSignal=nothing, 
                dt::Real = 1e-3, T::Real=10,
                onTime::Float64=0.0, offTime::Float64=5e7,
                onAmp::Float64=1.0, offAmp::Float64=5e-3,
                f::Float64=1.0, shift::Float64=0.0)

    period = 1/f
    timeVector = dt:dt:T
    offIndex = findall((timeVector .< onTime) .& (timeVector .>= offTime))

    if biasType == "Drain"
        signal = zeros(length(timeVector))
    elseif biasType == "DC"
        onIndex = findall((timeVector .>= onTime) .& (timeVector .<= offTime))
        signal = ones(length(timeVector)) .* offAmp
        signal[onIndex] .= onAmp
    elseif biasType == "AC"
        signal = onAmp .* sin.(2 .* pi .*f .* timeVector)
    elseif biasType == "Custom"
        if length(customSignal) >= length(timeVector)
            signal = customSignal
        else
	    len = length(timeVector)
            println("Signal length is shorter then simulation length. Current timeVector length is $len !")
            exit()
        end
    else
        println("Stimulus type error. Options are: ")
        println("| Drain | DC | AC | Custom |")
        exit()
    end

    if shift != 0
        signal .+= shift
    end

    if (biasType != "Custom") && (biasType != "Drain")
        signal[offIndex] .= offAmp
    end

    return Stimulus(biasType, dt, T, onTime, offTime, onAmp, offAmp, f, 
                    shift, signal)
end


" Contains general simulation details."
struct SimulationOptions
    "Time betweeen timesteps (s)."
    dt::Real
    "Duration of simulation (s)."
    T::Real
    "Vector with the time of each step in the simulation."
    timeVector
    "The number of timesteps in the simulation."
    numOfIterations::Int
    "The contact mode of the electrodes. Valid modes are:
    preSet, farthest, boundary."
    contactMode::String
    "The indices of all the electrode nodes."
    electrodes::Vector{Int}
    "The indices of only the drain electrode nodes."
    drains::Vector{Int}
    "The indices of only the source electrode nodes."
    sources::Vector{Int}
    "The external stimulus being applied to the electrodes."
    stimulus::Vector{Stimulus}
end


"""
    createSimulationOptions(stimulus::Vector{Stimulus},
                electrodes, timeVector;
                dt::Float64=1e-3, T::Float64=1e1,
                connectivity::Connectivity=nothing,
                contactMode::String="farthest")

Returns a SimulationOptions that contains general simulation details.

# Arguments:
- 'stimulus': Vector of the external stimulus that is applied to the electrodes.
- 'electrodes': The indicies of the nodes that are the electrodes.
- 'timeVector': Vector with the time of each step in the simulation.
"""
function createSimulationOptions(stimulus::Vector{Stimulus},
                electrodes, timeVector;
                dt::Real=1e-3, T::Real=1,
                connectivity::Connectivity=nothing,
                contactMode::String="farthest")

    numOfIterations = length(timeVector)

    if contactMode == "preSet" && electrodes == nothing
        println("Error: contactMode is preSet but no electrodes given.")
        exit()
    elseif contactMode == "farthest"
        electrodes = getFarthestPairing(connectivity)
    elseif contactMode == "boundary"
        electrodes = getBoundaryPairing(connectivity)
    end

    drains = []
    sources = []
    for i in 1:length(electrodes)
        if mean(stimulus[i].signal) != 0
            append!(sources, electrodes[i])
        else
            append!(drains, electrodes[i])
        end
    end

    if length(drains) == 0
        append!(drains, electrodes[1])
    end

    return SimulationOptions(dt, T, timeVector, numOfIterations, contactMode,
                             electrodes, drains, sources, stimulus)
end
