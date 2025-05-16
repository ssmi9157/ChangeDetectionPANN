"""Contains the activity of the network """
struct Network
    "The connectivity information of the network."
    connectivity::Connectivity
    "The current across all electrodes, sources and drains, at each timestep."
    electrodeCurrent::Matrix{Float64}
    "The filament state (lambda) value for each edge, at each timestep."
    filamentState::Matrix{Float64}
    "The voltage across each edge, at each timestep."
    edgeVoltage::Matrix{Float64}
    "The conductance of each edge, at each timestep."
    edgeConductance::Matrix{Float64}
    "Dict of the values to calculate the electrical elements in the network."
    edgeParams::Dict{String, Real}
    "The voltage across each node, at each timestep."
    nodeVoltage::Matrix{Float64}
    "The overall conductance of the Network"
    conductance::Vector{Float64}
end
