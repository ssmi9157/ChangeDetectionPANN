"""
    simulateNetworkBatch!(network::Network, edgeState::EdgeState,
                          lhs, rhs, simulationOptions::SimulationOptions;
                          enableProgBar=true)

Returns the updated network.

# Arguments:
- 'network': A pre-initialised network for the simulation.
- 'edgeState': The current electrical state of each edge in the network.
- 'lhs': The lhs of the main matrix being solved.
- 'rhs': The rhs of the main matrix being solved.
- 'simulationOptions': Contains general simulation details.
"""
function simulateNetworkBatch!(network::Network, edgeState::EdgeState, 
                          lhs, rhs, simulationOptions::SimulationOptions;
                          enableProgBar=true)

    niterations = simulationOptions.numOfIterations
    electrodes = simulationOptions.electrodes
    dt = simulationOptions.dt
    numOfElectrodes = length(electrodes)
    E = network.connectivity.numOfEdges
    V = size(network.connectivity.adjMat)[1]
    edgeList = network.connectivity.edgeList

    p = Progress(niterations, desc="Running Simulation ", enabled=enableProgBar)
    for thisTime in 1:niterations

        updateConductance!(edgeState)
        edgeConductance = edgeState.conductance

        fill!(lhs, 0.0)
        for i = 1:E
            lhs[edgeList[i,1], edgeList[i,2]] = -edgeConductance[i]
            lhs[edgeList[i,2], edgeList[i,1]] = -edgeConductance[i]
        end
        for i in 1:V
            lhs[i,i] = -sum(lhs[:,i])
        end

        fill!(rhs, 0.0)

        for (i, thisElec) in enumerate(electrodes)
            lhs[V+i, thisElec] = 1
            lhs[thisElec, V+i] = 1
            rhs[V+i] = simulationOptions.stimulus[i].signal[thisTime]
        end

	sol = lhs\rhs

        edgeState.voltage .= sol[edgeList[:,1]] - sol[edgeList[:,2]]

        updateEdgeState!(edgeState, dt)

        network.nodeVoltage[:,thisTime] .= sol[1:V]
        network.electrodeCurrent[:,thisTime] .= sol[V+1:end]
        network.edgeVoltage[:,thisTime] .= edgeState.voltage
        network.edgeConductance[:,thisTime] .= edgeState.conductance
        network.filamentState[:,thisTime] .= edgeState.filamentState

       next!(p)
    end

    totalVoltage = zeros(niterations)
    for stim in simulationOptions.stimulus
        totalVoltage .+= stim.signal
    end
    network.conductance .= network.electrodeCurrent[end,:] ./totalVoltage

    return network
end
