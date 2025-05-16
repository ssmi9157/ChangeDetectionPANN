"""
    getFarthestPairing(connectivity::Connectivity)

Returns the index of the nodes which are furthest apart in the network.

# Arguments:
- 'connectivity': The connectivity information of the network.
"""
function getFarthestPairing(connectivity::Connectivity)
    G = Graph(connectivity.adjMat)
    distMat = floyd_warshall_shortest_paths(G).dists
    f = findall(distMat .== maximum(distMat))

    # Convert from CartesianIndex to Matrix
    farthest = [[i[2] for i in f] [i[1] for i in f]]

    if length(farthest[:,1]) > 2
        xc0 = connectivity.xc[farthest[:,1]]
        yc0 = connectivity.yc[farthest[:,1]]
        xc1 = connectivity.xc[farthest[:,2]]
        yc1 = connectivity.yc[farthest[:,2]]
        idx = argmax((xc1.-xc0).^2 .+ (yc1.-yc0).^2)
    else
        idx = 1
    end
    return farthest[idx,:]
end


"""
    getBoundaryPairing(connectivity::Connectivity, [numOfPairs::Int=5])

Returns the set of boundary electrode pairs, for the network. 
Pairs are ordered as [a1, b1, c1, ..., c2, b2, a2].

# Arguments:
- 'connectivity': The connectivity information of the network.
- 'numOfPairs': The number of boundary pairs to be returned.
"""
function getBoundaryPairing(connectivity::Connectivity; numOfPairs::Int=5)
    centerX = connectivity.xc
    electrodes = zeros(Int, 2*numOfPairs)
    electrodes[1:numOfPairs] = sortperm(centerX,dims=2)[1:numOfPairs]
    electrodes[numOfPairs+1:end]=sortperm(centerX,dims=2)[end-numOfPairs+1:end]
    return electrodes
end
