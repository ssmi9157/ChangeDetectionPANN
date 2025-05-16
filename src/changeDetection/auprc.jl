"""
    getClassVectors(output, target, tileSz::Tuple{Int, Int})

Returns the change score of non-cloudy values of the outputs as a 
vector 'x' and the corresponding target vector 'y'.

# Arguments
- 'output': Matrix with all the output change scores.
- 'target': Matrix with the target values.
- 'tileSz': Tuple that contains the width and height of each tile.
"""
function getClassVectors(output, target, tileSz::Tuple{Int, Int})
    if size(output) == size(target)
        changesMap = output
    else
        changesSize = size(output).*Tuple(tileSz)
        changesMap = imresize(output, changesSize, method=BSpline(Constant()))
    end
    x = vec(changesMap)
    y = vec(target)
    keep = y.!= 2
    x = x[keep]
    y = y[keep]

    return x, y
end

"""
    calculateAUPRC(output, target, tileSz::Tuple{Int, Int})

Returns the auprc value for the given output and target values.

# Arguments
- 'output': Matrix with all the output change scores.
- 'target': Matrix with the target values.
- 'tileSz': Tuple that contains the width and height of each tile.
"""
function calculateAUPRC(output, target, tileSz::Tuple{Int, Int})
    vecs = getClassVectors(output, target, tileSz)

    precision, recall, thresholds = SK.precision_recall_curve(vecs[2],vecs[1])
    auprc = SK.auc(recall, precision)

    return auprc
end

"""
    calculateTotalAUPRC(outputs, targets, tileSz::Tuple{Int, Int})

Returns the auprc value for each of the outputs and targets.

# Arguments
- 'output': List of Matrices with all the output change scores.
- 'target': List of Matrices with the target values.
- 'tileSz': Tuple that contains the width and height of each tile.
"""
function calculateTotalAUPRC(outputs, targets, tileSz::Tuple{Int, Int})
    x = []
    y = []
    for (output, target) in zip(outputs, targets)
        vecs = getClassVectors(output, target, tileSz)

        append!(x, vecs[1])
        append!(y, vecs[2])
    end

    precision, recall, thresholds = SK.precision_recall_curve(y,x)
    auprc = SK.auc(recall, precision)
    return auprc
end

