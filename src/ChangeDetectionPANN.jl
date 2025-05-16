module ChangeDetectionPANN

using ArchGDAL; const AG = ArchGDAL
using DelimitedFiles
using Distances
using Graphs
using ImageTransformations
using Interpolations
using JSON3
using LinearAlgebra
using MAT
using ProgressMeter
using PyCall
using PyPlot
using Random
using Statistics
using UMAP

const SK = PyNULL()

function __init__()
    copy!(SK, pyimport_conda("sklearn.metrics", "scikit-learn"))
end

export 
    runScenes,
    evaluateChangeMaps,
    evaluateAllChangeMaps,
    evaluateErrors,
    evaluateAllErrors,
    getFeaturePositions,
    umapTransform,
    getPixelPositions,
    axTransform!,
    axTransformImages!

include("core/simSettings.jl")
include("core/networkStruct.jl")
include("core/getPairings.jl")
include("core/simulate.jl")

include("changeDetection/auprc.jl")
include("changeDetection/dataLoader.jl")
include("changeDetection/bandSelector.jl")
include("changeDetection/sceneOptions.jl")
include("changeDetection/pannModel.jl")
include("changeDetection/tileScene.jl")
include("changeDetection/evaluateChangeMap.jl")
include("changeDetection/featureSpacePlots.jl")

end # module
