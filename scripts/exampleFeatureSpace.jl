using Pkg
Pkg.activate("/path/to/ChangeDetectionPANN")

using ChangeDetectionPANN
using Distances
using PyPlot

function exampleFeatureSpacePlot(expName::String,location::String,path::String; 
                       mode::String="umap", distMetric=CorrDist(), 
                       neighbours::Int64=20, minDist::Float64=0.4, dims::Int=2)

    inputImgPath = joinpath(path, location, "S2")
    targetPath = joinpath(path, location, "changes")

    values = getFeaturePositions(expName, location)

    transform = umapTransform(values, dims, neighbours=neighbours,
                metric=distMetric, minDist=minDist)


    fig, axs = subplots(1,2, dpi=300, figsize=(14,7), 
                        gridspec_kw=Dict("width_ratios"=>[1,0.75]))

    # Only use the points corresponding to directly before and
    # after the natural disaster event for plotting
    transform = transform[:,4:5,:]

    axTransform!(axs[1], fig, transform, targetPath)
    axTransformImages!(axs[2], transform, inputImgPath)
    tight_layout()
    fig.savefig("exampleFeatureSpace.png")
end

name = "testExp"
location = "landslide_bigsur"
path = "/path/to/data/landslides"
exampleFeatureSpacePlot(name, location, path)
