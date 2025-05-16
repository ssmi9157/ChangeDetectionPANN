"""
    getFeaturePositions(expName::String, location::String)

Returns the output values for the requested natural disaster event in an 
experiment.

# Arguments
- 'expName': The name of the experiement.
- 'location': The location name for the natural disaster event.
"""
function getFeaturePositions(expName::String, location::String)
    metaDataFile = expName*"_"*location*"_metaData.json"
    metaData = JSON3.read(read(metaDataFile))

    dataShape = (metaData["dataShape"][1]*length(metaData["bands"]),
                 metaData["dataShape"][2], metaData["numTiles"])

    dataMap = Dict("Float64"=>Float64, "Float32"=>Float32,
                "Int64"=>Int64, "Int32"=>Int32)
    dataType = dataMap[metaData["dataType"]]
    data = Array{dataType}(undef, dataShape)

    dataFile = expName*"_"*location*"_output.dat"
    io = open(dataFile, "r")
    read!(io, data)
    close(io)

    return data
end

"""
    umapTransform(values, dims; neighbours::Int=20, metric=CorrDist(), 
                  minDist::Float64=0.4)

Returns the UMAP dimensionally reduced data of the given 'values'

# Arguments
- 'values': The data values to be reduced.
- 'dims': The number of dimensions to reduce the data to.
- 'neighbours': Hyperparamter for the UMAP reduction.
- 'metric': The distance metric to use in the UMAP process.
- 'minDist': Hyperparameter for the minimum distance in the UMAP reduction.
"""
function umapTransform(values, dims; neighbours::Int=20, metric=CorrDist(), 
                       minDist::Float64=0.4)
    sz = size(values)
    values = reshape(values, (sz[1], :))
    m = UMAP_(values, dims, n_neighbors=neighbours, metric=metric, 
              min_dist=minDist)
    return reshape(m.embedding, (dims, sz[2],sz[3]))
end

"""
    getTileImage!(buffer, handle, a::Int64, b::Int64; zoom::Float64=0.25)

Returns the RGB image of the given tile as an OffSetImage for plotting.

# Arguments
- 'buffer': The pre-initialised matrix the data is to be loaded into.
- 'handle': The ArchGDAL handle of the .tif file.
- 'a': The index of the first column in the window to be loaded.
- 'b': The index of the first row in the window to be loaded.
- 'zoom': Controls how the large the tile images will appear when plotted.
"""
function getTileImage!(buffer, handle, a::Int64, b::Int64; zoom::Float64=0.25)
    window = (b,a,32,32)
    AG.read!(handle, buffer, Int32.([4,3,2]),
                   window[1]:window[1]+window[3]-1,
                   window[2]:window[2]+window[4]-1)
    tile = permutedims(buffer, [2,1,3])

    vals = ((490,1810), (560,1635), (990,1810))

    # Clip values for better visualisation
    for i in 1:3
        band = tile[:,:,i]
        band[band .< vals[i][1]] .= vals[i][1]
        band[band .> vals[i][2]] .= vals[i][2]
        tile[:,:,i] .= band
    end

    for i in 1:3
        tile[:,:,i] .= ((tile[:,:,i] .- vals[i][1]) ./
                        (vals[i][2] - vals[i][1]))
    end

    return matplotlib.offsetbox.OffsetImage(tile, zoom=zoom)
end

"""
    tileChange(target, numTiles::Int)

Returns the change percentage for each tile in the target image.

# Arguments
- 'target': The target image that contains the which pixels are changes.
- 'numTiles': The number of tiles in the target image.
"""
function tileChange(target, numTiles::Int)
    tileSz = 32

    changePercent = Vector{Float64}(undef, numTiles)
    tile = 1
    for x in 1:tileSz:size(target)[1], y in 1:tileSz:size(target)[2]
        changePercent[tile] = sum(target[x:(x+tileSz-1), y:(y+tileSz-1)]) / 
	                           (tileSz^2)
        tile += 1
    end

    return changePercent .* 100
end

"""
    axTransform!(ax, fig, transform, path; legend::Bool=true,
            axisOff::Bool=true, cbarSize::Float64=0.5, cbarLoc::String="left",
	    legendLoc::String="lower left")

Adds the scatter points to the axis corresponding to the position of the tiles 
in the dimensionally reduced feature space.

# Arguments
- 'ax': The axis to be adding the sctter points to.
- 'fig': The figure the axis belongs to.
- 'transform': The UMAP dimensionally reduced values of the output values.
- 'path': The path to the target change map.
- 'legend': Determines if the legend is added to the axis.
- 'axisOff': Determines if the markers on the x- and y-axis are included.
- 'cbarSize': The relative size of the colour bar.
- 'cbarLoc': The location of the colour bar relative to the axis.
- 'legendLoc': The location of the legend on the axis.
"""
function axTransform!(ax, fig, transform, path; legend::Bool=true,
            axisOff::Bool=true, cbarSize::Float64=0.5, cbarLoc::String="left",
	    legendLoc::String="lower left")

    tileSz = 32
    file = readdir(path)[1]
    target = loadChangeMask(joinpath(path,file))
    targetSz = size(target)
    outputImgSize = [trunc(Int,targetSz[1]/tileSz),
                     trunc(Int,targetSz[2]/tileSz)]
    numTiles = outputImgSize[1]*outputImgSize[2]

    outputImgSize .*= tileSz
    target = target[1:outputImgSize[1], 1:outputImgSize[2]]
    target[target .== 2] .= 0
    
    changePercent = tileChange(target, numTiles)

    alpha = 0.7
    colours = ("tab:blue", "tab:green", "tab:purple", "tab:pink")
    for i in 1:(size(transform)[2]-1)
        time = size(transform)[2] - i
	ax.scatter(transform[1,i,:], transform[2,i,:], label="Before (t-$time)",
                   c=colours[time], alpha=alpha)
    end
    x = get_cmap("Oranges")
    startSpot = 0.4
    len = trunc(Int, x.N * startSpot)
    cvalues = x(range(startSpot, 1, length=len))
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("Upper Half", 
                                                                cvalues)
    a = ones(length(changePercent)) .* alpha
    map = ax.scatter(transform[1,end,:], transform[2,end,:], label="After", 
                     c=changePercent, cmap=cmap2, alpha=a, vmin=0, vmax=100)
    fig.colorbar(map, ax=ax, label="Tile Change Percentage", shrink=cbarSize, 
                 location=cbarLoc)
    ax.set_aspect(1)
    if axisOff
        ax.set_axis_off()
    end
    if legend
        ax.legend(loc=legendLoc)
    end
    return ax
end

"""
    axTransformImages!(ax, transform, path; axisOff::Bool=true)

Adds the RGB tile images to the axis corresponding to the position of the tiles
in the dimensionally reduced feature space.

# Arguments
- 'ax': The axis to be adding the sctter points to.
- 'transform': The UMAP dimensionally reduced values of the output values.
- 'path': The path to the directory containing the .tif files.
- 'axisOff': Determines if the markers on the x- and y-axis are included.
"""
function axTransformImages!(ax, transform, path; axisOff::Bool=true)

    files = [joinpath(path, fi) for fi in readdir(path)]
    files = files[(end-(size(transform)[2]-1)):end]

    tileSz = 32
    handle = AG.read(files[1])
    h = AG.height(handle)
    w = AG.width(handle)
    numTiles = trunc(Int,h/tileSz) * trunc(Int,w/tileSz)

    buffer = Array{Float32}(undef, (tileSz,tileSz,3))
    for i in 1:size(transform)[2]
        handle = AG.read(files[i])
        a, b = 1, 1

        ax.scatter(transform[1,i,:], transform[2,i,:])

        for j in 1:numTiles
            x = transform[1,i,j]
            y = transform[2,i,j]

            ab = matplotlib.offsetbox.AnnotationBbox(
                     getTileImage!(buffer, handle, a, b), (x,y), frameon=false)
            ax.add_artist(ab)

            a += tileSz
            if a > w - tileSz +1
                a = 1
                b += tileSz
            end
        end
    end
    ax.set_aspect(1)
    if axisOff
        ax.set_axis_off()
    end
    return ax
end

"""
    getPixelPositions(path::String)

Returns the location of the tiles in the pixel space.

# Arguments
- 'path': The path to the directory containing the .tif files.
"""
function getPixelPositions(path::String)
    files = [joinpath(path, fi) for fi in readdir(path)]

    tileSz = 32
    handle = AG.read(files[1])
    h = AG.height(handle)
    w = AG.width(handle)
    imgSz = [trunc(Int,h/tileSz), trunc(Int,w/tileSz)]
    numTiles = imgSz[1] * imgSz[2]
    imgSz .*= tileSz

    buffer = Array{Float32}(undef, (imgSz[2],imgSz[1],10))
    values = Array{Float32}(undef, (tileSz*tileSz*10, length(files), numTiles))
    bands = [2,3,4,5,6,7,8,9,12,13]

    for (f,fi) in enumerate(files)
        handle = AG.read(fi)
        dat = AG.read!(handle, buffer, bands, 1:imgSz[1], 1:imgSz[2])
        dat = permutedims(dat, [2,1,3])
        replace!(dat, NaN=>0)

        i = 1
        for y in 1:tileSz:imgSz[2], x in 1:32:imgSz[1]
            tile =  dat[x:x+tileSz-1, y:y+tileSz-1, :]

            values[:,f,i] .= reshape(tile, :)
            i += 1
        end
    end
    return values
end
