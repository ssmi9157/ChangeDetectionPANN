"Struct for loading in the data from the .tif files and creating the signals
to fed into the model"
struct DataLoader
    "The minimum value for normalising the data."
    minVal::Float32
    "The maximum value for normalising the data."
    maxVal::Float32
    "The band channel to be loaded."
    band::Int32
    "The location to store the input data when first read in."
    buffer::Array{Float32}
    "The location to store the signal data for the network."
    signals::Matrix{Float32}
    "The list of files that make up the scene."
    files::Vector{String}
    "The size of the pooling kernel. Assumed to be square."
    poolSz::Int32
    "The stride size for the max pooling"
    stride::Int32
end

"""
    createDL(minVal::Real, maxVal::Real, band::Int, tileSz::Tuple{Int,Int},
             history::Int, poolSz::Int, stride::Int)

Returns a new DataLoader struct.

# Arguments
- 'minVal': The minimum value used in the normalisation of the data.
- 'maxVal': The maximum value used in the normalisation of the data.
- 'band': The spectral band of the .tif files to be loaded.
- 'TileSz': The width and height of each tile.
- 'history': The number of .tif images before the natural disaster to be loaded.
- 'poolSz': The size of the square pooling kernel, for the signal creation.
- 'stride': The stride length of the pooling kernel.
"""
function createDL(minVal::Real, maxVal::Real, band::Int, tileSz::Tuple{Int,Int},
                  history::Int, poolSz::Int, stride::Int)

    numSignals = length(1:stride:tileSz[1])* length(1:stride:tileSz[2])

    buffer = Array{Float32}(undef, (tileSz[1],tileSz[2],1))
    signals = Matrix{Float32}(undef,(history+1,numSignals))

    return DataLoader(minVal, maxVal, band, buffer, signals, 
                      ["" for f in 1:(history+1)], poolSz, stride)
end

"""
    processData!(dat, dl::DataLoader)

Returns the normalised data.

# Arguments
- 'dat': Matrix of the data to be normalised
- 'dl': The DataLoader object used to load the data.
"""
function processData!(dat, dl::DataLoader)
    # These were found manually
    vals = ((7.3,7.6),  #B1
            (6.9,7.5),  #B2
            (6.5,7.4),  #B3
            (6.2,7.5),  #B4
            (6.1,7.5),  #B5
            (6.5,8),    #B6
            (6.5,8),    #B7
            (6.5,8),    #B8
            (6.5,8),    #B8A
            (6,7),      #B9
            (2.5,4.5),  #B10
            (6,8),      #B11
            (6,8),      #B12
           )

    dat[:,:,1] .= ((log.(dat[:,:,1]) .- vals[dl.band][1]) ./ 
                       (vals[dl.band][2] - vals[dl.band][1]))
                       .*(dl.maxVal-dl.minVal) .+ dl.minVal

    return dat
end

"""
    loadTIF!(dl::DataLoader, path::String, window::NTuple{4,Int32})

Returns the data from the .tif file given by path.

# Arguments
- 'dl': The DataLoader object used to load the data.
- 'path': The path to the .tif file to be loaded.
- 'window': The area of the window to be loaded. The first to points define the 
  left corner of the window and the last two points define how far beyond to 
  load in the respective direction.
"""
function loadTIF!(dl::DataLoader, path::String, window::NTuple{4,Int32})
    handle = AG.read(path)
    dat = AG.read!(handle, dl.buffer, [dl.band],
                   window[1]:window[1]+window[3]-1,
                   window[2]:window[2]+window[4]-1)
    return permutedims(dat, [2,1,3])
end

"""
    loadSignals!(dl::DataLoader, window::NTuple{4,Int32})

Converts the data in the images of the window location into the signals to be 
fed into the model.

# Arguments
- 'dl': The DataLoader object used to load the data.
- 'window': The area of the window to be loaded. The first to points define the
  left corner of the window and the last two points define how far beyond to
  load in the respective direction.
"""
function loadSignals!(dl::DataLoader, window::NTuple{4,Int32})
    images = [loadTIF!(dl, fi, window) for fi in dl.files]
    images = [processData!(img, dl) for img in images]

    imgs = cat(images..., dims=4)
    imgs = permutedims(imgs, [1,2,4,3])

    replace!(imgs, NaN=>0.005)
    fill!(dl.signals, 0)

    xEnd = window[3] - dl.poolSz + 1
    yEnd = window[4] - dl.poolSz + 1
    i = 1
    for y in 1:dl.stride:yEnd, x in 1:dl.stride:xEnd
        poolSlice = imgs[x:x+dl.poolSz-1, y:y+dl.poolSz-1, :, 1]
        dl.signals[:,i] = maximum(poolSlice, dims=1:2)
        i += 1
    end
end
