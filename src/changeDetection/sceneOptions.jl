"Struct that contains the information to process a natural disaster scene."
mutable struct SceneOptions
    "Name of the experiment and location of the natural disaster scene."
    name::String
    "Determines if the output values from the networks are saved."
    saveFiles::Bool
    "Pixel size of the tiles."
    tileSz::Tuple{Int32, Int32}
    "Number of images used before the natural disaster event."
    history::Int32
    "Seed used to run the model to process the scene."
    seed::Int32
    "Number of tiles per image."
    numTiles::Int32
    "Pixel width of the images in the scene."
    imgWidth::Int32
    "Pixel height of the images in the scene."
    imgHeight::Int32
    "List of .tif files used in the scene."
    files::Vector{String}
    "Size of the change map created for the scene."
    changeMapSz::Tuple{Int32, Int32}
end

"""
    updateLocation!(so::SceneOptions, name::String, location::String,
                    path::String)

Updates the values of the scene options for the new 'location'.

# Arguments
- 'so': The SceneOptions struct to be updated.
- 'name': Name of the experiment being run.
- 'location': New location to be updating the 'so' for.
- 'path': Path to the directory containing the .tif files for the 'location'.
"""
function updateLocation!(so::SceneOptions, name::String, location::String,
                        path::String)
    so.name = name*"_"*location

    files = readdir(path)
    files = [joinpath(path, fi) for fi in files]
    sort!(files)
    files = files[end-so.history:end]
    so.files .= files

    handle = AG.read(files[1])
    so.imgHeight = AG.height(handle)
    so.imgWidth = AG.width(handle)
    so.changeMapSz = (trunc(Int, so.imgHeight/so.tileSz[1]), 
                      trunc(Int, so.imgWidth/so.tileSz[2]))
    so.numTiles = so.changeMapSz[1] * so.changeMapSz[2]
end
