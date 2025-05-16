using Pkg
Pkg.activate("path/to/ChangeDetectionPANN")

using ChangeDetectionPANN

function main()

     path = "/path/to/Data/"
     connectivityFile = "803nw12279junctions256electrodes.mat"

     runScenes(connectivityFile, 16, 803, 400, "testExp", path, 
               enableProgBar=false, saveFiles=true)
end

main()
