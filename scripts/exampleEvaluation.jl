using Pkg
Pkg.activate("/path/to/ChangeDetectionPANN")

using ChangeDetectionPANN

name = "testExp"
path = "/path/to/Data/landslides"
evaluateChangeMaps(name, path)
