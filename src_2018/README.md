# Running the 2018 analysis script

To run the 2018 analysis script is very similar to the older version in the main repo. 

First installation

First change in the parent `EHTGRMHDCal` folder launch Julia by doing 

```
julia --project=@.
```

Then in the julia REPL do
```julia
julia> using Pkg
julia> Pkg.instantiate()
```

which will make sure all the dependencies are installed.

Next to use the script from the bash terminal we (assuming we are in the `src_2018` directory)
```
julia -p NCORES main.jl /path/to/file/with/list/of/image/paths --data /path/to/directory/with/data --pa "[0.0, 45.0, 90.0, 135.0]" --out mring_results.csv --stride 500
```

 -  first argument to the `main.jl` script is the path to the file that contains all the fits or hdf5 GRMHD snapshots.
 -  `--data` argument is the directory that contains all different uvfits coverages that we want to fit to
 -  `--pa` are the differnt PA's of the GRMHD simulation that we want to consider. For M87 for example, this could be limited to the 3mm M87 jet PA
 -  `--out` the CSV file that we will save the results to
 -  `--stride` the stride of the parallelization. This should be probably at least 2x the number of cores. 