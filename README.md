# EHTGRMHDCal


# Reproducing the environment
**Note I am assuming you have at least julia 1.6 installed**

Additionally you will require ehtim and dynesty to be installed in a system aware python. Currently I am just using the build in fft so 
`pip install ehtim`
`pip install dynesty`
should be enough. 

If for some reason Julia is looking at the wrong python installation please open a Julia REPL and type
```
ENV["PYTHON"] = "/path/to/python/binary"
using Pkg; Pkg.build()
```
This should allow you to use the custom python installation.



This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> EHTGRMHDCal

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

# Using the script
The main script you should look at in this repo is src/main.jl
To use it do something like
```
julia -p 2 main.jl filelist --data ../data/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s_preprocessed_snapshot_60_noisefrac0.05_scan252.uvfits --pa 90 --out test.csv  --stride 200
```

The only argument that isn't optional is `filelist`. This is a file that contains the paths of all the hdf5 file you
would like to analyze. For the other options please see the docstring of the main function.

**Note this will seem to hang at the begining. This is because Julia uses a JIT compiler which means it is called just before its first called**

# What about on clusters?
If you are using a single node then the -p option will work great. If you are using multiple nodes then it will fail! For a multiple node job there are two options. For a slurm cluster create a batch submision using

```
#!/bin/sh

#SBATCH .... #insert usual sbatch stuff

# create a host/machine file
srun hostname -s > hostfile

# now pass this machine file to julia
julia --machine-file ./hostfile main.jl filelist --data ../data/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s_preprocessed_snapshot_60_noisefrac0.05_scan252.uvfits --pa 90 --out test.csv  --stride 200

```
