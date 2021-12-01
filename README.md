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
julia -p 2 main.jl filelist --out test.csv  --stride 200
```
**Note that by default the script will load the data in data/"snapshot_fitting_scans_120s_noisefrac0.05"**

The only argument that isn't optional is `filelist`. This is a file that contains the paths of all the hdf5 file you
would like to analyze. For the other options please see the docstring of the main function.

**Note this will seem to hang at the begining. This is because Julia uses a JIT compiler which means it is called just before its first called**

# What about on clusters?
If you are using a single node then the -p option will work great. If you are using multiple nodes then it will fail! For a multiple node job there are two options. 


## Slurm
For a slurm cluster create a batch submission using

```
#!/bin/sh

#SBATCH .... #insert usual sbatch stuff

# create a host/machine file
srun hostname -s > hostfile

# now pass this machine file to julia
julia --machine-file ./hostfile main.jl filelist --out test.csv  --stride 200

```

## Harvard Hydra

On Harvard's hydra cluster I have had some issues getting everything to work. The best way forward I found was to use the ClusterManagers.jl package to directly ask for resources. This is partly automated using src/main_hydra.jl. You will need to change one thing in this script for your resources. Specifically 

```
addprocs_sge(80; qsub_flags=`-l s_cpu=24:00:00 -l mres=4G,h_data=4G,h_vmem=4G`, wd=pwd())
```

This will ask hydra for 80 cores throughout the cluster. Please change 80 to however many cores you want.
You may want to ask for a single node to lower communication. 


