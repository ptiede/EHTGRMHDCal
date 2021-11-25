#!/bin/bash -l
# -*- mode: julia -*-
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:05:00
#SBATCH --output=%x_%j.log
#SBATCH --job-name=precompile_everywhere

# Edit the SBATCH directives above as needed, then submit this script as a
# batch job with something like `sbatch precompile.jl`

#=
# This block will execute in bash.  The lines with "${BASH_SOURCE[0]}" will run
# the remainder of this file as a julia script with the given arguments,
# but should return and proceed with the other commands, eventually generating
# `sys_everywhere.so`.
# Echo each command in bash
set -x
# Exit on any error
set -e
# Prepare the project; I don't know why these are all used, but they are
julia -e 'using Pkg; Pkg.add(["Distributed", "ClusterManagers", "Sockets", "Serialization", "Logging", "LinearAlgebra", "REPL"])'
# Run @everywhere 1+2 on 2 processes and trace
julia --trace-compile=precompile01.jl "${BASH_SOURCE[0]}" everywhere
# Combine the precompilation files into one
echo "using Distributed, ClusterManagers, Sockets, Serialization, Logging, LinearAlgebra, REPL" > precompile_everywhere.jl
cat precompile01.jl precompile02.jl >> precompile_everywhere.jl
# Create the sysimage
julia "${BASH_SOURCE[0]}" precompile
# Finally, just print out instructions
exec cat << EOF
The sysimage `sys_everywhere.so` has now been built.  To use it, run julia as
    julia --sysimage $PWD/sys_everywhere.so [other switches] [program file] [args...]
You should be able to move the .so file around on this system, but you probably
won't be able to move it to different systems.
EOF
=#

using PackageCompiler
using Distributed
using ClusterManagers
using Sockets
using Serialization
using Logging
using LinearAlgebra
using REPL

if "everywhere" ∈ ARGS

    addprocs_slurm(
        1;
        partition=get(ENV, "SBATCH_PARTITION", "development"),
        time="00:04:50",
        exeflags="--trace-compile=precompile02.jl"
    )
    while nprocs() < 2
        sleep(1)
    end
    @everywhere 1+2
    map(rmprocs, workers())

elseif "precompile" ∈ ARGS

    create_sysimage(
        [:Distributed, :ClusterManagers, :Sockets, :Serialization, :Logging, :LinearAlgebra, :REPL],
        sysimage_path="sys_everywhere.so",
        precompile_execution_file="precompile_everywhere.jl"
    )

else

    ArgumentError("""Accepted args are "everywhere" or "precompile"; got $ARGS""")

end