using Distributed
@everywhere begin
    using Pkg; Pkg.activate(joinpath(@__DIR__, "../"))
end
Pkg.instantiate()
using Comonicon
using DrWatson
using DelimitedFiles
using DataFrames
using CSV
using ROSESoss

@everywhere include("rose_optimizer.jl")

"""
Run the mring optimizer on a list of hdf5 grmhd files

# Arguments

- `x`: A file that contains the paths to various GRMHd hdf5 files

# Options
- `--data <arg>`: The datafile you want to read in
- `--pa <arg>`: The position angles (deg) you want to rotate the images. To pass a list do e.g.  --pa "[0.0, 45.0, 90.0, 135.0]"
- `--out <arg>`: Where you want to save the output
- `--stride <arg>`: Checkpoint stride. This should be at least 2x the number of cores you are using.
"""
@main function main(x;
                    data=datadir("hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s_preprocessed_snapshot_60_noisefrac0.05_scan252.uvfits"),
                    pa::String = "[0.0]",
                    out=projectdir("_research/mring_grmhd.csv"),
                    stride::Int = 500
                   )

    println("\tParsed args:")
    println("\tfilelist => ", x)
    println("\tdata => ", data)
    println("\tpa => ", pa)
    println("\tstride => ", stride)
    println("\tout => ", out)
    println("Starting the run I currently have $(nworkers()) workers")

    # Parse in the pa angles
    pa_f = eval(Meta.parse(pa))
    println("Hello you are about use $pa_f pa angles")

    #Read in the file
    flist = open(x, "r") do io
        files = readlines(io)
    end

    println("I am about to analyze $(length(flist)) files")
    println("The first one is $(flist[1])")
    println("The last one is $(flist[end])")

    # Now I will construct an empty dataframe. This will be for checkpointing
    nfiles = length(flist)
    dftot = DataFrame()
    for pa in pa_f
        df = DataFrame(pa       = fill(pa, nfiles),
                       diam     = zeros(nfiles),
                       α        = zeros(nfiles),
                       ff       = zeros(nfiles),
                       fwhm_g    = zeros(nfiles),
                       amp1     = zeros(nfiles),
                       chi2_amp = zeros(nfiles),
                       chi2_cp = zeros(nfiles),
                       file     = flist
                       )

        pitr = Iterators.partition(eachindex(flist), stride)
        for p in pitr
            rows = pmap(p) do i
                fit_file(flist[i], data, pa)
            end
            @info "Checkpointing"
            df[p,1:end-1] = DataFrame(rows)
            push!(dftot, eachrow(df[p, 1:end])...)
            CSV.write(out, dftot)
        end
    end
    CSV.write(out ,dftot)

end
