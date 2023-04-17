
using Distributed
@everywhere begin
    using Pkg; Pkg.activate(joinpath(@__DIR__, "../"))
end
#try
#    Pkg.resolve()
#    Pkg.update()
#catch
#    println("I am guessing you are using Julia 1.5 using a hack to get around dependency issues")
#    Pkg.rm("Comrade")
#    Pkg.rm("ComradeSoss")
#    Pkg.instantiate()
#    Pkg.add(url="https://github.com/ptiede/Comrade.jl")
#    Pkg.add(url="https://github.com/ptiede/ComradeSoss.jl")
#end
Pkg.instantiate()
using Comonicon
using DelimitedFiles
using DataFrames
using CSV
using ComradeSoss

load_ehtim()

@everywhere include("comrade_optimizer.jl")

function rundata(flist, pa_f, data, out, stride)

    # Now I will construct an empty dataframe. This will be for checkpointing
    nfiles = length(flist)
    dftot = DataFrame()
    for pa in pa_f
        df = DataFrame(pa       = fill(pa, nfiles),
                       diam     = zeros(nfiles),
                       α        = zeros(nfiles),
                       ff       = zeros(nfiles),
                       fwhm_g   = zeros(nfiles),
                       amp1     = zeros(nfiles),
                       amp2     = zeros(nfiles),
                       amp3     = zeros(nfiles),
                       chi2_amp = zeros(nfiles),
                       chi2_cp  = zeros(nfiles),
                       logp     = zeros(nfiles),
                       file     = flist
                       )

        pitr = Iterators.partition(eachindex(flist), stride)
        for p in pitr
            rows = map(p) do i
                fit_file(flist[i], data, pa)
            end
            @info "Checkpointing"
            df[p,1:end-1] = DataFrame(rows)
            push!(dftot, eachrow(df[p, 1:end])...)
            CSV.write(out, dftot)
        end
    en
    CSV.write(out, dftot)

end


"""
Run the mring optimizer on a list of hdf5 grmhd files

# Arguments

- `x`: A file that contains the paths to various GRMHd hdf5 files

# Options
- `--data <arg>`: path to the directory where the uvfits files for each scan are. This should work by default if you are running this in src
- `--pa <arg>`: The position angles (deg) you want to rotate the images. To pass a list do e.g.  --pa "[0.0, 45.0, 90.0, 135.0]"
- `--out <arg>`: Where you want to save the output
- `--stride <arg>`: Checkpoint stride. This should be at least 2x the number of cores you are using.
"""
@main function main(x;
                    data=datadir("snapshot_fitting_scans_120s_noisefrac0.05"),
                    pa::String = "[0.0]",
                    out="mring_fits.csv",
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
    #Read in the list of data scans
    dlist = open(x, "r") do io
        readlines(io)
    end


    println("I am about to analyze $(length(flist)) files")
    println("The first one is $(flist[1])")
    println("The last one is $(flist[end])")

    # Now loop over the list of data
    dlist = filter(endswith(".uvfits"), readdir(data, join=true))
    for d in dlist
        scan = first(splitext(split(d, "_")[end]))
        println("On $scan")
        base, ext = splitext(out)
        outscan = base*"_"*scan*ext
        rundata(flist, pa_f, d, outscan, stride)
    end

end
