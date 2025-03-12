using Comrade
using Optimization
using OptimizationMetaheuristics
using Distributions
using CSV
using DataFrames
using NamedTupleTools
using VLBIImagePriors
using Pyehtim
const fwhmfac = 2*sqrt(2*log(2))


function mringwgfloor(θ, meta)
    (;diam, fwhm, ma, mp, floor, ϵ) = θ

    rad = diam/2
    σ   = fwhm/fwhmfac
    α = ma.*cos.(mp)
    β = ma.*sin.(mp)
    ring = smoothed(modify(MRing(α, β), Stretch(rad, rad)), σ)

    rg   = diam/fwhmfac*ϵ
    gaus = modify(Gaussian(), Stretch(rg, rg))

    m = ring*(1-floor) + gaus*floor
    return m
end

function create_prior(::typeof(mringwgfloor), nmodes)
    return (
            diam = Uniform(μas2rad(30.0), μas2rad(70.0)),
            fwhm = Uniform(μas2rad(1.0),  μas2rad(40.0)),
            ma   = ntuple(_->Uniform(0.0, 0.5), nmodes),
            mp   = ntuple(_->Uniform(0.0, 2π), nmodes),
            floor= Uniform(0.0, 1.0),
            ϵ   = Uniform(1.0, 5.0)
        )
end


function loaddata(imfile, datafile, pa; f0=0.6, ferr=0.0)
    # Load the image
    img = ehtim.image.load_image(imfile)
    img.imvec *= f0/img.total_flux()

    # Assign a particular position angle [need to iterate over some values]
    img.pa = pa

    # Load the observation
    obs = scan_average(ehtim.obsdata.load_uvfits(datafile))
    obs.add_fractional_noise(ferr)
    img.rf = obs.rf
    img.ra = obs.ra
    img.dec = obs.dec

    # Create a synthetic observation to fit
    obs_fit = img.observe_same(obs, ttype="fast", ampcal=true, phasecal=true, add_th_noise=true)
    dcp, damp = extract_table(obs_fit, ClosurePhases(;snrcut=3.0), VisibilityAmplitudes())
    return damp, dcp
end

function fit_file(imfile, datafile, pa; modes=4, model=mringwgfloor, maxevals=200_000)
    damp, dcp = loaddata(imfile, datafile, pa)
    sky = SkyModel(model, create_prior(model, modes), imagepixels(μas2rad(100.0), μas2rad(100.0), 256, 256))
    post = VLBIPosterior(sky, damp, dcp)

    xopt, sol = comrade_opt(post, ECA(); maxiters=30_000)

    ndim = dimension(post)
    chi2amp, chi2cp = chi2(post, xopt)
    rchi2   = (chi2amp + chi2cp)/(length(damp) + length(dcp) - ndim)


    df = (pa = pa, diam = xopt.sky.diam,
                   α  = xopt.sky.fwhm,
                   ff = xopt.sky.floor,
                   fwhm_g= xopt.sky.diam*xopt.sky.ϵ,
                   amp1  = xopt.sky.ma[1],
                   chi2_amp = chi2amp/length(damp),
                   chi2_cp = chi2cp/length(dcp),
                   chi2    = rchi2,
                   logp = -sol.minimum)
    return df
end
