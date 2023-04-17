using Comrade
using ComradeOptimization
using OptimizationMetaheuristics
import Distributions
const Dists = Distributions
using CSV
using DataFrames
using NamedTupleTools

const fwhmfac = 2*sqrt(2*log(2))


function mringwgfloor(θ)
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
            diam = Uniform(μas2rad(10.0), μas2rad(70.0)),
            fwhm = Uniform(μas2rad(1.0),  μas2rad(40.0)),
            ma   = ImageUniform(0.0, 0.5, nmodes, 1),
            mp   = DiagonalVonMises(fill(0.0, nmodes), fill(inv(π^2), nmodes)),
            floor= Uniform(0.0, 1.0),
            ϵ   = Uniform(1.0, 5.0)
        )
end

function create_post(model, data...)
    lklhd = RadioLikelihood(model, data...)
    prior = create_prior(model)
    return Posterior(lklhd, prior)
end


function loaddata(imfile, datafile, pa; f0=0.6, ferr=0.0)
    # Load the image
    img = ehtim.image.load_image(imfile)
    img.imvec *= f0/img.total_flux()

    # Assign a particular position angle [need to iterate over some values]
    img.pa = pa

    # Load the observation
    obs = ehtim.obsdata.load_uvfits(datafile)
    obs.add_fractional_noise(ferr)
    img.rf = obs.rf
    img.ra = obs.ra
    img.dec = obs.dec

    # Create a synthetic observation to fit
    obs_fit = img.observe_same(obs, ttype="fast", ampcal=true, phasecal=true, add_th_noise=true)
    #obs_fit.add_amp(debias=true)
    #obs_fit.add_cphase(count="min-cut0bl")
    damp = Comrade.extract_lcamp(obs_fit; snrcut=3.0)
    dcp = Comrade.extract_cphase(obs_fit; snrcut=3.0)
    return damp, dcp
end

function fit_file(imfile, datafile, pa; modes=1, model=mringwgfloor, maxevals=75_000)
    damp, dcp = loaddata(imfile, datafile, pa)
    post = create_post(model, damp, dcp)


    cpost = ascube(post)
    ndim = dimension(cpost)
    fopt  = OptimizationFunction(cpost)
    prob  = Optimization.OptimizationProblem(fopt, rand(ndim), nothing, lb=fill(0.001, ndim), ub = fill(0.999, ndim))
    sol   = solve(prob, ECA(), maxiters=maxevals)

    xopt = transform(cpost, sol.u)

    chi2amp = chi2(model(xopt), damp)/length(damp)
    chi2cp  = chi2(model(xopt), dcp)/length(dcp)


    df = (pa = pa, diam = xopt.diam,
                   α  = xopt.fwhm,
                   ff = xopt.floor,
                   fwhm_g= xopt.diam*xopt.ϵ,
                   amp1  = xopt.ma[1],
                   chi2_amp = chi2amp,
                   chi2_cp = chi2cp,
                   logp = -sol.minimum)
    return df
end
