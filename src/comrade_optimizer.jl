using ComradeSoss
using BlackBoxOptim
using Metaheuristics
import Distributions
const Dists = Distributions
using HypercubeTransform
using CSV
using DataFrames
using NamedTupleTools
import Metaheuristics
const MH = Metaheuristics

const fwhmfac = 2*sqrt(2*log(2))



mringwgfloor = @model N begin
    diam ~ Dists.Uniform(25.0, 85.0)
    fwhm ~ Dists.Uniform(1.0, 40.0)
    rad = diam/2
    σ = fwhm/fwhmfac

    ma ~ Dists.Uniform(0.0, 0.5) |> iid(N)
    mp ~ Dists.Uniform(-1π, 1π) |> iid(N)
    α = ma.*cos.(mp)
    β = ma.*sin.(mp)

    #Fraction of floor flux
    floor ~ Dists.Uniform(0.0, 1.0)
    dg ~ Dists.Uniform(40.0, 300.0)
    rg = dg/fwhmfac
    mring = smoothed(stretched(Comrade.MRing{N}(α, β), rad, rad), σ)
    g = stretched(Comrade.Gaussian(), rg, rg)
    img = mring*(1-floor) + g*floor
    return img
end


vacp = @model image, uamp, vamp, erramp,
               u1cp, v1cp, u2cp, v2cp, u3cp, v3cp, errcp begin

    img ~ image
    vamps = Comrade.amplitudes(img, uamp, vamp)
    amp ~ For(eachindex(vamps, erramp)) do i
            Dists.Normal(vamps[i], erramp[i])
    end

    cphases = Comrade.closure_phases(img, u1cp, v1cp, u2cp, v2cp, u3cp, v3cp)
    cphase ~ Comrade.CPVonMises{(:μ,:σ)}(cphases, errcp)

end



function create_joint_nog(model,
    ampobs::Comrade.EHTObservation{F,A},
    cpobs::Comrade.EHTObservation{F,P};
    ) where {F, A<:Comrade.EHTVisibilityAmplitudeDatum,P<:Comrade.EHTClosurePhaseDatum}
    uamp = Comrade.getdata(ampobs, :u)
    vamp = Comrade.getdata(ampobs, :v)
    bl = Comrade.getdata(ampobs, :baseline)
    #stations = Tuple(unique(vcat(s1,s2)))
    #gpriors = values(select(amppriors, stations))
    erramp = Comrade.getdata(ampobs, :error)
    amps = Comrade.getdata(ampobs, :amp)

    u1cp = Comrade.getdata(cpobs, :u1)
    v1cp = Comrade.getdata(cpobs, :v1)
    u2cp = Comrade.getdata(cpobs, :u2)
    v2cp = Comrade.getdata(cpobs, :v2)
    u3cp = Comrade.getdata(cpobs, :u3)
    v3cp = Comrade.getdata(cpobs, :v3)
    errcp = Comrade.getdata(cpobs, :error)
    cps = Comrade.getdata(cpobs, :phase)

    joint = vacp(
    image=model,
    uamp=μas2rad(uamp),
    vamp=μas2rad(vamp),
    erramp=erramp,
    u1cp = μas2rad(u1cp),
    v1cp = μas2rad(v1cp),
    u2cp = μas2rad(u2cp),
    v2cp = μas2rad(v2cp),
    u3cp = μas2rad(u3cp),
    v3cp = μas2rad(v3cp),
    errcp = errcp
    )
    conditioned = (amp = amps, cphase = cps,)
    return joint | conditioned
end


function loaddata(imfile, datafile, pa; ferr=0.0)
    # Load the image
    img = ehtim.image.load_image(imfile)
    img.imvec /= img.total_flux()

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
    damp = Comrade.extract_amp(obs_fit; debias=true)
    dcp = Comrade.extract_cphase(obs_fit; count="min-cut0bl")
    return damp, dcp
end

function fit_file(imfile, datafile, pa; model=mringwgfloor(N=3,), maxevals=75_000)
    damp, dcp = loaddata(imfile, datafile, pa)
    cmg = create_joint_nog(model, damp, dcp)
    opt, stats = ComradeSoss.optimize(ComradeSoss.MetaH(alg=MH.ECA(N=100,options=MH.Options(f_calls_limit=maxevals))), cmg)

    bl = Comrade.getdata(damp, :baseline)
    s1 = first.(bl)
    s2 = last.(bl)
    stations = Tuple(unique(vcat(s1,s2)))
    gains = NamedTuple{stations}(Tuple(ones(length(stations))))

    mopt = Soss.predict(cmg.argvals[:image], opt[:img])

    chi2amp = chi2(mopt, damp, gains)/length(damp)
    chi2cp = chi2(mopt, dcp)/length(dcp)


    df = (pa = pa, diam = opt.img.diam,
                   α=opt.img.fwhm,
                   ff=opt.img.floor,
                   fwhm_g=opt.img.dg,
                   amp1 = opt.img.ma[1],
                   amp2 = opt.img.ma[2],
                   amp3 = opt.img.ma[3],
                   chi2_amp = chi2amp,
                   chi2_cp = chi2cp,
                   logp = stats)
    return df
end
