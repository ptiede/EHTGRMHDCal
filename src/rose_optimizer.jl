using ROSESoss
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
    fwhm ~ Dists.Uniform(1.0, 50.0)
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
    mring = smoothed(renormed(ROSE.MRing{N}(rad, α, β), (1-floor)), σ)
    g = renormed(stretched(ROSE.Gaussian(), rg, rg), floor)
    img = mring + g
    return img
end


vacp = @model image, uamp, vamp, erramp,
               u1cp, v1cp, u2cp, v2cp, u3cp, v3cp, errcp begin

    img ~ image
    vamps = ROSE.visibility_amplitude.(Ref(img), uamp, vamp)
    amp ~ For(eachindex(vamps, erramp)) do i
            Dists.Normal(vamps[i], erramp[i])
    end

    cphases = ROSE.closure_phase.(Ref(img), u1cp, v1cp, u2cp, v2cp, u3cp, v3cp)
    cphase ~ For(eachindex(cphases, errcp)) do i
        ROSE.CPVonMises(cphases[i], errcp[i])
    end

end



function create_joint_nog(model,
    ampobs::ROSE.EHTObservation{F,A},
    cpobs::ROSE.EHTObservation{F,P};
    amppriors=(AA=0.1,AP=0.1,AZ=0.1,LM=0.2,JC=0.1,PV=0.1,SM=0.1, SP=0.1)
    ) where {F, A<:ROSE.EHTVisibilityAmplitudeDatum,P<:ROSE.EHTClosurePhaseDatum}
    uamp = ROSE.getdata(ampobs, :u)
    vamp = ROSE.getdata(ampobs, :v)
    bl = ROSE.getdata(ampobs, :baselines)
    #stations = Tuple(unique(vcat(s1,s2)))
    #gpriors = values(select(amppriors, stations))
    erramp = ROSE.getdata(ampobs, :error)
    amps = ROSE.getdata(ampobs, :amp)

    u1cp = ROSE.getdata(cpobs, :u1)
    v1cp = ROSE.getdata(cpobs, :v1)
    u2cp = ROSE.getdata(cpobs, :u2)
    v2cp = ROSE.getdata(cpobs, :v2)
    u3cp = ROSE.getdata(cpobs, :u3)
    v3cp = ROSE.getdata(cpobs, :v3)
    errcp = ROSE.getdata(cpobs, :error)
    cps = ROSE.getdata(cpobs, :phase)

    joint = vacp(
    image=model,
    uamp=uamp,
    vamp=vamp,
    erramp=erramp,
    u1cp = u1cp,
    v1cp = v1cp,
    u2cp = u2cp,
    v2cp = v2cp,
    u3cp = u3cp,
    v3cp = v3cp,
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
    obs_fit = img.observe_same(obs, ttype="fast", ampcal=true, phasecal=true, add_th_noise=false)
    obs_fit.add_amp(debias=true)
    obs_fit.add_cphase(count="min")
    damp = ROSESoss.extract_amps(obs_fit)
    dcp = ROSESoss.extract_cphase(obs_fit)
    return damp, dcp
end

function fit_file(imfile, datafile, pa; model=mringwgfloor(N=3,), maxevals=75_000)
    damp, dcp = loaddata(imfile, datafile, pa)
    cmg = create_joint_nog(model, damp, dcp)
    opt, stats = ROSESoss.optimize(ROSESoss.MetaH(alg=MH.ECA(N=100,options=MH.Options(f_calls_limit=maxevals))), cmg)

    bl = ROSE.getdata(damp, :baselines)
    s1 = first.(bl)
    s2 = last.(bl)
    stations = Tuple(unique(vcat(s1,s2)))
    gains = NamedTuple{stations}(Tuple(ones(length(stations))))

    mopt = Soss.predict(cmg.argvals[:image], opt[:img])

    chi2amp = chi2(mopt, damp, gains)/ROSE.nsamples(damp)
    chi2cp = chi2(mopt, dcp)/ROSE.nsamples(dcp)


    df = (pa = pa, diam = opt.img.diam,
                   α=opt.img.fwhm,
                   ff=opt.img.floor,
                   fwhm_g=opt.img.dg,
                   amp1 = opt.img.ma[1],
                   chi2_amp = chi2amp,
                   chi2_cp = chi2cp)
    return df
end
