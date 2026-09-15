using Distributions
using GW
using Random
using Statistics
using Test

const LAL_XPHM_ORACLE = joinpath(@__DIR__, "lal_xphm_reference.py")
const LAL_XPHM_PYTHONPATH = get(
    ENV,
    "LALSUITE_PYTHONPATH",
    joinpath(homedir(), "lalsuite-install"),
)

"""Draw a spin direction using the BBH tilt prior in `src/catalog.jl`."""
function sample_spin(rng, magnitude, tilt_prior)
    cos_tilt = rand(rng, tilt_prior)
    azimuth = 2pi * rand(rng)
    sin_tilt = sqrt(1-cos_tilt^2)
    return (
        magnitude*sin_tilt*cos(azimuth),
        magnitude*sin_tilt*sin(azimuth),
        magnitude*cos_tilt,
    )
end

"""
Draw reproducible validation events from broad BBH priors.

The spin-magnitude and tilt distributions match the BBH population prior in
`src/catalog.jl`; spin azimuth, orientation, phase, masses, and distance use
the stated broad validation priors.  Component masses satisfy m1 >= m2 and
mass ratio q >= 1/8.
"""
function sample_xphm_prior(number_events; seed=0x5850484d)
    rng = Xoshiro(seed)
    spin_magnitude_prior = Beta(1.6, 4.12)
    tilt_prior = MixtureModel(
        [truncated(Normal(0.0, 1.5), -1.0, 1.0), Uniform(-1.0, 1.0)],
        [0.66, 0.34],
    )

    return map(1:number_events) do event
        m1 = rand(rng, Uniform(10.0, 87.0))
        q = rand(rng, Uniform(max(5.1/m1, 1/8), 1.0))
        m2 = q*m1
        chi1 = sample_spin(rng, rand(rng, spin_magnitude_prior), tilt_prior)
        chi2 = sample_spin(rng, rand(rng, spin_magnitude_prior), tilt_prior)
        total_mass = m1+m2
        fmax = min(1024.0, 0.28/(total_mass*4.925491025543576e-6))
        frequencies = collect(10.0 .^ range(log10(20.0), log10(fmax), length=32))
        frequencies[1] = 20.0
        (
            id=event,
            m1=m1,
            m2=m2,
            chi1=chi1,
            chi2=chi2,
            distance=exp(rand(rng, Uniform(log(0.1), log(5.0)))),
            iota=acos(rand(rng, Uniform(-1.0, 1.0))),
            phiRef=rand(rng, Uniform(-pi, pi)),
            fRef=20.0,
            frequencies=frequencies,
        )
    end
end

function lal_xphm(case)
    arguments = string.((
        case.m1, case.m2,
        case.chi1...,
        case.chi2...,
        case.distance, case.iota, case.phiRef, case.fRef,
        case.frequencies...,
    ))
    command = addenv(
        `python3 $(LAL_XPHM_ORACLE) $arguments`,
        "PYTHONPATH" => LAL_XPHM_PYTHONPATH,
    )
    values = [parse.(Float64, split(row))
              for row in split(chomp(read(command, String)), '\n')]
    hp = ComplexF64[row[1] + im*row[2] for row in values]
    hc = ComplexF64[row[3] + im*row[4] for row in values]
    return hp, hc
end

function julia_xphm(case)
    total_mass = case.m1+case.m2
    eta = case.m1*case.m2/total_mass^2
    mc = total_mass*eta^(3/5)
    return hphc(
        PhenomXPHM(), case.frequencies, mc, eta,
        case.chi1, case.chi2, case.distance, case.iota;
        phiRef=case.phiRef, fRef=case.fRef,
    )
end

function envelope_differences(case)
    hp_lal, hc_lal = lal_xphm(case)
    hp_julia, hc_julia = julia_xphm(case)

    plus = abs.(hp_julia-hp_lal) ./ maximum(abs, hp_lal)
    cross = abs.(hc_julia-hc_lal) ./ maximum(abs, hc_lal)
    lal_envelope = sqrt.(abs2.(hp_lal) .+ abs2.(hc_lal))
    residual_envelope = sqrt.(abs2.(hp_julia-hp_lal) .+ abs2.(hc_julia-hc_lal))
    joint = residual_envelope ./ maximum(lal_envelope)
    return (; plus, cross, joint)
end

@testset "IMRPhenomXPHM 100-event prior sweep against LALSuite" begin
    cases = sample_xphm_prior(100)
    differences = Vector{NamedTuple}(undef, length(cases))

    # Keep this sequential: concurrent evaluation produced nondeterministic
    # comparison results and therefore is not suitable for validation.
    for event in eachindex(cases)
        differences[event] = envelope_differences(cases[event])
    end

    plus = reduce(vcat, getproperty.(differences, :plus))
    cross = reduce(vcat, getproperty.(differences, :cross))
    joint = reduce(vcat, getproperty.(differences, :joint))
    worst_event = argmax(maximum.(getproperty.(differences, :joint)))
    worst_frequency = argmax(differences[worst_event].joint)

    @test all(isfinite, plus)
    @test all(isfinite, cross)
    @test all(isfinite, joint)
    tolerance = parse(Float64, get(ENV, "GW_XPHM_SWEEP_TOLERANCE", "1e-3"))
    @test maximum(plus) < tolerance
    @test maximum(cross) < tolerance
    @test maximum(joint) < tolerance

    println("IMRPhenomXPHM prior sweep ($(length(cases)) events, $(length(joint)) samples, seed=0x5850484d):")
    println("  plus:  maximum=$(maximum(plus)), median=$(median(plus))")
    println("  cross: maximum=$(maximum(cross)), median=$(median(cross))")
    println("  joint: maximum=$(maximum(joint)), median=$(median(joint))")
    println("  worst joint sample: event=$(cases[worst_event].id), frequency=$(cases[worst_event].frequencies[worst_frequency]) Hz")
    println("  worst event parameters: $(cases[worst_event])")
end
