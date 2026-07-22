using LinearAlgebra

const _LAL_XPHM_ORACLE = joinpath(@__DIR__, "lal_xphm_reference.py")
const _LAL_XPHM_PYTHONPATH = get(
    ENV,
    "LALSUITE_PYTHONPATH",
    joinpath(homedir(), "lalsuite-install"),
)

function _lal_xphm(case, frequencies)
    arguments = string.((
        case.m1, case.m2,
        case.chi1...,
        case.chi2...,
        case.distance, case.iota, case.phiRef, case.fRef,
        frequencies...,
    ))
    command = addenv(
        `python3 $(_LAL_XPHM_ORACLE) $arguments`,
        "PYTHONPATH" => _LAL_XPHM_PYTHONPATH,
    )
    rows = split(chomp(read(command, String)), '\n')
    values = [parse.(Float64, split(row)) for row in rows]
    hp = ComplexF64[row[1] + im*row[2] for row in values]
    hc = ComplexF64[row[3] + im*row[4] for row in values]
    return hp, hc
end

function _julia_xphm(case, frequencies)
    total_mass = case.m1 + case.m2
    eta = case.m1 * case.m2 / total_mass^2
    mc = total_mass * eta^(3/5)
    return hphc(
        PhenomXPHM(), frequencies, mc, eta,
        case.chi1, case.chi2, case.distance, case.iota;
        phiRef=case.phiRef, fRef=case.fRef,
    )
end

@testset "IMRPhenomXPHM against LALSuite" begin
    @test _available_waveforms("PhenomXPHM") isa PhenomXPHM
    @test _npar(PhenomXPHM()) == 15
    @test_throws ArgumentError hphc(
        PhenomXPHM(), [20.0], 30.0, 0.24,
        (0.1, 0.2, 0.3), (0.0, 0.1, -0.2), 1.0, 0.7;
        prec_version=222,
    )

    cases = [
        (
            m1=40.0, m2=20.0,
            chi1=(0.30, 0.20, 0.40), chi2=(-0.10, 0.25, -0.20),
            distance=1.0, iota=0.70, phiRef=0.0, fRef=20.0,
            frequencies=[20.0, 30.0, 50.0, 80.0],
        ),
        (
            m1=35.0, m2=30.0,
            chi1=(0.05, -0.30, -0.40), chi2=(0.20, 0.15, 0.60),
            distance=1.4, iota=2.20, phiRef=-0.40, fRef=25.0,
            frequencies=[20.0, 25.0, 45.0, 75.0],
        ),
    ]

    for case in cases
        hp_lal, hc_lal = _lal_xphm(case, case.frequencies)
        hp_julia, hc_julia = _julia_xphm(case, case.frequencies)
        plus_relative_error = norm(hp_julia-hp_lal) / norm(hp_lal)
        cross_relative_error = norm(hc_julia-hc_lal) / norm(hc_lal)
        @test plus_relative_error < 1e-8
        @test cross_relative_error < 1e-8
    end
end
