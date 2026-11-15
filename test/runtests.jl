# using Pkg
# Pkg.activate("../.")
using GWInference
using Test

@testset "GWInference.jl" begin
    mc = 10.
    eta = 0.20
    chi1 = 0.1

    chi2 = 0.2
    dL = 10.
    iota = 1.

    theta = .5
    phi = 1.
    psi = 0.5

    tcoal = 0.3
    phiCoal = 0.5

    # reference values below computed with the Jan 2026 10 km ET curves (ETS), on Julia 1.10
    CE_1= CE1Id
    CE_2 = CE2NM
    ET = ETS
    network = [CE_1, CE_2, ET];
    snrD_network = SNR(PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal) 
    @test isapprox(snrD_network, 39.156938855608196, rtol = 1e-12 )             

    snrHM_network = SNR(PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal)
    @test isapprox(snrHM_network, 39.3804482334247, rtol = 1e-12 )

    fisherD_network = FisherMatrix(
    PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    cov = CovMatrix(fisherD_network)
    errors = Errors(cov)

    errors_tabulated= [ 0.003863803861676538
    0.0352320058894503
    0.6380118491347128
    1.6488131983913583
    0.7599165204162439
    0.031222495948098748
    0.009927687954261081
    0.06814555417938505
    0.08105308080024505
    0.008026468912306394
    1.758869875790778]

    @test isapprox(errors, errors_tabulated, rtol = 1e-8)


    fisherHM_network = FisherMatrix(
        PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    covHM = CovMatrix(fisherHM_network)
    errorsHM = Errors(covHM)

    errors_tabulatedHM= [0.001358919968470554
    0.014311145689930372
    0.27704103196699753
    0.7363589820753051
    0.7431691635712657
    0.030007509839817426
    0.009287802246134305
    0.06664432087699251
    0.08177991690258951
    0.0036155116993000025
    0.40857162217027904]

    @test isapprox(errorsHM, errors_tabulatedHM, rtol = 1e-8)
end
