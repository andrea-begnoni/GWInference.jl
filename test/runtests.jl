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

    CE_1= CE1Id
    CE_2 = CE2NM
    ET = ETS
    network = [CE_1, CE_2, ET];
    snrD_network = SNR(PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal) 
    @test isapprox(snrD_network, 37.695233406672514, rtol = 1e-12 )             

    snrHM_network = SNR(PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal)
    @test isapprox(snrHM_network, 37.92175372990184, rtol = 1e-12 )

    fisherD_network = FisherMatrix(
    PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    cov = CovMatrix(fisherD_network)
    errors = Errors(cov)

    errors_tabulated= [ 0.0039419796598880825
    0.03516266365602064
    0.6360722185017154
    1.644468273826038
    0.8085036923461478
    0.0335937819637836
    0.009880588299760065
    0.07298579934712374
    0.08559078006124807
    0.00800474002070679
    1.99701454428359]

    @test isapprox(errors, errors_tabulated, rtol = 1e-8)


    fisherHM_network = FisherMatrix(
        PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    covHM = CovMatrix(fisherHM_network)
    errorsHM = Errors(covHM)

    errors_tabulatedHM= [0.0013747374660239894
    0.014491816156939124
    0.28068667314616297
    0.7462978126898637
    0.7877776105708272
    0.0321076596579239
    0.009268865711735137
    0.07105591071390889
    0.08636127302050348
    0.0036638964937824718
    0.4224977988686636]

    @test isapprox(errorsHM, errors_tabulatedHM, rtol = 1e-8)
end
