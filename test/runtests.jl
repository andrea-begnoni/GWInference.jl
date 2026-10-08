
using GW
using Test

@testset "GW.jl" begin
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
    @test isapprox(snrD_network, 37.69523553190364, rtol = 1e-12 )             

    snrHM_network = SNR(PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal)
    @test isapprox(snrHM_network, 37.92175372990184, rtol = 1e-12 )

    fisherD_network = FisherMatrix(
    PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    cov = CovMatrix(fisherD_network)
    errors = Errors(cov)

    errors_tabulated= [0.003941783348926305
    0.03516108811431253
    0.6360420432561998
    1.644391008835697
    0.8085029521022669
    0.03359405656156905
    0.009880704766380628
    0.07298583788328453
    0.0855908247386635
    0.008003990135235535
    1.9969726122385671]

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

include("test_phenomxphm.jl")

@testset "PhenomXE" begin
    # Eccentric aligned-spin 22-mode model (port of phenomxpy IMRPhenomXE). The dynamics are integrated with an
    # adaptive ODE solver whose roundoff-level error estimate makes the accepted-step sequence (hence the waveform
    # at the ~1e-4 level) machine dependent, so the tolerances below are looser than for the other models.
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
    ecc = 0.1
    meanAno = 1.0
    network = [CE1Id, CE2NM, ETS]

    @test _npar(PhenomXE()) == 13
    @test _available_waveforms("PhenomXE") isa PhenomXE

    f = [10.0, 20.0, 50.0]
    hp, hc = hphc(PhenomXE(), f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno)
    hp_tabulated = [-1.5089773722822782e-24 - 1.5704287766696464e-24im, -2.964057299671978e-25 - 8.426930779708367e-25im, 2.0315512002982394e-25 + 2.7670630692530225e-25im]
    hc_tabulated = [-1.313551870884822e-24 + 1.2621521459175385e-24im, -7.048527673427448e-25 + 2.47922291620931e-25im, 2.3144512667308263e-25 - 1.6992479503649662e-25im]
    @test isapprox(hp, hp_tabulated, rtol = 1e-3)
    @test isapprox(hc, hc_tabulated, rtol = 1e-3)
    # the e -> 0 waveform has the same amplitude as IMRPhenomXAS
    hp0, hc0 = hphc(PhenomXE(), f, mc, eta, chi1, chi2, dL, iota, 0.0, 0.0)
    ampXAS = Ampl(PhenomXAS(), f, mc, eta, chi1, chi2, dL)
    @test isapprox(abs.(hp0), 0.5 .* (1 .+ cos(iota)^2) .* ampXAS, rtol = 1e-3)

    snrXE_network = SNR(PhenomXE(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, ecc, meanAno)
    @test isapprox(snrXE_network, 37.574351252981, rtol = 1e-3)

    fisherXE_network = FisherMatrix(
        PhenomXE(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, ecc, meanAno, coordinate_shift=false)
    @test size(fisherXE_network) == (13, 13)
    @test all(isfinite, fisherXE_network)
    errorsXE = Errors(CovMatrix(fisherXE_network))

    errors_tabulatedXE = [0.0013330044268134254
    0.0018860454311576797
    0.0917927950405101
    0.31881181456740976
    0.9092385615164901
    0.03447770082595458
    0.009399000345739641
    0.08164627135721421
    0.08437607203758284
    0.00016348878329178047
    7.497095635131257
    0.0004495179316747921
    2.8420654696997247]

    @test isapprox(errorsXE, errors_tabulatedXE, rtol = 1e-1)
end
