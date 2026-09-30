#! format: off

##########################################################
# IMRPhenomD_NRTidalv2 WAVEFORM
##########################################################

# All is taken from LALSimulation and arXiv:1508.07250, arXiv:1508.07253, arXiv:1905.06011

"""
Needs documentation 
"""
function PolAbs(
    model::PhenomD_NRTidal,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    fInsJoin_PHI = 0.018,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc 
)

    amp = Ampl(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        Lambda1,
        Lambda2,
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc
    )

    # take into account inclination 
    hp = @. 0.5 * (1.0 + (cos(iota))^2) .* amp
    hc = @. cos(iota) .* amp

    return [hp, hc]
    
end

"""
ToDo: Need documentation
"""
function Pol(
    model::PhenomD_NRTidal,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    fInsJoin_PHI = 0.018,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc 
)

    hp, hc = PolAbs(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota,
        Lambda1,
        Lambda2;
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        fInsJoin_PHI = fInsJoin_PHI,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc 
    )

    return [hp, 1im .* hc]

end


"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomD_NRTidal(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    mc = 3.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    Lambda1 = 1000.
    Lambda2 = 2000.
    fcut = _fcut(PhenomD_NRTidal(), mc, eta, Lambda1, Lambda2)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomD_NRTidal(), f, mc, eta, dL, chi1, chi2, Lambda1, Lambda2)
```
"""
function Phi(model::PhenomD_NRTidal,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1,
    Lambda2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    interpolation = false,
)


    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    pi2 = pi * pi

    fInsJoin = fInsJoin_PHI

        # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
        # The relation between the two is given in arxiv:1608.02582 eq. (15) with coefficients from third row of Table I
        # We also extend the range to 0 <= Lam < 1, as done in LALSimulation in LALSimUniversalRelations.c line 123
    QuadMon1 = ifelse(Lambda1 < 1., 1. + Lambda1*(0.427688866723244 + Lambda1*(-0.324336526985068 + Lambda1*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda1) + 0.04812 * log(Lambda1)^2 -4.283e-3 * log(Lambda1)^3 + 1.245e-4 * log(Lambda1)^4))
    QuadMon2 = ifelse(Lambda2 < 1., 1. + Lambda2*(0.427688866723244 + Lambda2*(-0.324336526985068 + Lambda2*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda2) + 0.04812 * log(Lambda2)^2 -4.283e-3 * log(Lambda2)^3 + 1.245e-4 * log(Lambda2)^4))
    

    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2  = chi1*chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)       
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a
    chi_sdotchi_a  = chi_s*chi_a
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M*GMsun_over_c3.*f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = - 1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    
    if interpolation == true
        # Get the path to the directory of this file
        PACKAGE_DIR = @__DIR__

        # Go one step back in the path (from ""GW.jl/src" to "GW.jl")
        PARENT_DIR = dirname(dirname(PACKAGE_DIR))
        
        # Construct the path to the "useful_files" folder from the parent directory
        USEFUL_FILES_DIR = joinpath(PARENT_DIR, "useful_files/WFfiles/")

        QNMgrid_a = _readQNMgrid_a(USEFUL_FILES_DIR)
        QNMgrid_fring = _readQNMgrid_fring(USEFUL_FILES_DIR)
        QNMgrid_fdamp = _readQNMgrid_fdamp(USEFUL_FILES_DIR)
        fring = linear_interpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
        fdamp = linear_interpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)
    else
        fring = ((0.05947169566573468 - 0.14989771215394762*aeff + 0.09535606290986028*aeff*aeff + 0.02260924869042963*aeff*aeff*aeff - 0.02501704155363241*aeff*aeff*aeff*aeff - 0.005852438240997211*(aeff^5) + 0.0027489038393367993*(aeff^6) + 0.0005821983163192694*(aeff^7))/(1 - 2.8570126619966296*aeff + 2.373335413978394*aeff*aeff - 0.6036964688511505*aeff*aeff*aeff*aeff + 0.0873798215084077*(aeff^6)))/(1. - Erad)
        fdamp = ((0.014158792290965177 - 0.036989395871554566*aeff + 0.026822526296575368*aeff*aeff + 0.0008490933750566702*aeff*aeff*aeff - 0.004843996907020524*aeff*aeff*aeff*aeff - 0.00014745235759327472*(aeff^5) + 0.0001504546201236794*(aeff^6))/(1 - 2.5900842798681376*aeff + 1.8952576220623967*aeff*aeff - 0.31416610693042507*aeff*aeff*aeff*aeff + 0.009002719412204133*(aeff^6)))/(1. - Erad)
    end

    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a structure (spin effects are included up to 3.5PN)
    TF2OverallAmpl = 3 / (128.0 * eta)

    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2 *
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )


    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )

    #Now compute the coefficients to align the three parts

    fMRDJoin = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms

    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoin)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoin) / 3.0)
            ) * ((pi * fInsJoin)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoin)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoin)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoin) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoin)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoin)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoin)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoin^(1.0 / 3.0)) +
            sigma3 * (fInsJoin^(2.0 / 3.0)) +
            sigma4 * fInsJoin
        ) * etaInv

    # This is the first derivative of the Intermediate phase computed at fInsJoin
    DPhiInt = (beta1 + beta3 / (fInsJoin^4) + beta2 / fInsJoin) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoin^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoin^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoin^(1.0 / 3.0)) * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoin^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoin^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoin +
        PhiInspcoeffs.min_four_thirds * (fInsJoin^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoin^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoin +
            PhiInspcoeffs.four_thirds * (fInsJoin^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoin^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoin * fInsJoin
        ) * etaInv
    # This is the Intermediate phase computed at fInsJoin
    PhiIntJoin =
        beta1 * fInsJoin - beta3 / (3.0 * fInsJoin * fInsJoin * fInsJoin) +
        beta2 * log(fInsJoin)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoin

    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoin - beta3 / (3.0 * fMRDJoin^3) + beta2 * log(fMRDJoin)) * etaInv +
        C1Int +
        C2Int * fMRDJoin
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoin^4) + beta2 / fMRDJoin) * etaInv
    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoin^2) +
            alpha3 / (fMRDJoin^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) / (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoin) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoin^(0.75))) +
        alpha1 * fMRDJoin +
        alpha4 * atan((fMRDJoin - alpha5 * fring) / fdamp)

    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoin
    # Time shift so that peak amplitude is approximately at t=0
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi^2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi^2
        ) * xi
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs.(fring - (fdamp * gamma3) / gamma2),
        abs.(fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2^2)) * gamma3) / gamma2),
    )
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak^2) +
            alpha3 / (fpeak^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv

    # LAL sets fRef as the minimum frequency, do the same
    fRef = fgrid[1]  

    phiRef = ifelse(
        fRef < fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fRef^(2.0 / 3.0) +
        PhiInspcoeffs.third * fRef^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fRef^(1.0 / 3.0) * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.log * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.min_third * fRef^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fRef^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fRef +
        PhiInspcoeffs.min_four_thirds * fRef^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fRef^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fRef +
            PhiInspcoeffs.four_thirds * fRef^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fRef^(5.0 / 3.0) +
            PhiInspcoeffs.two * fRef * fRef
        ) * etaInv,
        ifelse(
            fRef < fMRDJoin,
            (beta1 * fRef - beta3 / (3.0 * fRef * fRef * fRef) + beta2 * log(fRef)) *
            etaInv +
            C1Int +
            C2Int * fRef,
            ifelse(
                fRef < fcutPar,
                (
                    -(alpha2 / fRef) +
                    (4.0 / 3.0) * (alpha3 * (fRef^(3.0 / 4.0))) +
                    alpha1 * fRef +
                    alpha4 * atan((fRef - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fRef,
                0.0,
            ),
        ),
    )
    
    phis = @. ifelse.(
        fgrid .< fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fgrid^(2.0 / 3.0) +
        PhiInspcoeffs.third * fgrid^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fgrid^(1.0 / 3.0) * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.log * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.min_third * fgrid^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fgrid^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fgrid +
        PhiInspcoeffs.min_four_thirds * fgrid^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fgrid^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fgrid +
            PhiInspcoeffs.four_thirds * fgrid^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fgrid^(5.0 / 3.0) +
            PhiInspcoeffs.two * fgrid * fgrid
        ) * etaInv,
        ifelse(
            fgrid < fMRDJoin,
            (beta1 * fgrid - beta3 / (3.0 * fgrid * fgrid * fgrid) + beta2 * log(fgrid)) * etaInv +
            C1Int +
            C2Int * fgrid,
            ifelse(
                fgrid < fcutPar,
                (
                    -(alpha2 / fgrid) +
                    (4.0 / 3.0) * (alpha3 * (fgrid^(3.0 / 4.0))) +
                    alpha1 * fgrid +
                    alpha4 * atan((fgrid - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fgrid,
                0.0,
            ),
        ),
    )

    # Add the tidal contribution to the phase, as in arXiv:1905.06011
    # Compute the tidal coupling constant, arXiv:1905.06011 eq. (8) using Lambda = 2/3 k_2/C^5 (eq. (10))

    kappa2T = (3.0/13.0) * ((1.0 + 12.0*m2ByM/m1ByM)*(m1ByM^5)*Lambda1 + (1.0 + 12.0*m1ByM/m2ByM)*(m2ByM^5)*Lambda2)
    
    c_Newt   = 2.4375
    n_1      = -12.615214237993088
    n_3over2 =  19.0537346970349
    n_2      = -21.166863146081035
    n_5over2 =  90.55082156324926
    n_3      = -60.25357801943598
    d_1      = -15.11120782773667
    d_3over2 =  22.195327350624694
    d_2      =   8.064109635305156

    numTidal = @.  1.0 + (n_1 * ((pi*fgrid)^(2. /3.))) + (n_3over2 * pi*fgrid) + (n_2 * ((pi*fgrid)^(4. /3.))) + (n_5over2 * ((pi*fgrid)^(5. /3.))) + (n_3 * (pi*fgrid)^2)
    denTidal = @. 1.0 + (d_1 * ((pi*fgrid)^(2. /3.))) + (d_3over2 * pi*fgrid) + (d_2 * ((pi*fgrid)^(4. /3.)))
    
    tidal_phase = @. - kappa2T * c_Newt / (m1ByM * m2ByM) * ((pi*fgrid)^(5. /3.)) * numTidal / denTidal
    
    # In the NRTidalv2 extension also 3.5PN spin-squared and 3.5PN spin-cubed terms are included, see eq. (27) of arXiv:1905.06011
    # This is needed to account for spin-induced quadrupole moments
    # Compute octupole moment and emove -1 to account for BBH baseline
    
    OctMon1 = - 1. + exp(0.003131 + 2.071 * log(QuadMon1)  - 0.7152 * log(QuadMon1)^2 + 0.2458 * log(QuadMon1)^3 - 0.03309 * log(QuadMon1)^4)
    OctMon2 = - 1. + exp(0.003131 + 2.071 * log(QuadMon2)  - 0.7152 * log(QuadMon2)^2 + 0.2458 * log(QuadMon2)^3 - 0.03309 * log(QuadMon2)^4)

    SS_3p5PN  = - 400. *pi*(QuadMon1-1.)*chi12*m1ByM^2 - 400. *pi*(QuadMon2-1.)*chi22*m2ByM^2
    SSS_3p5PN = 10. *((m1ByM^2+308. /3. *m1ByM)*chi1+(m2ByM^2-89. /3. *m2ByM)*chi2)*(QuadMon1-1.)*m1ByM^2*chi12 + 10. *((m2ByM^2+308. /3. *m2ByM)*chi2+(m1ByM^2-89. /3. *m1ByM)*chi1)*(QuadMon2-1.)*m2ByM^2*chi22 - 440. *OctMon1*m1ByM^3*chi12*chi1 - 440. *OctMon2*m2ByM^3*chi22*chi2

    return @. phis + ifelse(fgrid < fcutPar, - t0*(fgrid - fRef) - phiRef + tidal_phase + (SS_3p5PN + SSS_3p5PN)*TF2OverallAmpl*((pi*fgrid)^(2. /3.)), 0.)
end  

"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the source, in Gpc.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.
#### Optional arguments:
-  `fInsJoin`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is dimensionless.

#### Example:
```julia
    mc = 3.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    Lambda1 = 1000.
    Lambda2 = 2000.
    fcut = _fcut(PhenomD_NRTidal(), mc, eta, Lambda1, Lambda2)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2)
```
"""
function Ampl(model::PhenomD_NRTidal,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    ampl_notidal = Ampl(PhenomD(),
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = fcutPar,
    fInsJoin_Ampl = fInsJoin_Ampl
    )  
    
    M = mc/eta^(3. /5.)
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.) 
    fgrid = M * GMsun_over_c3 .* f
    Overallamp =
    2.0 * sqrt(5.0 / (64.0 * pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    # Now add the tidal amplitude as in arXiv:1905.06011
    # Compute the tidal coupling constant, arXiv:1905.06011 eq. (8) using Lambda = 2/3 k_2/C^5 (eq. (10))
    Xa = 0.5 * (1.0 + Seta)
    Xb = 0.5 * (1.0 - Seta)
    kappa2T = (3.0/13.0) * ((1.0 + 12.0*Xb/Xa)*(Xa^5)*Lambda1 + (1.0 + 12.0*Xa/Xb)*(Xb^5)*Lambda2)
    
    # Now compute the amplitude modification as in arXiv:1905.06011 eq. (24)
    xTidal = @. (pi * fgrid)^(2. /3.)
    n1T    = 4.157407407407407
    n289T  = 2519.111111111111
    dTidal = 13477.8073677
    polyTidal = @. (1.0 + n1T*xTidal + n289T*(xTidal^(2.89)))/(1. +dTidal*(xTidal^4))
    ampTidal = @. -9.0*kappa2T*(xTidal^3.25)*polyTidal
    
    # # Compute the dimensionless merger frequency (Mf) for the Planck taper filtering

    f_end_taper = _fcut(model, mc, eta, Lambda1, Lambda2) * GMsun_over_c3 * M
    f_merger = f_end_taper / 1.2
    planck_taper = @. ifelse(fgrid < f_merger, 1., ifelse(fgrid > f_end_taper, 0., 1. - 1. /(exp((f_end_taper - f_merger)/(fgrid - f_merger) + (f_end_taper - f_merger)/(fgrid - f_end_taper)) + 1.)))

    return  @. (ampl_notidal + 2*sqrt(pi/5.)*ampTidal*Overallamp)*planck_taper
end 



"""
Compute the cut frequency of the waveform as a function of the events parameters, in `Hz`.
Valid for PhenomD_NRTidal.
We cut the waveform slightly before the end of the Planck taper filter, for numerical stability.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.

#### Return:
-  (float) Cut frequency of the waveform, in Hz.

"""
function _fcut(model::PhenomD_NRTidal, mc, eta, Lambda1, Lambda2; GMsun_over_c3=uc.GMsun_over_c3)

    M = mc/eta^(3. /5.)
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    Xa = 0.5 * (1.0 + Seta)
    Xb = 0.5 * (1.0 - Seta)
    kappa2T = (3.0/13.0) * ((1.0 + 12.0*Xb/Xa)*(Xa^5)*Lambda1 + (1.0 + 12.0*Xa/Xb)*(Xb^5)*Lambda2)
    
    # Compute the dimensionless merger frequency (Mf) for the Planck taper filtering
    a_0 = 0.3586
    n_1 = 3.35411203e-2
    n_2 = 4.31460284e-5
    d_1 = 7.54224145e-2
    d_2 = 2.23626859e-4
    
    numPT = 1.0 + n_1*kappa2T + n_2*kappa2T^2
    denPT = 1.0 + d_1*kappa2T + d_2*kappa2T^2
    Q_0 = a_0 / sqrt(q)
    f_merger = Q_0 * (numPT / denPT) / (2. *pi)
    # Terminate the waveform at 1.2 times the merger frequency
    f_end_taper = 1.2*f_merger

    return f_end_taper/(M*GMsun_over_c3)
end

