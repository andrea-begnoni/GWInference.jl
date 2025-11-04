#! format: off

#############################################################
# IMRPhenomD WAVEFORM
#############################################################

#function _readQNMgrid_a(pathWF::String)
#    return readdlm(pathWF * "QNMData_a.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
#end
#
#function _readQNMgrid_fring(pathWF::String)
#    return readdlm(pathWF * "QNMData_fring.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
#end
#
#function _readQNMgrid_fdamp(pathWF::String)
#    return readdlm(pathWF * "QNMData_fdamp.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
#end

"""
ToDo: Need documentation
"""
function PolAbs(model::PhenomD,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc,
        container = container
    )

    # take into account inclination 
    hp = @. 0.5 * (1.0 + (cos(iota))^2) .* amp
    hc = @. cos(iota) .* amp

    # return plus and cross polarization absolute values
    return [hp, hc]
    
end

"""
ToDo: Need documentation
"""
function Pol(model::PhenomD,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing
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
        Lambda2,
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc,
        container = container 
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

## # Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude
# fInsJoin_Ampl = 0.014
# # Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase
# fInsJoin_PHI = 0.018
# # Dimensionless frequency (Mf) at which we define the end of the waveform
# fcutPar = 0.2
"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomD(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomD(), f, mc, eta, chi1, chi2)
```
"""
function Phi(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    container = nothing,
    fRef_hz = 20.0
)
    println("Using $fRef_hz Hz as reference frequency for PhenomD phase computation.")

    # Get the path to the directory of this file
    PACKAGE_DIR = @__DIR__

    # Go one step back in the path (from ""GW.jl/src" to "GW.jl")
    PARENT_DIR = dirname(dirname(PACKAGE_DIR))
    
    # Construct the path to the "useful_files" folder from the parent directory
    USEFUL_FILES_DIR = joinpath(PARENT_DIR, "useful_files/WFfiles/")

    QNMgrid_a = _readQNMgrid_a(USEFUL_FILES_DIR)
    QNMgrid_fring = _readQNMgrid_fring(USEFUL_FILES_DIR)
    QNMgrid_fdamp = _readQNMgrid_fdamp(USEFUL_FILES_DIR)


    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    pi2 = pi * pi

    fInsJoin = fInsJoin_PHI
    QuadMon1, QuadMon2 = 1.0, 1.0

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)  
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f 
    fRef = M * GMsun_over_c3 * fRef_hz
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = linear_interpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
    fdamp = linear_interpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)

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

    # For 3PN coeff we use chi1 and chi2 so to have the quadrupole moment explicitly appearing


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
    #fRef = fgrid[1] 

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
    if typeof(eta) == Float64
       phi = @. phis + ifelse(fgrid < fcutPar, -t0 * (fgrid - fRef) - phiRef, 0.0)
       if container !== nothing
        container .= [phi[i] for i in eachindex(phi)]
        end
    else 
        phi = @. phis + ifelse(fgrid .< fcutPar, -t0 * (fgrid - fRef) - phiRef, ForwardDiff.Dual{typeof(eta).parameters[1]}(0.,zeros(typeof(eta).parameters[3])...))
        if container !== nothing
            container .= [phi[i].value for i in eachindex(phi)]
        end
    end


    return phi
end


"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the binary, in Gpc.
#### Optional arguments:
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is dimensionless.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)
```
"""
function Ampl(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing,
)

    # Get the path to the directory of this file
    PACKAGE_DIR = @__DIR__

    # Go two steps back in the path (from ""GW.jl/src/waveforms" to "GW.jl")
    PARENT_DIR = dirname(dirname(PACKAGE_DIR))
    
    # Construct the path to the "useful_files" folder from the parent directory
    USEFUL_FILES_DIR = joinpath(PARENT_DIR, "useful_files/WFfiles/")

    QNMgrid_a = _readQNMgrid_a(USEFUL_FILES_DIR)
    QNMgrid_fring = _readQNMgrid_fring(USEFUL_FILES_DIR)
    QNMgrid_fdamp = _readQNMgrid_fdamp(USEFUL_FILES_DIR)

    # Useful quantities
    M = mc / (eta^(3.0 / 5.0))
    eta2 = eta * eta # This can speed up a bit, we call it multiple times
    pi2 = pi * pi
    chi12, chi22 = chi1 * chi1, chi2 * chi2

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = linear_interpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
    fdamp = linear_interpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)
    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
    xi2 = xi * xi
    gamma1 =
        0.006927402739328343 +
        0.03020474290328911 * eta +
        (
            0.006308024337706171 - 0.12074130661131138 * eta +
            0.26271598905781324 * eta2 +
            (
                0.0034151773647198794 - 0.10779338611188374 * eta +
                0.27098966966891747 * eta2
            ) * xi +
            (
                0.0007374185938559283 - 0.02749621038376281 * eta +
                0.0733150789135702 * eta2
            ) * xi2
        ) * xi
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi2
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
            ) * xi2
        ) * xi
    # Compute fpeak, from arXiv:1508.07253 eq. (20), we remove the square root term in case it is complex
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs(fring - (fdamp * gamma3) / gamma2),
        fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2 * gamma2)) * gamma3) / gamma2,
    )
    # Compute coefficients rho appearing in arXiv:1508.07253 eq. (30), the numerical coefficients are in Tab. 5
    rho1 =
        3931.8979897196696 - 17395.758706812805 * eta +
        (
            3132.375545898835 + 343965.86092361377 * eta - 1.2162565819981997e6 * eta2 +
            (-70698.00600428853 + 1.383907177859705e6 * eta - 3.9662761890979446e6 * eta2) *
            xi +
            (-60017.52423652596 + 803515.1181825735 * eta - 2.091710365941658e6 * eta2) *
            xi2
        ) * xi
    rho2 =
        -40105.47653771657 +
        112253.0169706701 * eta +
        (
            23561.696065836168 - 3.476180699403351e6 * eta +
            1.137593670849482e7 * eta2 +
            (754313.1127166454 - 1.308476044625268e7 * eta + 3.6444584853928134e7 * eta2) *
            xi +
            (596226.612472288 - 7.4277901143564405e6 * eta + 1.8928977514040343e7 * eta2) *
            xi2
        ) * xi
    rho3 =
        83208.35471266537 - 191237.7264145924 * eta +
        (
            -210916.2454782992 + 8.71797508352568e6 * eta - 2.6914942420669552e7 * eta2 +
            (
                -1.9889806527362722e6 + 3.0888029960154563e7 * eta -
                8.390870279256162e7 * eta2
            ) * xi +
            (
                -1.4535031953446497e6 + 1.7063528990822166e7 * eta -
                4.2748659731120914e7 * eta2
            ) * xi2
        ) * xi
    # Compute coefficients delta appearing in arXiv:1508.07253 eq. (21)
    f1Interm = fInsJoin_Ampl
    f3Interm = fpeak
    dfInterm = 0.5 * (f3Interm - f1Interm)
    f2Interm = f1Interm + dfInterm
    # First write the inspiral coefficients, we put them in a structure and label with the power in front of which they appear
    amp0 = sqrt(2.0 * eta / 3.0) * (pi^(-1.0 / 6.0))
    Acoeffs = AcoeffsStructure(
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / 672.0,
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0,
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) * (pi^(4.0 / 3.0))
        ) / 8.128512e6,
        (
            (pi^(5.0 / 3.0)) * (
                chi2 * (
                    -285197.0 * (-1.0 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1.0 + 4.0 * eta) * pi
            )
        ) / 32256.0,
        -(
            (pi^2.0) * (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368.0 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320.0 * pi2
                )
            )
        ) / 6.0085960704e10,
        rho1,
        rho2,
        rho3,
    )

    # v1 is the inspiral model evaluated at f1Interm
    v1 =
        1.0 +
        (f1Interm^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (f1Interm^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (f1Interm^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (f1Interm^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (f1Interm^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        f1Interm *
        (Acoeffs.one + f1Interm * Acoeffs.two + f1Interm * f1Interm * Acoeffs.three)
    # d1 is the derivative of the inspiral model evaluated at f1
    d1 =
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / (1008.0 * (f1Interm^(1.0 / 3.0))) +
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0 +
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) *
            (f1Interm^(1.0 / 3.0)) *
            (pi^(4.0 / 3.0))
        ) / 6.096384e6 +
        (
            5.0 *
            (f1Interm^(2.0 / 3.0)) *
            (pi^(5.0 / 3.0)) *
            (
                chi2 * (
                    -285197.0 * (-1 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1 + 4 * eta) * pi
            )
        ) / 96768.0 -
        (
            f1Interm *
            pi2 *
            (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320 * pi2
                )
            )
        ) / 3.0042980352e10 +
        (7.0 / 3.0) * (f1Interm^(4.0 / 3.0)) * rho1 +
        (8.0 / 3.0) * (f1Interm^(5.0 / 3.0)) * rho2 +
        3.0 * (f1Interm * f1Interm) * rho3
    # v3 is the merger-ringdown model (eq. (19) of arXiv:1508.07253) evaluated at f3
    v3 =
        exp(-(f3Interm - fring) * gamma2 / (fdamp * gamma3)) * (fdamp * gamma3 * gamma1) /
        ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
    # d2 is the derivative of the merger-ringdown model evaluated at f3
    d2 =
        (
            (-2.0 * fdamp * (f3Interm - fring) * gamma3 * gamma1) /
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3) -
            (gamma2 * gamma1)
        ) / (
            exp((f3Interm - fring) * gamma2 / (fdamp * gamma3)) *
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
        )
    # v2 is the value of the amplitude evaluated at f2. They come from the fit of the collocation points in the intermediate region
    v2 =
        0.8149838730507785 +
        2.5747553517454658 * eta +
        (
            1.1610198035496786 - 2.3627771785551537 * eta +
            6.771038707057573 * eta2 +
            (0.7570782938606834 - 2.7256896890432474 * eta + 7.1140380397149965 * eta2) *
            xi +
            (0.1766934149293479 - 0.7978690983168183 * eta + 2.1162391502005153 * eta2) *
            xi2
        ) * xi
    # Now some definitions to speed up
    f1 = f1Interm
    f2 = f2Interm
    f3 = f3Interm
    f12 = f1Interm * f1Interm
    f13 = f1Interm * f12
    f14 = f1Interm * f13
    f15 = f1Interm * f14
    f22 = f2Interm * f2Interm
    f23 = f2Interm * f22
    f24 = f2Interm * f23
    f32 = f3Interm * f3Interm
    f33 = f3Interm * f32
    f34 = f3Interm * f33
    f35 = f3Interm * f34
    # Finally conpute the deltas
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2.0 * d1 * f13 * f22 * f33 - 2.0 * d2 * f13 * f22 * f33 -
            d1 * f12 * f23 * f33 + d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 -
            d1 * f13 * f2 * f34 - d1 * f12 * f22 * f34 +
            2.0 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4.0 * f12 * f23 * f32 * v1 -
            3.0 * f1 * f24 * f32 * v1 - 8.0 * f12 * f22 * f33 * v1 +
            4.0 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4.0 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2.0 * f23 * f34 * v1 - 2.0 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3.0 * f14 * f33 * v2 -
            3.0 * f13 * f34 * v2 + f12 * f35 * v2 - f15 * f22 * v3 +
            2.0 * f14 * f23 * v3 - f13 * f24 * v3 + 2.0 * f15 * f2 * f3 * v3 -
            f14 * f22 * f3 * v3 - 4.0 * f13 * f23 * f3 * v3 +
            3.0 * f12 * f24 * f3 * v3 - 4.0 * f14 * f2 * f32 * v3 +
            8.0 * f13 * f22 * f32 * v3 - 4.0 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2 * d1 * f13 * f22 * f33 - 2 * d2 * f13 * f22 * f33 - d1 * f12 * f23 * f33 +
            d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 - d1 * f13 * f2 * f34 -
            d1 * f12 * f22 * f34 +
            2 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4 * f12 * f23 * f32 * v1 -
            3 * f1 * f24 * f32 * v1 - 8 * f12 * f22 * f33 * v1 +
            4 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2 * f23 * f34 * v1 - 2 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3 * f14 * f33 * v2 - 3 * f13 * f34 * v2 +
            f12 * f35 * v2 - f15 * f22 * v3 + 2 * f14 * f23 * v3 - f13 * f24 * v3 +
            2 * f15 * f2 * f3 * v3 - f14 * f22 * f3 * v3 - 4 * f13 * f23 * f3 * v3 +
            3 * f12 * f24 * f3 * v3 - 4 * f14 * f2 * f32 * v3 +
            8 * f13 * f22 * f32 * v3 - 4 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta1 = -(
        (
            -(d2 * f15 * f22) + 2.0 * d2 * f14 * f23 - d2 * f13 * f24 -
            d2 * f14 * f22 * f3 +
            2.0 * d1 * f13 * f23 * f3 +
            2.0 * d2 * f13 * f23 * f3 - 2 * d1 * f12 * f24 * f3 - d2 * f12 * f24 * f3 +
            d2 * f15 * f32 - 3 * d1 * f13 * f22 * f32 - d2 * f13 * f22 * f32 +
            2 * d1 * f12 * f23 * f32 - 2 * d2 * f12 * f23 * f32 +
            d1 * f1 * f24 * f32 +
            2 * d2 * f1 * f24 * f32 - d2 * f14 * f33 +
            d1 * f12 * f22 * f33 +
            3 * d2 * f12 * f22 * f33 - 2 * d1 * f1 * f23 * f33 -
            2 * d2 * f1 * f23 * f33 +
            d1 * f24 * f33 +
            d1 * f13 * f34 +
            d1 * f1 * f22 * f34 - 2 * d1 * f23 * f34 - d1 * f12 * f35 + d1 * f22 * f35 -
            8 * f12 * f23 * f3 * v1 +
            6 * f1 * f24 * f3 * v1 +
            12 * f12 * f22 * f32 * v1 - 8 * f1 * f23 * f32 * v1 - 4 * f12 * f34 * v1 +
            2 * f1 * f35 * v1 +
            2 * f15 * f3 * v2 - 4 * f14 * f32 * v2 + 4 * f12 * f34 * v2 -
            2 * f1 * f35 * v2 - 2 * f15 * f3 * v3 + 8 * f12 * f23 * f3 * v3 -
            6 * f1 * f24 * f3 * v3 + 4 * f14 * f32 * v3 - 12 * f12 * f22 * f32 * v3 +
            8 * f1 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta2 = -(
        (
            d2 * f15 * f2 - d1 * f13 * f23 - 3 * d2 * f13 * f23 +
            d1 * f12 * f24 +
            2.0 * d2 * f12 * f24 - d2 * f15 * f3 + d2 * f14 * f2 * f3 -
            d1 * f12 * f23 * f3 +
            d2 * f12 * f23 * f3 +
            d1 * f1 * f24 * f3 - d2 * f1 * f24 * f3 - d2 * f14 * f32 +
            3 * d1 * f13 * f2 * f32 +
            d2 * f13 * f2 * f32 - d1 * f1 * f23 * f32 + d2 * f1 * f23 * f32 -
            2 * d1 * f24 * f32 - d2 * f24 * f32 - 2 * d1 * f13 * f33 +
            2 * d2 * f13 * f33 - d1 * f12 * f2 * f33 - 3 * d2 * f12 * f2 * f33 +
            3 * d1 * f23 * f33 +
            d2 * f23 * f33 +
            d1 * f12 * f34 - d1 * f1 * f2 * f34 + d1 * f1 * f35 - d1 * f2 * f35 +
            4 * f12 * f23 * v1 - 3 * f1 * f24 * v1 + 4 * f1 * f23 * f3 * v1 -
            3 * f24 * f3 * v1 - 12 * f12 * f2 * f32 * v1 +
            4 * f23 * f32 * v1 +
            8 * f12 * f33 * v1 - f1 * f34 * v1 - f35 * v1 - f15 * v2 - f14 * f3 * v2 +
            8 * f13 * f32 * v2 - 8 * f12 * f33 * v2 +
            f1 * f34 * v2 +
            f35 * v2 +
            f15 * v3 - 4 * f12 * f23 * v3 +
            3 * f1 * f24 * v3 +
            f14 * f3 * v3 - 4 * f1 * f23 * f3 * v3 + 3 * f24 * f3 * v3 -
            8 * f13 * f32 * v3 + 12 * f12 * f2 * f32 * v3 - 4 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta3 = -(
        (
            -2.0 * d2 * f14 * f2 + d1 * f13 * f22 + 3 * d2 * f13 * f22 - d1 * f1 * f24 -
            d2 * f1 * f24 + 2 * d2 * f14 * f3 - 2.0 * d1 * f13 * f2 * f3 -
            2 * d2 * f13 * f2 * f3 + d1 * f12 * f22 * f3 - d2 * f12 * f22 * f3 +
            d1 * f24 * f3 +
            d2 * f24 * f3 +
            d1 * f13 * f32 - d2 * f13 * f32 - 2 * d1 * f12 * f2 * f32 +
            2 * d2 * f12 * f2 * f32 +
            d1 * f1 * f22 * f32 - d2 * f1 * f22 * f32 + d1 * f12 * f33 -
            d2 * f12 * f33 +
            2 * d1 * f1 * f2 * f33 +
            2 * d2 * f1 * f2 * f33 - 3 * d1 * f22 * f33 - d2 * f22 * f33 -
            2 * d1 * f1 * f34 + 2 * d1 * f2 * f34 - 4 * f12 * f22 * v1 +
            2 * f24 * v1 +
            8 * f12 * f2 * f3 * v1 - 4 * f1 * f22 * f3 * v1 - 4 * f12 * f32 * v1 +
            8 * f1 * f2 * f32 * v1 - 4 * f22 * f32 * v1 - 4 * f1 * f33 * v1 +
            2 * f34 * v1 +
            2 * f14 * v2 - 4 * f13 * f3 * v2 + 4 * f1 * f33 * v2 - 2 * f34 * v2 -
            2 * f14 * v3 + 4 * f12 * f22 * v3 - 2 * f24 * v3 + 4 * f13 * f3 * v3 -
            8 * f12 * f2 * f3 * v3 +
            4 * f1 * f22 * f3 * v3 +
            4 * f12 * f32 * v3 - 8 * f1 * f2 * f32 * v3 + 4 * f22 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta4 = -(
        (
            d2 * f13 * f2 - d1 * f12 * f22 - 2 * d2 * f12 * f22 +
            d1 * f1 * f23 +
            d2 * f1 * f23 - d2 * f13 * f3 +
            2.0 * d1 * f12 * f2 * f3 +
            d2 * f12 * f2 * f3 - d1 * f1 * f22 * f3 + d2 * f1 * f22 * f3 -
            d1 * f23 * f3 - d2 * f23 * f3 - d1 * f12 * f32 + d2 * f12 * f32 -
            d1 * f1 * f2 * f32 - 2 * d2 * f1 * f2 * f32 +
            2 * d1 * f22 * f32 +
            d2 * f22 * f32 +
            d1 * f1 * f33 - d1 * f2 * f33 + 3 * f1 * f22 * v1 - 2 * f23 * v1 -
            6 * f1 * f2 * f3 * v1 +
            3 * f22 * f3 * v1 +
            3 * f1 * f32 * v1 - f33 * v1 - f13 * v2 + 3 * f12 * f3 * v2 -
            3 * f1 * f32 * v2 +
            f33 * v2 +
            f13 * v3 - 3 * f1 * f22 * v3 + 2 * f23 * v3 - 3 * f12 * f3 * v3 +
            6 * f1 * f2 * f3 * v3 - 3 * f22 * f3 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )

    # Defined as in LALSimulation - LALSimIMRPhenomD.c line 332. Final units are correctly Hz^-1
    Overallamp =
        2.0 * sqrt(5.0 / (64.0 * pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    amplitudeIMR = @. ifelse(
        fgrid < fInsJoin_Ampl,
        1.0 +
        (fgrid^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (fgrid^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (fgrid^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (fgrid^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (fgrid^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        fgrid * (Acoeffs.one + fgrid * Acoeffs.two + fgrid * fgrid * Acoeffs.three),
        ifelse(
            fgrid < fpeak,
            delta0 +
            fgrid * delta1 +
            fgrid * fgrid * (delta2 + fgrid * delta3 + fgrid * fgrid * delta4),
            ifelse(
                fgrid < fcutPar,
                exp(-(fgrid - fring) * gamma2 / (fdamp * gamma3)) *
                (fdamp * gamma3 * gamma1) /
                ((fgrid - fring) * (fgrid - fring) + fdamp * gamma3 * fdamp * gamma3),
                0.0,
            ),
        ),
    )
    if typeof(eta) == Float64
        ampl = Overallamp * amp0 .* (fgrid .^ (-7.0 / 6.0)) .* amplitudeIMR
        if container !== nothing
            container .= [ampl[i] for i in eachindex(ampl)]
         end
    else 
        ampl = Overallamp * amp0 .* (fgrid .^ (-7.0 / 6.0)) .* amplitudeIMR

        if container !== nothing
            container .= [ampl[i].value for i in eachindex(ampl)]
        end
    end

    return ampl
end

