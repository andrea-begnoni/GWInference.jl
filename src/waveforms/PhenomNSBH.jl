#! format: off

#######################################################
# IMRPhenomNSBH WAVEFORM
#######################################################

"""
ToDo: Need documentation
"""
function PolAbs(model::PhenomNSBH,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda,
    Lambda2=0.0; #Actually not required for this waveform
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc
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
        Lambda,
        fcutPar = fcutPar,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc
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
function Pol(model::PhenomNSBH,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda,
    Lambda2=0.0; #Actually not required for this waveform
    fcutPar = 0.2,
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
        Lambda,
        fcutPar = fcutPar,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

"""
IMRPhenomNSBH waveform model

The inputs labelled as 1 refer to the BH 
The inputs labelled as 2 refer to the NS
There is only one Lambda, the one of the NS
Ref: arXiv:1508.07250, arXiv:1508.07253, arXiv:1509.00512, arXiv:1905.06011
"""
    # Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase
    # Dimensionless frequency (Mf) at which we define the end of the waveform
    

"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomNSBH(), f, mc, eta, chi1, chi2, Lambda)

Note that the second spin is assumed to be zero.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `Lambda`: Dimensionless tidal deformability of the NS.
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
    Phi(PhenomNSBH(), f, mc, eta, chi1, chi2)
```
"""
function Phi(model::PhenomNSBH,
        f,
        mc,
        eta,
        chi1,
        chi2, # tuned only for chi2 = 0
        Lambda;
        fInsJoin = 0.018,
        fcutPar = 0.2,
        GMsun_over_c3 = uc.GMsun_over_c3,
        interpolation = false,
    )


    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    eta3 = eta * eta2
    etaInv = 1. / eta

    pi2 = pi * pi
    
    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2 = chi1*chi2

    # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
    # Taken from arXiv:1303.1528 eq. (54) and Tab. I
    QuadMon1 =1.
    QuadMon2 = ifelse(Lambda<1e-5, 1., exp(0.194 + 0.0936*log(Lambda) + 0.0474*log(Lambda)^2 - 0.00421*log(Lambda)^3 + 0.000123*log(Lambda)^4))
    
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)        
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a
    chi1dotchi2 = chi1*chi2
    chi_sdotchi_a = chi_s*chi_a
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M*GMsun_over_c3.*f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = - 1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy for IMRPhenomNSBH, the rest is equivalent to IMRPhenomD_NRTidalv2
    # Get remnant spin for assumed aligned spin system, from arXiv:1903.11622 Table I and eq. (4), (5) and (6)
    
    p1_remSp = ((-5.44187381e-03*chi1 + 7.91165608e-03) + (2.33362046e-02*chi1 + 2.47764497e-02)*eta)*eta
    p2_remSp = ((-8.56844797e-07*chi1 - 2.81727682e-06) + (6.61290966e-06*chi1 + 4.28979016e-05)*eta)*eta
    p3_remSp = ((-3.04174272e-02*chi1 + 2.54889050e-01) + (1.47549350e-01*chi1 - 4.27905832e-01)*eta)*eta
    
    modelRemSp = (1. + Lambda * p1_remSp + Lambda^2 * p2_remSp) / ((1. + Lambda*p3_remSp^2)*(1. + Lambda*p3_remSp^2))

    modelRemSp = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemSp)
    modelRemSp = ifelse(chi1 < -0.5, 1., modelRemSp)
    modelRemSp = ifelse(modelRemSp > 1., 1., modelRemSp)
    
    #del p1_remSp, p2_remSp, p3_remSp
    
    # Work with spin variables weighted on square of the BH mass over total mass
    S1BH = chi1 * m1ByM^2
    Shat = S1BH / (m1ByM^2 + m2ByM^2) # this would be = (chi1*m1*m1 + chi2*m2*m2)/(m1*m1 + m2*m2), but chi2=0 by assumption
    
    # Compute fit to L_orb in arXiv:1611.00332 eq. (16)
    Lorb = (2. *sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta3)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16. *4.409160174224525 - 4. *0.5118334706832706)*eta3) 
    + 0.0851*0.7877509372255369*Shat^2*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16. *8.77367320110712 - 4. *(-32.060648277652994))*eta3) + 0.00954*0.6540138407185817*Shat^3*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16. *22.830033250479833 - 4. *(-153.83722669033995))*eta3))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0. *eta2 + (64. - 64. *1.8804718791591157 - 16. *(-4.770246856212403) - 4. *0.)*eta3)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta3*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta3*chi12;
    
    chif = (Lorb + S1BH)*modelRemSp

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
        fring = linear_interpolation(QNMgrid_a, QNMgrid_fring)(chif) / (1.0 - Erad)
        fdamp = linear_interpolation(QNMgrid_a, QNMgrid_fdamp)(chif) / (1.0 - Erad)
    else
        fring = ((0.05947169566573468 - 0.14989771215394762*chif + 0.09535606290986028*chif*chif + 0.02260924869042963*chif*chif*chif - 0.02501704155363241*chif*chif*chif*chif - 0.005852438240997211*(chif^5) + 0.0027489038393367993*(chif^6) + 0.0005821983163192694*(chif^7))/(1 - 2.8570126619966296*chif + 2.373335413978394*chif*chif - 0.6036964688511505*chif*chif*chif*chif + 0.0873798215084077*(chif^6)))/(1. - Erad)
        fdamp = ((0.014158792290965177 - 0.036989395871554566*chif + 0.026822526296575368*chif*chif + 0.0008490933750566702*chif*chif*chif - 0.004843996907020524*chif*chif*chif*chif - 0.00014745235759327472*(chif^5) + 0.0001504546201236794*(chif^6))/(1 - 2.5900842798681376*chif + 1.8952576220623967*chif*chif - 0.31416610693042507*chif*chif*chif*chif + 0.009002719412204133*(chif^6)))/(1. - Erad)
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
        m1ByM^2*
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
    
    fMRDJoin = 0.5*fring
    
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

    fpeak = fgrid[end] # In LAL the maximum of the grid is used to rescale
    
    t0 = (alpha1 + alpha2/(fpeak*fpeak) + alpha3/(fpeak^(1. /4.)) + alpha4/(fdamp*(1. + (fpeak - alpha5*fring)*(fpeak - alpha5*fring)/(fdamp*fdamp))))/eta
    
    # LAL sets fRef as the minimum frequency, do the same
    fRef   = fgrid[1]

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

    kappa2T = (3.0/13.0) * ((1.0 + 12.0*m1ByM/m2ByM)*(m2ByM^5)*Lambda)
    
    c_Newt   = 2.4375
    n_1      = -12.615214237993088
    n_3over2 =  19.0537346970349
    n_2      = -21.166863146081035
    n_5over2 =  90.55082156324926
    n_3      = -60.25357801943598
    d_1      = -15.11120782773667
    d_3over2 =  22.195327350624694
    d_2      =   8.064109635305156

    numTidal = @. 1.0 + (n_1 * ((pi*fgrid)^(2. /3.))) + (n_3over2 * pi*fgrid) + (n_2 * ((pi*fgrid)^(4. /3.))) + (n_5over2 * ((pi*fgrid)^(5. /3.))) + (n_3 * (pi*fgrid)^2)
    denTidal = @. 1.0 + (d_1 * ((pi*fgrid)^(2. /3.))) + (d_3over2 * pi*fgrid) + (d_2 * ((pi*fgrid)^(4. /3.)))
    
    tidal_phase = @. - kappa2T * c_Newt / (m1ByM * m2ByM) * ((pi*fgrid)^(5. /3.)) * numTidal / denTidal
    # This pi factor is needed to include LAL fRef rescaling, so to end up with the exact same waveform
    return @. phis + ifelse(fgrid < fcutPar, - t0*(fgrid - fRef) - phiRef + pi +  tidal_phase, 0.)
end



"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomNSBH(), f, mc, eta,  chi1, chi2, dL, Lambda)
Note that the second spin is assumed to be zero.
#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the binary, in Gpc.
-  `Lambda`: Dimensionless tidal deformability of the NS.
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
    fcut = _fcut(PhenomNSBH(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomNSBH(), f, mc, eta, chi1, chi2, dL, Lambda)
```
"""
function Ampl(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, #assumed to be chi2=0
    dL,
    Lambda;
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # Useful quantities

    M = mc / (eta^(3.0 / 5.0))
    eta2 = eta * eta # This can speed up a bit, we call it multiple times
    eta3 = eta * eta2
    pi2 = pi * pi
    chi12, chi22 = chi1*chi1, chi2*chi2
    
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # As in arXiv:0909.2867
    chieff = m1ByM * chi1 + m2ByM * chi2
    chisum = 2. *chieff
    chiprod = chieff^2
    
    # compute needed IMRPhenomC attributes
    # First the SPA part, LALSimIMRPhenomC_internals.c line 38
    # Frequency-domain Amplitude coefficients
    xdotaN = 64. *eta/5.
    xdota2 = -7.43/3.36 - 11. *eta/4.
    xdota3 = 4. *pi - 11.3*chieff/1.2 + 19. *eta*chisum/6.
    xdota4 = 3.4103/1.8144 + 5*chiprod + eta*(13.661/2.016 - chiprod/8.) + 5.9*eta2/1.8
    xdota5 = -pi*(41.59/6.72 + 189. *eta/8.) - chieff*(31.571/1.008 - 116.5*eta/2.4) + chisum*(21.863*eta/1.008 - 79. *eta2/6.) - 3*chieff*chiprod/4. + 9. *eta*chieff*chiprod/4.
    xdota6 = 164.47322263/1.39708800 - 17.12*MathConstants.eulergamma/1.05 + 16. *pi*pi/3 - 8.56*log(16.)/1.05 + eta*(45.1*pi*pi/4.8 - 561.98689/2.17728) + 5.41*eta2/8.96 - 5.605*eta3/2.592 - 80. *pi*chieff/3. + eta*chisum*(20. *pi/3. - 113.5*chieff/3.6) + chiprod*(64.153/1.008 - 45.7*eta/3.6) - chiprod*(7.87*eta/1.44 - 30.37*eta2/1.44)
    xdota6log = -856. /105.
    xdota7 = -pi*(4.415/4.032 - 358.675*eta/6.048 - 91.495*eta2/1.512) - chieff*(252.9407/2.7216 - 845.827*eta/6.048 + 415.51*eta2/8.64) + chisum*(158.0239*eta/5.4432 - 451.597*eta2/6.048 + 20.45*eta3/4.32 + 107. *eta*chiprod/6. - 5. *eta2*chiprod/24.) + 12. *pi*chiprod - chiprod*chieff*(150.5/2.4 + eta/8.) + chieff*chiprod*(10.1*eta/2.4 + 3. *eta2/8.)
    # Time-domain amplitude coefficients, which also enters the fourier amplitude in this model
    AN = 8. *eta*sqrt(pi/5.)
    A2 = (-107. + 55. *eta)/42.
    A3 = 2. *pi - 4. *chieff/3. + 2. *eta*chisum/3.
    A4 = -2.173/1.512 - eta*(10.69/2.16 - 2. *chiprod) + 2.047*eta2/1.512
    A5 = -10.7*pi/2.1 + eta*(3.4*pi/2.1)
    A5imag = -24. *eta
    A6 = 270.27409/6.46800 - 8.56*MathConstants.eulergamma/1.05 + 2. *pi*pi/3. + eta*(4.1*pi*pi/9.6 - 27.8185/3.3264) - 20.261*eta2/2.772 + 11.4635*eta3/9.9792 - 4.28*log(16.)/1.05
    A6log = -428. /105.
    A6imag = 4.28*pi/1.05
    
    z701, z702, z711, z710, z720 = 4.149e+00, -4.070e+00, -8.752e+01, -4.897e+01, 6.665e+02
    z801, z802, z811, z810, z820 = -5.472e-02, 2.094e-02, 3.554e-01, 1.151e-01, 9.640e-01
    z901, z902, z911, z910, z920 = -1.235e+00, 3.423e-01, 6.062e+00, 5.949e+00, -1.069e+01
    
    g1 = z701 * chieff + z702 * chiprod + z711 * eta * chieff + z710 * eta + z720 * eta2
    g1 = ifelse(g1 < 0., 0., g1)
    
    del1 = z801 * chieff + z802 * chiprod + z811 * eta * chieff + z810 * eta + z820 * eta2
    del2 = z901 * chieff + z902 * chiprod + z911 * eta * chieff + z910 * eta + z920 * eta2
    del1 = ifelse(del1 < 0., 0., del1)
    del2 = ifelse(del2 < 1.0e-4, 1.0e-4, del2)
    
    d0 = 0.015
    
    # All the other coefficients from IMRPhenomC are not needed
    
    # Now compute NSBH coefficients
    # Get NS compactness and baryonic mass, see arXiv:1608.02582 eq. (78)
    a0Comp = 0.360
    a1Comp = -0.0355
    a2Comp = 0.000705
    
    comp = ifelse(Lambda > 1., a0Comp + a1Comp*log(Lambda) + a2Comp*log(Lambda)^2, 0.5 + (3. *a0Comp-a1Comp-1.5)*Lambda^2 + (-2. *a0Comp+a1Comp+1.)*Lambda^3)
    
    # Get baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass,
    # see arXiv:1509.00512 eq. (11)
    alphaTor = 0.296
    betaTor = 0.171
    
    if typeof(eta) == Float64
        xiTide = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2, 20., Order1(),[q, q*comp, chi1], maxiters=100)^2
    else    # else we are doing the Fisher matrix, i.e. we are using FowardDiff.Duals
        
        # if we are doing the Fisher matrix we need to derive the roots of the following equation w.r.t. q, q*comp, chi1
        # to do that we need to re-pack the Duals in a way that ForwardDiff can handle
        tmp = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2
                        ,20., Order1(),[q, q*comp, chi1], maxiters=100)^2
        partials = [(tmp.partials)[i].value for i in eachindex(tmp.partials)]
        xiTide = ForwardDiff.Dual{typeof(chi1).parameters[1]}((tmp.value).value, partials...)



    end

    Z1_ISCO = 1.0 + ((1.0 - chi1^2)^(1. /3.))*((1.0+chi1)^(1. /3.) + (1.0-chi1)^(1. /3.))
    Z2_ISCO = sqrt(3.0*chi1^2 + Z1_ISCO^2)
    r_ISCO  = ifelse(chi1>0., 3.0 + Z2_ISCO - sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)), 3.0 + Z2_ISCO + sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)))
    
    tmpMtorus = alphaTor * xiTide * (1.0-2.0*comp) - betaTor * q*comp * r_ISCO
    
    Mtorus = ifelse(tmpMtorus>0., tmpMtorus, 0.)
    
    
    # Get remnant spin for assumed aligned spin system, from arXiv:1903.11622 Table I and eq. (4), (5) and (6)
    
    p1_remSp = ((-5.44187381e-03*chi1 + 7.91165608e-03) + (2.33362046e-02*chi1 + 2.47764497e-02)*eta)*eta
    p2_remSp = ((-8.56844797e-07*chi1 - 2.81727682e-06) + (6.61290966e-06*chi1 + 4.28979016e-05)*eta)*eta
    p3_remSp = ((-3.04174272e-02*chi1 + 2.54889050e-01) + (1.47549350e-01*chi1 - 4.27905832e-01)*eta)*eta
    
    modelRemSp = (1. + Lambda * p1_remSp + Lambda^2 * p2_remSp) / ((1. + Lambda*p3_remSp^2)*(1. + Lambda*p3_remSp^3))

    modelRemSp = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemSp)
    modelRemSp = ifelse(chi1 < -0.5, 1., modelRemSp)
    modelRemSp = ifelse(modelRemSp > 1., 1., modelRemSp)
        
    # Work with spin variables weighted on square of the BH mass over total mass
    S1BH = chi1 * m1ByM^2 
    Shat = S1BH / (m1ByM^2 + m2ByM^2) # this would be = (chi1*m1*m1 + chi2*m2*m2)/(m1*m1 + m2*m2), but chi2=0 by assumption # ANDREA what? chi2 was not put to zero
    
    # Compute fit to L_orb in arXiv:1611.00332 eq. (16)
    
    Lorb = (2. *sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta3)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16. *4.409160174224525 - 4. *0.5118334706832706)*eta3) + 0.0851*0.7877509372255369*Shat^2*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16. *8.77367320110712 - 4. *(-32.060648277652994))*eta3) + 0.00954*0.6540138407185817*Shat^3*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16. *22.830033250479833 - 4. *(-153.83722669033995))*eta3))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0. *eta2 + (64. - 64. *1.8804718791591157 - 16. *(-4.770246856212403) - 4. *0.)*eta3)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta3*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta3*chi12;
    
    
    chif = (Lorb + S1BH)*modelRemSp
    
    # Get remnant mass scaled to a total (initial) mass of 1
    
    p1_remM = ((-1.83417425e-03*chi1 + 2.39226041e-03) + (4.29407902e-03*chi1 + 9.79775571e-03)*eta)*eta
    p2_remM = ((2.33868869e-07*chi1 - 8.28090025e-07) + (-1.64315549e-06*chi1 + 8.08340931e-06)*eta)*eta
    p3_remM = ((-2.00726981e-02*chi1 + 1.31986011e-01) + (6.50754064e-02*chi1 - 1.42749961e-01)*eta)*eta

    modelRemM = (1. + Lambda * p1_remM + Lambda^2 * p2_remM) / ((1. + Lambda*p3_remM^2)*(1. + Lambda*p3_remM^2))
    modelRemM = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemM)
    modelRemM = ifelse(chi1 < -0.5, 1., modelRemM)
    modelRemM = ifelse(modelRemM > 1., 1., modelRemM)
    
    # Compute the radiated-energy fit from arXiv:1611.00332 eq. (27)
    #ANDREA why Erad in Phi is different?
    EradNSBH = (((1. + -2.0/3.0*sqrt(2.))*eta + 0.5609904135313374*eta2 + (-0.84667563764404)*eta3 + 3.145145224278187*eta2*eta2)*(1. + 0.346*(-0.2091189048177395)*Shat*(1.8083565298668276 + 15.738082204419655*eta + (16. - 16. *1.8083565298668276 - 4. *15.738082204419655)*eta2) + 0.211*(-0.19709136361080587)*Shat^2*(4.271313308472851 + 0. *eta + (16. - 16. *4.271313308472851 - 4. *0.)*eta2) + 0.128*(-0.1588185739358418)*Shat^3*(31.08987570280556 + (-243.6299258830685)*eta + (16. - 16. *31.08987570280556 - 4. *(-243.6299258830685))*eta2)))/(1. + (-0.212)*2.9852925538232014*Shat*(1.5673498395263061 + (-0.5808669012986468)*eta + (16. - 16. *1.5673498395263061 - 4. *(-0.5808669012986468))*eta2)) + (-0.09803730445895877)*Seta*eta2*(1. + (-3.2283713377939134)*eta)*chi1 + (-0.01978238971523653)*Shat*Seta*eta*(1. + (-4.91667749015812)*eta)*chi1 + 0.01118530335431078*eta3*chi12
    finalMass = (1. -EradNSBH)*modelRemM
    
    # Compute 22 quasi-normal mode dimensionless frequency
    kappaOm = sqrt(log(2. -chif)/log(3.))
    omega_tilde = (1.0 + kappaOm*(1.5578*exp(1im*2.9031) + 1.9510*exp(1im*5.9210)*kappaOm + 2.0997*exp(1im*2.7606)*kappaOm^2 + 1.4109*exp(1im*5.9143)*kappaOm^3 + 0.4106*exp(1im*2.7952)*(kappaOm^4)))
    
    fring = 0.5*real(omega_tilde)/pi/finalMass
    
    rtide = xiTide * (1.0 - 2.0 * comp) / (q*comp)
    
    q_factor = 0.5*real(omega_tilde)/imag(omega_tilde)
    
    ftide = abs(1.0/(pi*(chi1 + sqrt(rtide^3)))*(1.0 + 1.0 / q))
    
    # Now compute last amplitude quantities
    fring_tilde = 0.99 * 0.98 * fring
    
    gamma_correction = ifelse(Lambda > 1.0, 1.25, 1.0 + 0.5*Lambda - 0.25*Lambda^2)
    delta_2_prime = ifelse(Lambda > 1.0, 1.62496*0.25*(1. + tanh(4.0*((ftide/fring_tilde - 1.)-0.0188092)/0.338737)), del2 - 2. *(del2 - 0.81248)*Lambda + (del2 - 0.81248)*Lambda^2)
    
    sigma = delta_2_prime * fring / q_factor
    
    # Determine the type of merger we see and determine coefficients
    epsilon_tide = ifelse(ftide < fring, 0., 2. *0.25*(1 + tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.571505*comp - 0.00508451*chi1)+0.0796251)/0.0801192)))
    
    epsilon_ins  = ifelse(ftide < fring, ifelse(1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1)>1., 1., 1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1)), ifelse(Mtorus > 0., 1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1), 1.))
    
    sigma_tide   = ifelse(ftide < fring, ifelse(Mtorus>0., 0.137722 - 0.293237*(Mtorus - 0.132754*comp + 0.576669*sqrt(eta) - 0.0603749*chi1 - 0.0601185*chi1^2 - 0.0729134*chi1^3), 0.5*(0.137722 - 0.293237*(Mtorus - 0.132754*comp + 0.576669*sqrt(eta) - 0.0603749*chi1 - 0.0601185*chi1^2 - 0.0729134*chi1^3) + 0.5*(1. - tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.657424*comp - 0.0259977*chi1)+0.206465)/0.226844)))),0.5*(1. - tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.657424*comp - 0.0259977*chi1)+0.206465)/0.226844)))
    
    f0_tilde_PN  = ifelse(ftide < fring, ifelse(Mtorus>0., ftide / (M*GMsun_over_c3), ((1.0 - 1.0 / q) * fring_tilde + epsilon_ins * ftide / q)/(M*GMsun_over_c3)), ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))
    
    f0_tilde_PM  = ifelse(ftide < fring, ifelse(Mtorus>0., ftide / (M*GMsun_over_c3), ((1.0 - 1.0 / q) * fring_tilde + ftide/q)/(M*GMsun_over_c3)), ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))
    
    f0_tilde_RD  = ifelse(ftide < fring, 0., ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))

    # This can be used to output the merger type if needed
    
    v = @. (fgrid*pi)^(1. /3.)
    v2 = v .* v
    xdot = @. xdotaN*(v^10)*(1. + xdota2*v2 + xdota3 * fgrid*pi + xdota4 * fgrid*pi*v + xdota5 * v2*fgrid*pi + (xdota6 + xdota6log * 2. *log(v)) * (fgrid*pi)^2 + xdota7 * v*(fgrid*pi)^2)
    ampfacTime = @. sqrt(abs(pi / (1.5 * v * xdot)))
    
    AmpPNre = @. ampfacTime * AN * v2 * (1. + A2*v2 + A3 * fgrid*pi + A4 * v*fgrid*pi + A5 * v2*fgrid*pi + (A6 + A6log * 2. *log(v)) * (fgrid*pi)^2)
    AmpPNim = @. ampfacTime * AN * v2 * (A5imag * v2*fgrid*pi + A6imag * (fgrid*pi)^2)
    
    aPN = @. sqrt(AmpPNre * AmpPNre + AmpPNim * AmpPNim)
    aPM = @. (gamma_correction * g1 * (fgrid^(5. /6.)))

    LRD = @. sigma^2 / ((fgrid - fring) * (fgrid - fring) + sigma^2*0.25)
    aRD = @. epsilon_tide * del1 * LRD * (fgrid^(-7. /6.))
    
    wMinusf0_PN = @. 0.5 * (1. - tanh(4. *(fgrid - (epsilon_ins * f0_tilde_PN)*M*GMsun_over_c3)/(d0 + sigma_tide)))
    wMinusf0_PM = @. 0.5 * (1. - tanh(4. *(fgrid - f0_tilde_PM*M*GMsun_over_c3)/(d0 + sigma_tide)))
    wPlusf0     = @. 0.5 * (1. + tanh(4. *(fgrid - f0_tilde_RD*M*GMsun_over_c3)/(d0 + sigma_tide)))
    
    amplitudeIMR = @. ifelse(fgrid < fcutPar, (aPN * wMinusf0_PN + aPM * wMinusf0_PM + aRD * wPlusf0), 0.)
    
    # Defined as in LALSimulation - LALSimIMRPhenomD.c line 332. Final units are correctly Hz^-1
    Overallamp = 2. * sqrt(5. /(64. *pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL
    return Overallamp.*amplitudeIMR


end
