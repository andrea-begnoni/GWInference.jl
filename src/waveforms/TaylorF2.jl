#! format: off

####################################################
# TAYLORF2 WAVEFORM
####################################################

"""
ToDo: Need documentation
"""
function PolAbs(model::TaylorF2,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0;
    clightGpc = uc.clightGpc, 
    GMsun_over_c3 = uc.GMsun_over_c3
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        mc,
        dL,
        clightGpc = clightGpc, 
        GMsun_over_c3 = GMsun_over_c3
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
function Pol(model::TaylorF2,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0;
    clightGpc = uc.clightGpc, 
    GMsun_over_c3 = uc.GMsun_over_c3
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
        clightGpc = clightGpc, 
        GMsun_over_c3 = GMsun_over_c3
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(TaylorF2(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    f = 10 .^range(1, stop=3, length=100)
    Phi(TaylorF2(), f, 30., 0.25, 0.5, 0.5)
```
"""
function Phi(model::TaylorF2,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1=0.,
    Lambda2=0.;
    GMsun_over_c3 = uc.GMsun_over_c3,
    is_tidal = false,
    use_QuadMonTid = false,
    is_eccentric = false,
    phiref_vlso = false,
    use_3p5PN_SpinHO = false,
    fRef_ecc = nothing,
)

    # From A. Buonanno, B. Iyer, E. Ochsner, Y. Pan, B.S. Sathyaprakash - arXiv:0907.0700 - eq. (3.18) plus spins as in arXiv:1107.1267 eq. (5.3) up to 2.5PN and PhysRevD.93.084054 eq. (6) for 3PN and 3.5PN
    Mtot_sec = mc * GMsun_over_c3/eta^(0.6)
    v = @. (pi*Mtot_sec*f)^(1. /3.)
    eta2 = eta*eta

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)        
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    
    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2  = chi1*chi2
    chi_s, chi_a   = 0.5*(chi1 + chi2), 0.5*(chi1 - chi2)
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a

    vlso = 1. /sqrt(6.)
    if Lambda1 > 0. || Lambda2 > 0.
        is_tidal = true
        use_QuadMonTid = true
    end
    
    
    if is_tidal && use_QuadMonTid
        # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
        # The relation between the two is given in arxiv:1608.02582 eq. (15) with coefficients from third row of Table I
        # We also extend the range to 0 <= Lam < 1, as done in LALSimulation in LALSimUniversalRelations.c line 123
        QuadMon1 = ifelse(Lambda1 < 1., 1. + Lambda1*(0.427688866723244 + Lambda1*(-0.324336526985068 + Lambda1*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda1) + 0.04812 * log(Lambda1)^2 -4.283e-3 * log(Lambda1)^3 + 1.245e-4 * log(Lambda1)^4))
        QuadMon2 = ifelse(Lambda2 < 1., 1. + Lambda2*(0.427688866723244 + Lambda2*(-0.324336526985068 + Lambda2*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda2) + 0.04812 * log(Lambda2)^2 -4.283e-3 * log(Lambda2)^3 + 1.245e-4 * log(Lambda2)^4))
    else
        QuadMon1, QuadMon2 = 1.,1.
    end
    TF2OverallAmpl = 3. /(128. * eta)
    TF2_5coeff_tmp = (732985. /2268. - 24260. *eta/81. - 340. *eta2/9.)*chi_s + (732985. /2268. + 140. *eta/9.)*Seta*chi_a

    if phiref_vlso
        TF2coeffs_five = (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)*(1. -3. *log(vlso))
        phiR = 0.
    else
        TF2coeffs_five = (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)
        # This pi factor is needed to include LAL fRef rescaling, so to end up with the exact same waveform
        phiR = pi
    end

    if use_3p5PN_SpinHO
        # This part includes SS and SSS contributions at 3.5PN, which are not included in LAL
        TF2coeffs_seven = 77096675. *pi/254016. + 378515. *pi*eta/1512. - 74045. *pi*eta2/756. + (-25150083775. /3048192. + 10566655595. *eta/762048. - 1042165. *eta2/3024. + 5345. *eta^3/36. + (14585. /8. - 7270. *eta + 80. *eta2)*chi_a2)*chi_s + (14585. /24. - 475. *eta/6. + 100. *eta2/3.)*chi_s2*chi_s + Seta*((-25150083775. /3048192. + 26804935. *eta/6048. - 1985. *eta2/48.)*chi_a + (14585. /24. - 2380. *eta)*chi_a2*chi_a + (14585. /8. - 215. *eta/2.)*chi_a*chi_s2)
    else
        TF2coeffs_seven = 77096675. *pi/254016. + 378515. *pi*eta/1512. - 74045. *pi*eta2/756. + (-25150083775. /3048192. + 10566655595. *eta/762048. - 1042165. *eta2/3024. + 5345. *eta^3/36.)*chi_s + Seta*((-25150083775. /3048192. + 26804935. *eta/6048. - 1985. *eta2/48.)*chi_a)
    end

    TF2coeffs = TF2coeffsStructure(
        1.,
        0.,
        3715. /756. + (55. *eta)/9.,
        -16. *pi + (113. *Seta*chi_a)/3. + (113. /3. - (76. *eta)/3.)*chi_s,
        5. *(3058.673/7.056 + 5429. /7. *eta+617. *eta2)/72. + 247. /4.8*eta*chi1dotchi2 -721. /4.8*eta*chi1dotchi2 + (-720. /9.6*QuadMon1 + 1. /9.6)*m1ByM^2*chi12 + (-720. /9.6*QuadMon2 + 1. /9.6)*m2ByM^2*chi22 + (240. /9.6*QuadMon1 - 7. /9.6)*m1ByM^2*chi12 + (240. /9.6*QuadMon2 - 7. /9.6)*m2ByM^2*chi22,
        TF2coeffs_five,
        (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)*3,
        11583.231236531/4.694215680 - 640. /3. *pi^2 - 684.8/2.1* MathConstants.eulergamma + eta*(-15737.765635/3.048192 + 225.5/1.2*pi^2) + eta2*76.055/1.728 - eta^3*127.825/1.296 - log(4.)*684.8/2.1 + pi*chi1*m1ByM*(1490. /3. + m1ByM*260.) + pi*chi2*m2ByM*(1490. /3. + m2ByM*260.) + (326.75/1.12 + 557.5/1.8*eta)*eta*chi1dotchi2 + (4703.5/8.4+2935. /6. *m1ByM-120. *m1ByM^2)*m1ByM^2*QuadMon1*chi12 + (-4108.25/6.72-108.5/1.2*m1ByM+125.5/3.6*m1ByM^2)*m1ByM^2*chi12 + (4703.5/8.4+2935. /6. *m2ByM-120. *m2ByM^2)*m2ByM^2*QuadMon2*chi22 + (-4108.25/6.72-108.5/1.2*m2ByM+125.5/3.6*m2ByM^2)*m2ByM^2*chi22,
        -(6848. /21.),
        TF2coeffs_seven,
    )
    
    if is_eccentric
        # These are the eccentricity dependent coefficients up to 3 PN order, in the low-eccentricity limit, from arXiv:1605.00304
        
        if isnothing(fRef_ecc)
            v0ecc = amin(v)
        else
            v0ecc = (pi*Mtot_sec*fRef_ecc)^(1. /3.)
        end  
        TF2EccCoeffs = TF2EccCoeffsStruct(
            1.,
            0.,
            29.9076223/8.1976608 + 18.766963/2.927736*eta,
            2.833/1.008 - 19.7/3.6*eta,
            -28.19123/2.82600*pi,
            37.7/7.2*pi,
            16.237683263/3.330429696 + 241.33060753/9.71375328*eta+156.2608261/6.9383952*eta2,
            84.7282939759/8.2632420864-7.18901219/3.68894736*eta-36.97091711/1.05398496*eta2,
            -1.193251/3.048192 - 66.317/9.072*eta +18.155/1.296*eta2,
            -28.31492681/1.18395270*pi - 115.52066831/2.70617760*pi*eta,
            -79.86575459/2.84860800*pi + 55.5367231/1.0173600*pi*eta,
            112.751736071/5.902315776*pi + 70.75145051/2.10796992*pi*eta,
            76.4881/9.0720*pi - 94.9457/2.2680*pi*eta,
            -436.03153867072577087/1.32658535116800000 + 53.6803271/1.9782000*MathConstants.eulergamma + 157.22503703/3.25555200*pi^2 +(2991.72861614477/6.89135247360 - 15.075413/1.446912*pi^2)*eta +345.5209264991/4.1019955200*eta2 + 506.12671711/8.78999040*eta^3 + 384.3505163/5.9346000*log(2.) - 112.1397129/1.7584000*log(3.),
            46.001356684079/3.357073133568 + 253.471410141755/5.874877983744*eta - 169.3852244423/2.3313007872*eta2 - 307.833827417/2.497822272*eta^3,
            -106.2809371/2.0347200*pi^2,
            -3.56873002170973/2.49880440692736 - 260.399751935005/8.924301453312*eta + 15.0484695827/3.5413894656*eta2 + 340.714213265/3.794345856*eta^3,
            265.31900578691/1.68991764480 - 33.17/1.26*MathConstants.eulergamma + 12.2833/1.0368*pi^2 + (91.55185261/5.48674560 - 3.977/1.152*pi^2)*eta - 5.732473/1.306368*eta2 - 30.90307/1.39968*eta^3 + 87.419/1.890*log(2.) - 260.01/5.60*log(3.),
        )
        
        TF2EccOverallAmpl = -2.355/1.462*ecc^2*((v0ecc/v)^(19. /3.))
        phi_Ecc = TF2EccOverallAmpl * (TF2EccCoeffs.zero + TF2EccCoeffs.one*v + (TF2EccCoeffs.twoV*v^2 + TF2EccCoeffs.twoV0*v0ecc^2) + (TF2EccCoeffs.threeV*v^3 + TF2EccCoeffs.threeV0*v0ecc^3) + (TF2EccCoeffs.fourV4*v^4 + TF2EccCoeffs.fourV2V02*v^2*v0ecc^2 + TF2EccCoeffs.fourV04*v0ecc^4) + (TF2EccCoeffs.fiveV5*v^5 + TF2EccCoeffs.fiveV3V02*v^3*v0ecc^2 + TF2EccCoeffs.fiveV2V03*v^2*v0ecc^3 + TF2EccCoeffs.fiveV05*v0ecc^5) + ((TF2EccCoeffs.sixV6 + 53.6803271/3.9564000*log(16. *v^2))*(v^6)+ TF2EccCoeffs.sixV4V02*v^4*v0ecc^2 + TF2EccCoeffs.sixV3V03*v^3*v0ecc^3 + TF2EccCoeffs.sixV2V04*v^2*v0ecc^4 + (TF2EccCoeffs.sixV06 - 33.17/2.52*log(16. *v0ecc^2))*v0ecc^6))
    else
        phi_Ecc = 0.
    end    
    if is_tidal
        # Add tidal contribution if needed, as in PhysRevD.89.103012
        Lam_t, delLam    = uc.Lamt_delLam_from_Lam12(eta, Lambda1, Lambda2)
        phi_Tidal = @. (-0.5*39. *Lam_t)*(v^10.) + (-3115. /64. *Lam_t + 6595. /364. *Seta*delLam)*(v^12.)
    
    else
        phi_Tidal = 0.
    end
    phase = @. TF2OverallAmpl*(TF2coeffs.zero + TF2coeffs.one*v + TF2coeffs.two*v^2 + TF2coeffs.three*v^3 + TF2coeffs.four*v^4 + (TF2coeffs.five + TF2coeffs.five_log*log(v))*v^5 + (TF2coeffs.six + TF2coeffs.six_log*log(v))*v^6 + TF2coeffs.seven*v^7 + phi_Tidal + phi_Ecc)/(v^5.)
    
    return phase .+ phiR .- pi*0.25
end


"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(TaylorF2(), f, mc, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `dL`: Luminosity distance to the source, in Gpc.

#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is given as a strain and is dimensionless.

#### Example:
```julia
    f = 10 .^range(1, stop=3, length=100)
    Ampl(TaylorF2(), f, 30., 8.)
```
"""
function Ampl(model::TaylorF2, f, mc, dL; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)

    # In the restricted PN approach the amplitude is the same as for the Newtonian approximation, so this term is equivalent
    amplitude = @. sqrt(5. /24.) * (pi^(-2. /3.)) * clightGpc/dL * (GMsun_over_c3*mc)^(5. /6.) * (f^(-7. /6.))
    return amplitude
end


"""
Compute the cut frequency of the waveform as a function of the events parameters, in Hz.

This can be approximated as 2 f_ISCO for inspiral only waveforms:

    - if Schwarzschild ISCO for a non-rotating final BH is used (default), use _fcut(model::TaylorF2, mc, eta) to compute the cut frequency of the waveform for the chosen events, in Hz.
    - if Kerr ISCO for a rotating final BH is computed, use _fcut(model::TaylorF2, mc, eta, chi1, chi2) to compute the cut frequency of the waveform for the chosen events, in Hz. As in `arXiv:2108.05861 <https://arxiv.org/abs/2108.05861>`_ (see in particular App. C). NOTE: this is pushing the validity of the model to the limit, and is not the default option.
#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `model`: Model of the waveform.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Return:
-  Cut frequency (float) of the waveform for the chosen event, in Hz.

"""
function _fcut(model::TaylorF2, mc, eta, chi1, chi2, Lambda1, Lambda2; GMsun_over_c3 = uc.GMsun_over_c3)

    println("Using fcut of Kerr orbit since the two angular momenta 'chi1' and 'chi2' were given")
    eta2 = eta*eta
    Mtot = mc/(eta^(3. /5.))
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1^2 * chi1 + m2^2 * chi2) / (m1^2 + m2^2)
    atot = (chi1 + chi2*(m2/m1)*(m2/m1))/((1. +m2/m1)*(1. +m2/m1))
    aeff = atot + 0.41616*eta*(chi1 + chi2)

    r_ISCO = _r_ISCO_of_chi(aeff)
    
    EradNS = eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta^3)
    EradTot = (EradNS * (1. + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s)) / (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s)
    
    Mfin = Mtot*(1. -EradTot)
    L_ISCO = 2. /(3. *sqrt(3.))*(1. + 2. *sqrt(3. *r_ISCO - 2.))
    E_ISCO = sqrt(1. - 2. /(3. *r_ISCO))
    
    chif = atot + eta*(L_ISCO - 2. *atot*(E_ISCO - 1.)) + (-3.821158961 - 1.2019*aeff - 1.20764*aeff^2)*eta2 + (3.79245 + 1.18385*aeff + 4.90494*aeff^3)*eta^3
    
    Om_ISCO = 1. /(((_r_ISCO_of_chi(chif))^(3. /2.))+chif)

    return Om_ISCO/(pi*Mfin*GMsun_over_c3)
end

function _r_ISCO_of_chi(chi)
    Z1_ISCO = 1.0 + ((1.0 - chi^2)^(1. /3.))*((1.0+chi)^(1. /3.) + (1.0-chi)^(1. /3.))
    Z2_ISCO = sqrt(3.0*chi*chi + Z1_ISCO^2)
    return ifelse(chi>0., 3.0 + Z2_ISCO - sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)), 3.0 + Z2_ISCO + sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)))
end
