#! format: off

####################################################
#   IMRPHENOM_HM WAVEFORM
####################################################

"""
ToDo: Need documentation
"""
function PolAbs(model::PhenomHM,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0
)
    hp, hc = waveform.hphc(model, f, mc, eta, chi1, chi2, dL, iota)

    return [abs.(hp), abs.(hc)]
    
end

"""
ToDo: Need documentation
"""
function Pol(model::PhenomHM,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    Lambda1=0.0,
    Lambda2=0.0
)
    hp, hc = waveform.hphc(model, f, mc, eta, chi1, chi2, dL, iota)
    
    # Return polarization with correct relative phase.
    # ATTENTION: For this waveform, the complete phase is already included 
    # in hphc() function. 
    return [hp, hc]

end

"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomHM(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns a structure containg the 6 modes of the phase. The phase is given in radians.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomHM(), f, mc, eta, chi1, chi2)
```
"""
function Phi(model::PhenomHM,
    f,
    mc,
    eta,
    chi1,
    chi2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)
    return f .* 0.0
end


"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomHM(), f, mc, eta, chi1, chi2, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the binary, in Gpc.
#### Optional arguments:
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns a structure containg the 6 modes of the amplitude. 
           To access the amplitude of the (2,2) mode, for example, you can use `Ampl(PhenomHM(), f, mc, eta, dL, chi1, chi2).two_two`. See also the julia function `fieldnames` to see all the available fields.
The amplitude is dimensionless. 
#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)
```
"""
function Ampl(model::PhenomHM,
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
)


    # Useful quantities
    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta
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
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    finMass = 1.0 - Erad

    fring, fdamp = _RDfreqCalc(model, finMass, aeff, 2, 2)

    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
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


    # Defined as in LALSimulation - LALSimIMRPhenomUtils.c line 70. Final units are correctly Hz^-1
    # there is a 2 * sqrt(5/(64*pi)) missing w.r.t the standard coefficient, which comes from the (2,2) shperical harmonic
    Overallamp = M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    tmpResults = Vector{Any}(undef, 6)
    i = 1
    for ell in (2:4)
        for mm in (ell - 1, ell)
            # Compute ringdown and damping frequencies from fits for the various modes
            fringlm, fdamplm = _RDfreqCalc(model, finMass, aeff, ell, mm)
            Rholm, Taulm = fring / fringlm, fdamplm / fdamp

            # Scale input frequencies according to PhenomHM model
            # Compute mapping coefficinets
            Map_fl = fInsJoin_Ampl
            Map_fi = Map_fl / Rholm
            Map_fr = fringlm

            Map_ai, Map_bi = 2.0 / mm, 0.0

            Map_Trd = Map_fr - fringlm + fring
            Map_Ti = 2.0 * Map_fi / mm
            Map_am = (Map_Trd - Map_Ti) / (Map_fr - Map_fi)
            Map_bm = Map_Ti - Map_fi * Map_am

            Map_ar, Map_br = 1.0, -Map_fr + fring

            # Now scale as f -> f*a+b for each regime
            fgridScaled = @. ifelse(
                fgrid < Map_fi,
                fgrid * Map_ai + Map_bi,
                ifelse(fgrid < Map_fr, fgrid * Map_am + Map_bm, fgrid * Map_ar + Map_br),
            )

            # Map the ampliude's range
            # We divide by the leading order l=m=2 behavior, and then scale in the expected PN behavior for the multipole of interest.

            beta_term1 = _onePointFiveSpinPN_Ampl(fgrid, ell, mm, chi_s, chi_a, eta, Seta)
            beta_term2 = _onePointFiveSpinPN_Ampl(
                2.0 .* fgrid ./ mm,
                ell,
                mm,
                chi_s,
                chi_a,
                eta,
                Seta,
            )

            HMamp_term1 =
                _onePointFiveSpinPN_Ampl(fgridScaled, ell, mm, chi_s, chi_a, eta, Seta)
            HMamp_term2 = _onePointFiveSpinPN_Ampl(fgridScaled, 2, 2, 0.0, 0.0, eta, Seta)

            # The (3,3) and (4,3) modes vanish if eta=0.25 (equal mass case) and the (2,1) mode vanishes if both eta=0.25 and chi1z=chi2z
            completeAmpl = @. Overallamp *
               amp0 *
               (fgridScaled^(-7.0 / 6.0)) *
               ifelse(
                   fgridScaled < fInsJoin_Ampl,
                   1.0 +
                   (fgridScaled^(2.0 / 3.0)) * Acoeffs.two_thirds +
                   (fgridScaled^(4.0 / 3.0)) * Acoeffs.four_thirds +
                   (fgridScaled^(5.0 / 3.0)) * Acoeffs.five_thirds +
                   (fgridScaled^(7.0 / 3.0)) * Acoeffs.seven_thirds +
                   (fgridScaled^(8.0 / 3.0)) * Acoeffs.eight_thirds +
                   fgridScaled * (
                       Acoeffs.one +
                       fgridScaled * Acoeffs.two +
                       fgridScaled * fgridScaled * Acoeffs.three
                   ),
                   ifelse(
                       fgridScaled < fpeak,
                       delta0 +
                       fgridScaled * delta1 +
                       fgridScaled *
                       fgridScaled *
                       (delta2 + fgridScaled * delta3 + fgridScaled * fgridScaled * delta4),
                       ifelse(
                           fgridScaled < fcutPar,
                           exp(-(fgridScaled - fring) * gamma2 / (fdamp * gamma3)) *
                           (fdamp * gamma3 * gamma1) / (
                               (fgridScaled - fring) * (fgridScaled - fring) +
                               fdamp * gamma3 * fdamp * gamma3
                           ),
                           0.0,
                       ),
                   ),
               )

            # This results in NaNs having 0/0, correct for this with replace!
            tmpResults[i] = replace!(
                completeAmpl .* (beta_term1 ./ beta_term2) .* HMamp_term1 ./ HMamp_term2,
                NaN => 0.0,
            )
            i += 1
        end
    end
    ampllm = AmplslmpStructure(
        tmpResults[1],
        tmpResults[2],
        tmpResults[3],
        tmpResults[4],
        tmpResults[5],
        tmpResults[6],
    )
    return ampllm
end




"""
Compute plus and cross polarisations of the GW as a function of frequency, given the events parameters. The function gives as output h_+ and h_x, which result from
the sum of the contributions of the different modes of the waveform. The function avoids for loops over the modes, dealing with both phase and amplitude in a single function.

    hphc(PhenomHM(), f, mc, eta, chi1, chi2, dL, iota)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `dL`: Luminosity distance to the source, in Gpc.
-  `iota`: Inclination angle of the binary. In radians.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate phase. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  (::Tuple{Vector{ComplexF64}, Vector{ComplexF64}}, so a tuple of the two polarisation) GW polarisations h_+ and h_x, as a function of frequency.
#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    iota = pi/3
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    hp, hc = hphc(PhenomHM(), f, mc, eta, chi1, chi2, dL, iota)
```
"""
function hphc(model::PhenomHM,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota;
    fcutPar = 0.2,
    fInsJoin_PHI = 0.018,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container= nothing,
    call_number = 1,
    optimization = false
)

    # This function retuns directly the full plus and cross polarisations, avoiding for loops over the modes

    complShiftm = [0.0, pi * 0.5, 0.0, -pi * 0.5, pi, pi * 0.5, 0.0]

    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta
    pi2 = pi * pi

    QuadMon1, QuadMon2 = 1.0, 1.0

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    q = 0.5 * (1.0 + Seta - 2.0 * eta) * etaInv
    chi_s2, chi_a2 = chi_s * chi_s, chi_a * chi_a
    chi1dotchi2 = chi1 * chi2
    chi_sdotchi_a = chi_s * chi_a
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # This is MfRef, needed to recover LAL, which sets fRef to f_min if fRef=0
    fRef = fgrid[1] 
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin, radiated energy and mass
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    finMass = 1.0 - Erad


    # Compute the real and imag parts of the complex ringdown frequency for the (l,m) mode as in LALSimIMRPhenomHM.c line 189
    # These are all fits of the different modes. We directly exploit the fact that the relevant HM in this WF are 6
    modes = [21, 22, 32, 33, 43, 44]
    ells = [2, 2, 3, 3, 4, 4]
    mms = [1, 2, 2, 3, 3, 4]
    # Domain mapping for dimensionless BH spin
    alphaRDfr = log(2.0 - aeff) / log(3.0)
    betaRDfr = [1.0 / 3.0, 0.5, 1.0 / 3.0, 0.5, 1.0 / 3.0, 0.5]
    kappaRDfr = alphaRDfr .^ betaRDfr
    kappaRDfr2 = kappaRDfr .* kappaRDfr
    kappaRDfr3 = kappaRDfr .* kappaRDfr2
    kappaRDfr4 = kappaRDfr .* kappaRDfr3

    tmpRDfr = @. ifelse(
        modes == 21,
        0.589113 * exp(0.043525 * 1im) +
        0.18896353 * exp(2.289868 * 1im) * kappaRDfr +
        1.15012965 * exp(5.810057 * 1im) * kappaRDfr2 +
        6.04585476 * exp(2.741967 * 1im) * kappaRDfr3 +
        11.12627777 * exp(5.844130 * 1im) * kappaRDfr4 +
        9.34711461 * exp(2.669372 * 1im) * kappaRDfr4 * kappaRDfr +
        3.03838318 * exp(5.791518 * 1im) * kappaRDfr4 * kappaRDfr2,
        ifelse(
            modes == 22,
            1.0 +
            kappaRDfr * (
                1.557847 * exp(2.903124 * 1im) +
                1.95097051 * exp(5.920970 * 1im) * kappaRDfr +
                2.09971716 * exp(2.760585 * 1im) * kappaRDfr2 +
                1.41094660 * exp(5.914340 * 1im) * kappaRDfr3 +
                0.41063923 * exp(2.795235 * 1im) * kappaRDfr4
            ),
            ifelse(
                modes == 32,
                1.022464 * exp(0.004870 * 1im) +
                0.24731213 * exp(0.665292 * 1im) * kappaRDfr +
                1.70468239 * exp(3.138283 * 1im) * kappaRDfr2 +
                0.94604882 * exp(0.163247 * 1im) * kappaRDfr3 +
                1.53189884 * exp(5.703573 * 1im) * kappaRDfr4 +
                2.28052668 * exp(2.685231 * 1im) * kappaRDfr4 * kappaRDfr +
                0.92150314 * exp(5.841704 * 1im) * kappaRDfr4 * kappaRDfr2,
                ifelse(
                    modes == 33,
                    1.5 +
                    kappaRDfr * (
                        2.095657 * exp(2.964973 * 1im) +
                        2.46964352 * exp(5.996734 * 1im) * kappaRDfr +
                        2.66552551 * exp(2.817591 * 1im) * kappaRDfr2 +
                        1.75836443 * exp(5.932693 * 1im) * kappaRDfr3 +
                        0.49905688 * exp(2.781658 * 1im) * kappaRDfr4
                    ),
                    ifelse(
                        modes == 43,
                        1.5 +
                        kappaRDfr * (
                            0.205046 * exp(0.595328 * 1im) +
                            3.10333396 * exp(3.016200 * 1im) * kappaRDfr +
                            4.23612166 * exp(6.038842 * 1im) * kappaRDfr2 +
                            3.02890198 * exp(2.826239 * 1im) * kappaRDfr3 +
                            0.90843949 * exp(5.915164 * 1im) * kappaRDfr4
                        ),
                        2.0 +
                        kappaRDfr * (
                            2.658908 * exp(3.002787 * 1im) +
                            2.97825567 * exp(6.050955 * 1im) * kappaRDfr +
                            3.21842350 * exp(2.877514 * 1im) * kappaRDfr2 +
                            2.12764967 * exp(5.989669 * 1im) * kappaRDfr3 +
                            0.60338186 * exp(2.830031 * 1im) * kappaRDfr4
                        ),
                    ),
                ),
            ),
        ),
    )

    fringlm = real(tmpRDfr) ./ (2.0 * pi * finMass)
    fdamplm = imag(tmpRDfr) ./ (2.0 * pi * finMass)

    

    fring = fringlm[2]
    fdamp = fdamplm[2]
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

    fInsJoinPh = fInsJoin_PHI
    fMRDJoinPh = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms
    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoinPh)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoinPh) / 3.0)
            ) * ((pi * fInsJoinPh)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoinPh)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoinPh)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoinPh) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoinPh)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoinPh)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoinPh)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoinPh^(1.0 / 3.0)) +
            sigma3 * (fInsJoinPh^(2.0 / 3.0)) +
            sigma4 * fInsJoinPh
        ) * etaInv
    # This is the first derivative of the Intermediate phase computed at fMRDJoin
    DPhiInt = (beta1 + beta3 / (fInsJoinPh^4) + beta2 / fInsJoinPh) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoinPh^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoinPh^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoinPh^(1.0 / 3.0)) * log(pi * fInsJoinPh) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoinPh) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoinPh^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoinPh^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoinPh +
        PhiInspcoeffs.min_four_thirds * (fInsJoinPh^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoinPh^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoinPh +
            PhiInspcoeffs.four_thirds * (fInsJoinPh^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoinPh^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoinPh * fInsJoinPh
        ) * etaInv
    # This is the Intermediate phase computed at fMRDJoinPh
    PhiIntJoin =
        beta1 * fInsJoinPh - beta3 / (3.0 * fInsJoinPh * fInsJoinPh * fInsJoinPh) +
        beta2 * log(fInsJoinPh)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoinPh
    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoinPh - beta3 / (3.0 * fMRDJoinPh^3) + beta2 * log(fMRDJoinPh)) *
        etaInv +
        C1Int +
        C2Int * fMRDJoinPh
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoinPh^4) + beta2 / fMRDJoinPh) * etaInv

    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoinPh^2) +
            alpha3 / (fMRDJoinPh^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoinPh - alpha5 * fring) * (fMRDJoinPh - alpha5 * fring) /
                    (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoinPh) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoinPh^(0.75))) +
        alpha1 * fMRDJoinPh +
        alpha4 * atan((fMRDJoinPh - alpha5 * fring) / fdamp)
    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoinPh

    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
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

    # Defined as in LALSimulation - LALSimIMRPhenomUtils.c line 70. Final units are correctly Hz^-1
    # there is a 2 * sqrt(5/(64*pi)) missing w.r.t the standard coefficient, which comes from the (2,2) shperical harmonic
    Overallamp = M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL




    # Time shift so that peak amplitude is approximately at t=0
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak * fpeak) +
            alpha3 / (fpeak^(1.0 / 4.0)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv

    phiRef = _completePhase(
        fRef,
        C1MRD,
        C2MRD,
        1.0,
        1.0,
        C1Int,
        C2Int,
        PhiInspcoeffs,
        fInsJoin_PHI,
        alpha1,
        alpha2,
        alpha3,
        alpha4,
        alpha5,
        beta1,
        beta2,
        beta3,
        etaInv,
        fMRDJoinPh,
        fcutPar,
        fring,
        fdamp,
    )

    phi0 = 0.5 * phiRef

    # Now compute all the modes, they are 6

    Rholm, Taulm = fring ./ fringlm, fdamplm ./ fdamp
    # Rholm and Taulm only figure in the MRD part, the rest of the coefficients is the same, recompute only this
    DPhiMRDVal = @. (
        alpha1 +
        alpha2 / (fMRDJoinPh * fMRDJoinPh) +
        alpha3 / (fMRDJoinPh^(1.0 / 4.0)) +
        alpha4 / (
            fdamp *
            Taulm *
            (
                1.0 +
                (fMRDJoinPh - alpha5 * fring) * (fMRDJoinPh - alpha5 * fring) /
                (fdamp * Taulm * Rholm * fdamp * Taulm * Rholm)
            )
        )
    ) * etaInv
    PhiMRJoinTemp = @. -(alpha2 / fMRDJoinPh) +
       (4.0 / 3.0) * (alpha3 * (fMRDJoinPh^(3.0 / 4.0))) +
       alpha1 * fMRDJoinPh +
       alpha4 * Rholm * atan((fMRDJoinPh - alpha5 * fring) / (fdamp * Rholm * Taulm))
    C2MRDHM = @. DPhiIntTempVal - DPhiMRDVal
    C1MRDHM = @. (PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRDHM * fMRDJoinPh)

    # Scale input frequencies according to PhenomHM model
    # Compute mapping coefficinets
    Map_flPhi = fInsJoin_PHI
    Map_fiPhi = Map_flPhi ./ Rholm
    Map_flAmp = fInsJoin_Ampl
    Map_fiAmp = Map_flAmp ./ Rholm
    Map_fr = fringlm

    Map_ai, Map_bi = 2.0 ./ mms, 0.0


    Map_TrdAmp = @. Map_fr - fringlm + fring
    Map_TiAmp = @. 2.0 * Map_fiAmp / mms
    Map_amAmp = @. (Map_TrdAmp - Map_TiAmp) / (Map_fr - Map_fiAmp)
    Map_bmAmp = @. Map_TiAmp - Map_fiAmp * Map_amAmp

    Map_TrdPhi = @. Map_fr * Rholm
    Map_TiPhi = @. 2.0 * Map_fiPhi / mms
    Map_amPhi = @. (Map_TrdPhi - Map_TiPhi) / (Map_fr - Map_fiPhi)
    Map_bmPhi = @. Map_TiPhi - Map_fiPhi * Map_amPhi

    Map_arAmp, Map_brAmp = 1.0, -Map_fr .+ fring
    Map_arPhi, Map_brPhi = Rholm, 0.0

    # Here we do a reshape of the input frequencies, to make the operations easier
    # We convert them to matrices (1-dim matrices) because we need to broadcast the operations
    # For example a 1-dim matrix of 1000x1 elements and a 1-dim matrix of 1x6 elements, when multiplied, will give a 1000x6 matrix
    # It can be hard to follow but remember that the matrices (1000,1) are frequencies or quantities independent of the modes, 
    # and the matrices (1,6) are the coefficients of the modes
    # The result is a matrix (1000,6) with the values of the modes for each frequency
    fgrid = reshape(fgrid, :, 1)
    Map_ai, Map_amAmp, Map_bmAmp, Map_brAmp = reshape(Map_ai, 1, :),
    reshape(Map_amAmp, 1, :),
    reshape(Map_bmAmp, 1, :),
    reshape(Map_brAmp, 1, :)

    Map_amPhi, Map_bmPhi, Map_arPhi =
        reshape(Map_amPhi, 1, :), reshape(Map_bmPhi, 1, :), reshape(Map_arPhi, 1, :)

    Map_fr = reshape(Map_fr, 1, :)
    mms = reshape(mms, 1, :)
    modes = reshape(modes, 1, :)
    Map_fiAmp = reshape(Map_fiAmp, 1, :)
    Map_fiPhi = reshape(Map_fiPhi, 1, :)
    fgridScaled = @. ifelse(
        fgrid < Map_fiAmp,
        fgrid * Map_ai + Map_bi,
        ifelse(
            fgrid < Map_fr,
            fgrid * Map_amAmp + Map_bmAmp,
            fgrid * Map_arAmp + Map_brAmp,
        ),
    )

    # Map the ampliude's range
    # We divide by the leading order l=m=2 behavior, and then scale in the expected PN behavior for the multipole of interest.

    beta_term1 = _onePointFiveSpinPN_hphc(fgrid, chi_s, chi_a, modes, mms, eta, Seta)
    beta_term2 =
        _onePointFiveSpinPN_hphc(2.0 * fgrid ./ mms, chi_s, chi_a, modes, mms, eta, Seta)
    HMamp_term1 = _onePointFiveSpinPN_hphc(fgridScaled, chi_s, chi_a, modes, mms, eta, Seta)


    HMamp_term2 = pi * sqrt(eta * 2.0 / 3.0) .* ((pi .* fgridScaled) .^ (-7.0 / 6.0))


    # The (3,3) and (4,3) modes vanish if eta=0.25 (equal mass case) and the (2,1) mode vanishes if both eta=0.25 and chi1z=chi2z
    # This results in NaNs having 0/0, correct for this using replace!()
    completeAmpl = @. Overallamp *
       amp0 *
       (fgridScaled^(-7.0 / 6.0)) *
       ifelse(
           fgridScaled < fInsJoin_Ampl,
           1.0 +
           (fgridScaled^(2.0 / 3.0)) * Acoeffs.two_thirds +
           (fgridScaled^(4.0 / 3.0)) * Acoeffs.four_thirds +
           (fgridScaled^(5.0 / 3.0)) * Acoeffs.five_thirds +
           (fgridScaled^(7.0 / 3.0)) * Acoeffs.seven_thirds +
           (fgridScaled^(8.0 / 3.0)) * Acoeffs.eight_thirds +
           fgridScaled * (
               Acoeffs.one +
               fgridScaled * Acoeffs.two +
               fgridScaled * fgridScaled * Acoeffs.three
           ),
           ifelse(
               fgridScaled < fpeak,
               delta0 +
               fgridScaled * delta1 +
               fgridScaled *
               fgridScaled *
               (delta2 + fgridScaled * delta3 + fgridScaled * fgridScaled * delta4),
               ifelse(
                   fgridScaled < fcutPar,
                   exp(-(fgridScaled - fring) * gamma2 / (fdamp * gamma3)) *
                   (fdamp * gamma3 * gamma1) / (
                       (fgridScaled - fring) * (fgridScaled - fring) +
                       fdamp * gamma3 * fdamp * gamma3
                   ),
                   0.0,
               ),
           ),
       )

    AmplsAllModes = replace!(
        completeAmpl .* (beta_term1 ./ beta_term2) .* HMamp_term1 ./ HMamp_term2,
        NaN => 0.0,
    )

    C1MRDHM, C2MRDHM, Rholm, Taulm = reshape(C1MRDHM, 1, :),
    reshape(C2MRDHM, 1, :),
    reshape(Rholm, 1, :),
    reshape(Taulm, 1, :)

    tmpMf = @. Map_amPhi * Map_fiPhi + Map_bmPhi
    tmpMf = reshape(tmpMf, 1, :)
    PhDBconst =
        1.0 ./ Map_amPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_arPhi * Map_fr + Map_brPhi
    tmpMf = reshape(tmpMf, 1, :)
    PhDCconst =
        1.0 ./ Map_arPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_ai * Map_fiPhi + Map_bi
    tmpMf = reshape(tmpMf, 1, :)
    PhDBAterm =
        1.0 ./ Map_ai .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_amPhi * Map_fr + Map_bmPhi
    tmpMf = reshape(tmpMf, 1, :)
    tmpphaseC =
        -PhDBconst .+ PhDBAterm .+
        1.0 ./ Map_amPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )



    PhisAllModes =
        ifelse.(
            fgrid .< Map_fiPhi,
            _completePhase(
                fgrid .* Map_ai .+ Map_bi,
                C1MRDHM,
                C2MRDHM,
                Rholm,
                Taulm,
                C1Int,
                C2Int,
                PhiInspcoeffs,
                fInsJoin_PHI,
                alpha1,
                alpha2,
                alpha3,
                alpha4,
                alpha5,
                beta1,
                beta2,
                beta3,
                etaInv,
                fMRDJoinPh,
                fcutPar,
                fring,
                fdamp,
            ) ./ Map_ai,
            ifelse.(
                fgrid .< Map_fr,
                -PhDBconst .+ PhDBAterm .+
                _completePhase(
                    fgrid .* Map_amPhi .+ Map_bmPhi,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin_PHI,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoinPh,
                    fcutPar,
                    fring,
                    fdamp,
                ) ./ Map_amPhi,
                -PhDCconst .+ tmpphaseC .+
                _completePhase(
                    fgrid .* Map_arPhi .+ Map_brPhi,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin_PHI,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoinPh,
                    fcutPar,
                    fring,
                    fdamp,
                ) ./ Map_arPhi,
            ),
        )
    PhisAllModes = @. PhisAllModes - t0 * (fgrid - fRef) - mms * phi0 + complShiftm[mms+1]
    Y, Ymstar = _spinWeighted_SphericalHarmonic(iota, modes)
    Ymstar = conj(Ymstar)
    Y, Ymstar = reshape(Y, 1, :), reshape(Ymstar, 1, :)
    ells = reshape(ells, 1, :)


    hp = vec(sum(
        AmplsAllModes .* exp.(-1im .* PhisAllModes) .*
        (0.5 .* (Y .+ ((-1) .^ ells) .* Ymstar)),
        dims = 2,
    ))  # vec() to obtain an array from a nx1 matrix
    hc = -vec(sum(
        AmplsAllModes .* exp.(-1im .* PhisAllModes) .*
        (-1im .* 0.5 .* (Y .- ((-1) .^ ells) .* Ymstar)),
        dims = 2,
    ))

    if typeof(eta) !== Float64 && !isnothing(container)

        for i in eachindex(fgrid)
            container[i] = real(hp[i]).value + 1im * imag(hp[i]).value
            container[i + length(fgrid)] = real(hc[i]).value + 1im * imag(hc[i]).value
        end    

    end
    
    if optimization == true

        if call_number == 1
            return [real(hp); imag(hp); real(hc); imag(hc)]
        else call_number == 2 || call_number == 3
            return [hp; hc]

        end
    else
        return hp, hc
    end

end

"""
useful functions for hphc()
"""
function _onePointFiveSpinPN_Ampl(infreqs, l, m, chi_s, chi_a, eta, Seta) 
    # PN amplitudes function, needed to scale

    v = @. (2.0 * pi * infreqs / m)^(1.0 / 3.0)
    v2 = v .* v
    v3 = v2 .* v

    if (2 == l) && (2 == m)
        Hlm = 1.0    # is it an array or a scalar? infreqs is an array, so is v

    elseif (2 == l) && (1 == m)
        Hlm = @. (sqrt(2.0) / 3.0) * (
            v * Seta - v2 * 1.5 * (chi_a + Seta * chi_s) +
            v3 * Seta * ((335.0 / 672.0) + (eta * 117.0 / 56.0)) +
            v3 *
            v *
            (
                chi_a * (3427.0 / 1344.0 - eta * 2101.0 / 336.0) +
                Seta * chi_s * (3427.0 / 1344 - eta * 965 / 336) +
                Seta * (-1im * 0.5 - pi - 2 * 1im * 0.69314718056)
            )
        )

    elseif (3 == l) && (3 == m)
        Hlm = @. 0.75 * sqrt(5.0 / 7.0) * (v * Seta)

    elseif (3 == l) && (2 == m)
        Hlm = @. (1.0 / 3.0) * sqrt(5.0 / 7.0) * (v2 * (1.0 - 3.0 * eta))

    elseif (4 == l) && (4 == m)
        Hlm = @. (4.0 / 9.0) * sqrt(10.0 / 7.0) * v2 * (1.0 - 3.0 * eta)

    elseif (4 == l) && (3 == m)
        Hlm = @. 0.75 * sqrt(3.0 / 35.0) * v3 * Seta * (1.0 - 2.0 * eta)

    else
        error("Mode not present in IMRPhenomHM waveform model.")
    end
    # Compute the final PN Amplitude at Leading Order in Mf

    return @. pi * sqrt(eta * 2.0 / 3.0) .* (v .^ (-3.5)) * abs.(Hlm)
end

"""
useful functions for hphc()
"""
function _completePhase(
    infreqs,
    C1MRDuse,
    C2MRDuse,
    RhoUse,
    TauUse,
    C1Int,
    C2Int,
    PhiInspcoeffs,
    fInsJoin_PHI,
    alpha1,
    alpha2,
    alpha3,
    alpha4,
    alpha5,
    beta1,
    beta2,
    beta3,
    etaInv,
    fMRDJoinPh,
    fcutPar,
    fring,
    fdamp,
)
    return @. ifelse(
        infreqs < fInsJoin_PHI,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (infreqs^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (infreqs^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (infreqs^(1.0 / 3.0)) * log(pi * infreqs) / 3.0 +
        PhiInspcoeffs.log * log(pi * infreqs) / 3.0 +
        PhiInspcoeffs.min_third * (infreqs^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (infreqs^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / infreqs +
        PhiInspcoeffs.min_four_thirds * (infreqs^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (infreqs^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * infreqs +
            PhiInspcoeffs.four_thirds * (infreqs^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (infreqs^(5.0 / 3.0)) +
            PhiInspcoeffs.two * infreqs * infreqs
        ) * etaInv,
        ifelse(
            infreqs < fMRDJoinPh,
            (beta1 * infreqs - beta3 / (3.0 * infreqs^3) + beta2 * log(infreqs)) * etaInv +
            C1Int +
            C2Int * infreqs,
            ifelse(
                infreqs < fcutPar,
                (
                    -(alpha2 / infreqs) +
                    (4.0 / 3.0) * (alpha3 * (infreqs^(3.0 / 4.0))) +
                    alpha1 * infreqs +
                    alpha4 *
                    RhoUse *
                    atan((infreqs - alpha5 * fring) / (fdamp * RhoUse * TauUse))
                ) * etaInv +
                C1MRDuse +
                C2MRDuse * infreqs,
                0.0,
            ),
        ),
    )
end

"""
useful functions for hphc()
"""
function _onePointFiveSpinPN_hphc(infreqs, chi_s, chi_a, modes, mms, eta, Seta) #this one is different from the one in Ampl since it computes ones for all the multipoles, to distinguish look at the number of inputs
    # PN amplitudes function, needed to scale

    v = @. (2.0 * pi * infreqs / mms)^(1.0 / 3.0)
    v2 = v .* v
    v3 = v2 .* v
    Hlm = @. ifelse(
        modes == 21,
        (sqrt(2.0) / 3.0) * (
            v * Seta - v2 * 1.5 * (chi_a + Seta * chi_s) +
            v3 * Seta * ((335.0 / 672.0) + (eta * 117.0 / 56.0)) +
            v3 *
            v *
            (
                chi_a * (3427.0 / 1344.0 - eta * 2101.0 / 336.0) +
                Seta * chi_s * (3427.0 / 1344 - eta * 965 / 336) +
                Seta * (-1im * 0.5 - pi - 2 * 1im * 0.69314718056)
            )
        ),
        ifelse(
            modes == 22,
            1.0,
            ifelse(
                modes == 32,
                (1.0 / 3.0) * sqrt(5.0 / 7.0) * (v2 * (1.0 - 3.0 * eta)),
                ifelse(
                    modes == 33,
                    0.75 * sqrt(5.0 / 7.0) * (v * Seta),
                    ifelse(
                        modes == 43,
                        0.75 * sqrt(3.0 / 35.0) * v3 * Seta * (1.0 - 2.0 * eta),
                        (4.0 / 9.0) * sqrt(10.0 / 7.0) * v2 * (1.0 - 3.0 * eta),
                    ),
                ),
            ),
        ),
    )

    # Compute the final PN Amplitude at Leading Order in Mf

    return @. pi * sqrt(eta * 2.0 / 3.0) * (v^(-3.5)) * abs(Hlm)
end


"""
useful functions for hphc()
"""
function _RDfreqCalc(model::PhenomHM, finalmass, finalspin, l, m)

    # Domain mapping for dimenionless BH spin
    alpha = log(2.0 - finalspin) / log(3.0)
    beta = 1.0 / (2.0 + l - abs(m))
    kappa = alpha^beta
    kappa2 = kappa * kappa
    kappa3 = kappa * kappa2
    kappa4 = kappa * kappa3

    if (2 == l) && (2 == m)
        res =
            1.0 +
            kappa * (
                1.557847 * exp(2.903124 * 1im) +
                1.95097051 * exp(5.920970 * 1im) * kappa +
                2.09971716 * exp(2.760585 * 1im) * kappa2 +
                1.41094660 * exp(5.914340 * 1im) * kappa3 +
                0.41063923 * exp(2.795235 * 1im) * kappa4
            )

    elseif (3 == l) && (2 == m)
        res =
            1.022464 * exp(0.004870 * 1im) +
            0.24731213 * exp(0.665292 * 1im) * kappa +
            1.70468239 * exp(3.138283 * 1im) * kappa2 +
            0.94604882 * exp(0.163247 * 1im) * kappa3 +
            1.53189884 * exp(5.703573 * 1im) * kappa4 +
            2.28052668 * exp(2.685231 * 1im) * kappa4 * kappa +
            0.92150314 * exp(5.841704 * 1im) * kappa4 * kappa2

    elseif (4 == l) && (4 == m)
        res =
            2.0 +
            kappa * (
                2.658908 * exp(3.002787 * 1im) +
                2.97825567 * exp(6.050955 * 1im) * kappa +
                3.21842350 * exp(2.877514 * 1im) * kappa2 +
                2.12764967 * exp(5.989669 * 1im) * kappa3 +
                0.60338186 * exp(2.830031 * 1im) * kappa4
            )

    elseif (2 == l) && (1 == m)
        res =
            0.589113 * exp(0.043525 * 1im) +
            0.18896353 * exp(2.289868 * 1im) * kappa +
            1.15012965 * exp(5.810057 * 1im) * kappa2 +
            6.04585476 * exp(2.741967 * 1im) * kappa3 +
            11.12627777 * exp(5.844130 * 1im) * kappa4 +
            9.34711461 * exp(2.669372 * 1im) * kappa4 * kappa +
            3.03838318 * exp(5.791518 * 1im) * kappa4 * kappa2

    elseif (3 == l) && (3 == m)
        res =
            1.5 +
            kappa * (
                2.095657 * exp(2.964973 * 1im) +
                2.46964352 * exp(5.996734 * 1im) * kappa +
                2.66552551 * exp(2.817591 * 1im) * kappa2 +
                1.75836443 * exp(5.932693 * 1im) * kappa3 +
                0.49905688 * exp(2.781658 * 1im) * kappa4
            )

    elseif (4 == l) && (3 == m)
        res =
            1.5 +
            kappa * (
                0.205046 * exp(0.595328 * 1im) +
                3.10333396 * exp(3.016200 * 1im) * kappa +
                4.23612166 * exp(6.038842 * 1im) * kappa2 +
                3.02890198 * exp(2.826239 * 1im) * kappa3 +
                0.90843949 * exp(5.915164 * 1im) * kappa4
            )

    else
        error("Mode not present in IMRPhenomHM waveform model.")
    end

    if m < 0
        res = -conj(res)
    end
    fring = real(res) / (2.0 * pi * finalmass)

    fdamp = imag(res) / (2.0 * pi * finalmass)

    return fring, fdamp
end
