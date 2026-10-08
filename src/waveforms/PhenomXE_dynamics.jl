#! format: off
#######################################################
# IMRPhenomXE: eccentric quasi-Keplerian dynamics
#######################################################
#
# Port of phenomxpy/eccentricity (secular_eqs_eob.py, secular_eqs_pn.py, eccentric_dynamics.py, compiled_dop853.py)
# and of the dynamics part of phenomxpy/phenomxe (spa_xe_num.py, internals_xe.py), restricted to the default
# IMRPhenomXE configuration (rhs_eqs = "eob_phenomT"):
#   * 3PN EOB quasi-Keplerian secular equations for (x, e, l, lambda) without the quasi-circular x-dot
#     (_rhs_3pn_eob_noqc) and conversions between the EOB and the PN (harmonic) eccentricities;
#   * a DOP853 integrator with the same step controller as SciPy;
#   * not-a-knot / PCHIP cubic splines with SciPy conventions.
#
# The functions _rhs_3pn_eob_noqc, _eEOB_from_et, _et_from_eEOB, _compute_dldt_3pn and _compute_dldtnorm_3pn were
# transliterated automatically from the Python sources and checked numerically against them.

const _XE_PI2 = pi * pi
const _XE_LOG2 = log(2.0)
const _XE_LOG3 = log(3.0)
const _XE_LOG5 = log(5.0)
const _XE_EULERGAMMA = MathConstants.eulergamma

# DOP853 Runge-Kutta tableau (Hairer, Norsett & Wanner), identical to SciPy's DOP853 implementation.

const DOP853_A = [
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.05260015195876773 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0197250569845379 0.0591751709536137 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.02958758547680685 0.0 0.08876275643042054 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.2413651341592667 0.0 -0.8845494793282861 0.924834003261792 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.037037037037037035 0.0 0.0 0.17082860872947386 0.12546768756682242 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.037109375 0.0 0.0 0.17025221101954405 0.06021653898045596 -0.017578125 0.0 0.0 0.0 0.0 0.0 0.0;
    0.03709200011850479 0.0 0.0 0.17038392571223998 0.10726203044637328 -0.015319437748624402 0.008273789163814023 0.0 0.0 0.0 0.0 0.0;
    0.6241109587160757 0.0 0.0 -3.3608926294469414 -0.868219346841726 27.59209969944671 20.154067550477894 -43.48988418106996 0.0 0.0 0.0 0.0;
    0.47766253643826434 0.0 0.0 -2.4881146199716677 -0.590290826836843 21.230051448181193 15.279233632882423 -33.28821096898486 -0.020331201708508627 0.0 0.0 0.0;
    -0.9371424300859873 0.0 0.0 5.186372428844064 1.0914373489967295 -8.149787010746927 -18.52006565999696 22.739487099350505 2.4936055526796523 -3.0467644718982196 0.0 0.0;
    2.273310147516538 0.0 0.0 -10.53449546673725 -2.0008720582248625 -17.9589318631188 27.94888452941996 -2.8589982771350235 -8.87285693353063 12.360567175794303 0.6433927460157636 0.0
]

const DOP853_B = [0.054293734116568765, 0.0, 0.0, 0.0, 0.0, 4.450312892752409, 1.8915178993145003, -5.801203960010585, 0.3111643669578199, -0.1521609496625161, 0.20136540080403034, 0.04471061572777259]

const DOP853_C = [0.0, 0.05260015195876773, 0.0789002279381516, 0.1183503419072274, 0.2816496580927726, 0.3333333333333333, 0.25, 0.3076923076923077, 0.6512820512820513, 0.6, 0.8571428571428571, 1.0]

const DOP853_E3 = [-0.18980075407240762, 0.0, 0.0, 0.0, 0.0, 4.450312892752409, 1.8915178993145003, -5.801203960010585, -0.4226823213237919, -0.1521609496625161, 0.20136540080403034, 0.02265179219836082, 0.0]

const DOP853_E5 = [0.01312004499419488, 0.0, 0.0, 0.0, 0.0, -1.2251564463762044, -0.4957589496572502, 1.6643771824549864, -0.35032884874997366, 0.3341791187130175, 0.08192320648511571, -0.022355307863886294, 0.0]

function _rhs_3pn_eob_noqc(x, e, eta, delta, chiS, chiA)
    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    e10 = e8 * e2
    eta2 = eta * eta
    eta3 = eta2 * eta
    chiA2 = chiA * chiA
    chiS2 = chiS * chiS
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x
    xsqrt = sqrt(x)
    x32 = x * xsqrt
    x52 = x2 * xsqrt
    e2m1 = 1.0 - e2
    e2m12 = e2m1 * e2m1
    e2m13 = e2m12 * e2m1
    e2m1sqrt = e2m1^0.5
    e2m132 = e2m1 * e2m1sqrt
    e2m152 = e2m12 * e2m1sqrt
    e2m172 = e2m152 * e2m1
    loge2m1p1 = log((1.0 + e2m1sqrt) / (2.0 * xsqrt * e2m1))
    log1xsqrt = log(1.0 / xsqrt)
    loge2m1sqrt2 = log((1.0 + e2m1sqrt)^2)
    LOG2xe2m1 = log(2.0 * xsqrt * e2m1)
    fact0 = 64.0 * eta * x5 / (5.0 * e2m172)
    facte0 = 64.0 * eta * x5 / 5.0
    xdot0PN = 1.0 + 73.0 * e2 / 24.0 + 37.0 * e4 / 96.0
    xdot0PN_e0 = 1.0
    xdot1PN = -0.00018601190476190475 * ((16.0 * (743.0 + 924.0 * eta) + e6 * (6931.0 + 2072.0 * eta) + e4 * (99106.0 + 51660.0 * eta) + e2 * (123288.0 + 89264.0 * eta)) * x) / e2m1
    xdot1PN_e0 = -0.002976190476190476 * ((743.0 + 924.0 * eta) * x)
    xdot15PN = (-86784.0 * chiA * delta + 768.0 * chiS * (-113.0 + 76.0 * eta) + 36864.0 * pi + e6 * (11952.0 * chiA * delta - 144.0 * chiS * (-83.0 + 8.0 * eta) + 10007.0 * pi) + e4 * (41632.0 * chiA * delta + chiS * (41632.0 + 65152.0 * eta) + 188880.0 * pi) + e2 * (-200576.0 * chiA * delta + 128.0 * chiS * (-1567.0 + 1670.0 * eta) + 264000.0 * pi)) * x32 / (9216.0 * e2m132)
    xdot15PN_e0 = (-86784.0 * chiA * delta + 768.0 * chiS * (-113.0 + 76.0 * eta) + 36864.0 * pi) * x32 / 9216.0
    xdot2PN = (5878656.0 * chiA * chiS * delta - 36288.0 * chiS2 * (-81.0 + 4.0 * eta) - 36288.0 * chiA2 * (-81.0 + 320.0 * eta) + e2m1sqrt * (-290304.0 * (-5.0 + 2.0 * eta) - 1512000.0 * e2 * (-5.0 + 2.0 * eta) + 1324512.0 * e4 * (-5.0 + 2.0 * eta) + 477792.0 * e6 * (-5.0 + 2.0 * eta)) + 3.0 * e8 * (734703.0 + 290664.0 * eta + 58016.0 * eta2) + 32.0 * (-11257.0 + 141093.0 * eta + 59472.0 * eta2) + e6 * (-902664.0 * chiA * chiS * delta - 2268.0 * chiS2 * (199.0 + 36.0 * eta) + 2268.0 * chiA2 * (-199.0 + 832.0 * eta) + 6.0 * (5905155.0 + 6204657.0 * eta + 1292312.0 * eta2)) + e2 * (10463040.0 * chiA * chiS * delta - 6048.0 * chiS2 * (-865.0 + 228.0 * eta) - 6048.0 * chiA2 * (-865.0 + 3232.0 * eta) + 16.0 * (-2678686.0 + 5601690.0 * eta + 1331295.0 * eta2)) + e4 * (-5954256.0 * chiA * chiS * delta - 1512.0 * chiS2 * (1969.0 + 828.0 * eta) + 1512.0 * chiA2 * (-1969.0 + 8704.0 * eta) + 12.0 * (896914.0 + 11637378.0 * eta + 2585233.0 * eta2))) * x2 / (580608.0 * e2m12)
    xdot2PN_e0 = (5878656.0 * chiA * chiS * delta - 290304.0 * (-5.0 + 2.0 * eta) - 36288.0 * chiS2 * (-81.0 + 4.0 * eta) - 36288.0 * chiA2 * (-81.0 + 320.0 * eta) + 32.0 * (-11257.0 + 141093.0 * eta + 59472.0 * eta2)) * x2 / 580608.0
    xdot25PN = -2.583498677248677e-06 * ((-384.0 * chiA * delta * (-31319.0 + 48678.0 * eta) + 384.0 * chiS * (31319.0 - 91900.0 * eta + 26544.0 * eta2) + e2m132 * (e4 * (-1300992.0 * chiA * delta * (-3.0 + eta) + 650496.0 * chiS * (6.0 - 8.0 * eta + eta2)) + e2 * (-3268608.0 * chiA * delta * (-3.0 + eta) + 1634304.0 * chiS * (6.0 - 8.0 * eta + eta2))) + 576.0 * (4159.0 + 15876.0 * eta) * pi + e8 * (27.0 * chiA * delta * (125467.0 + 27916.0 * eta) + 27.0 * chiS * (125467.0 + 41528.0 * eta - 8008.0 * eta2) + 7.0 * (-151281.0 + 10007.0 * eta) * pi) + e2 * (-96.0 * chiA * delta * (642396.0 + 696941.0 * eta) + 96.0 * chiS * (-642396.0 - 125759.0 * eta + 632758.0 * eta2) + 576.0 * (115991.0 + 171104.0 * eta) * pi) + e6 * (12.0 * chiA * delta * (2448332.0 + 1062285.0 * eta) + 12.0 * chiS * (2448332.0 + 3133943.0 * eta + 337722.0 * eta2) + 4.0 * (17257369.0 + 7471473.0 * eta) * pi) + e4 * (-672.0 * chiA * delta * (80924.0 + 13833.0 * eta) + 1344.0 * chiS * (-40462.0 + 79817.0 * eta + 36774.0 * eta2) + 9.0 * (18599341.0 + 14964748.0 * eta) * pi)) * x52) / e2m152
    xdot25PN_e0 = -2.583498677248677e-06 * ((-384.0 * chiA * delta * (-31319.0 + 48678.0 * eta) + 384.0 * chiS * (31319.0 - 91900.0 * eta + 26544.0 * eta2) + 576.0 * (4159.0 + 15876.0 * eta) * pi) * x52)
    xdot3PN = -1.8639961596310804e-11 * ((-79833600.0 * chiS2 * (55817.0 - 117201.0 * eta + 33460.0 * eta2) - 79833600.0 * chiA2 * (55817.0 - 243029.0 * eta + 59136.0 * eta2) + 175.0 * e10 * (3116030391.0 + 1005041664.0 * eta + 405402624.0 * eta2 + 58347520.0 * eta3) + (-874720788480.0 - 12391877836800.0 * e2 - 23558235402240.0 * e4 - 7977271357440.0 * e6 - 253703197440.0 * e8) * loge2m1p1 - 8941363200.0 * chiS * (-225.0 + 148.0 * eta) * pi + chiA * (159667200.0 * chiS * delta * (-55817.0 + 68481.0 * eta) + 2011806720000.0 * delta * pi) + 640.0 * (-15399771333.0 + 22047056185.0 * eta + 501236505.0 * eta2 + 181265700.0 * eta3 + 1366751232.0 * _XE_EULERGAMMA + 2733502464.0 * _XE_LOG2 - 436590.0 * (1024.0 + 1845.0 * eta) * _XE_PI2) + e8 * (-1247400.0 * chiS2 * (503031.0 - 226408.0 * eta + 57792.0 * eta2) - 1247400.0 * chiA2 * (503031.0 - 2069104.0 * eta + 112000.0 * eta2) - 1140955200.0 * chiS * (-2.0 + eta) * pi + chiA * (2494800.0 * chiS * delta * (-503031.0 + 141694.0 * eta) + 2281910400.0 * delta * pi) + 14.0 * (358019945973.0 + 906482016000.0 * eta + 351807667200.0 * eta2 + 50935808000.0 * eta3 + 18121656960.0 * _XE_EULERGAMMA + 36243313920.0 * _XE_LOG2 + 92619450.0 * (-64.0 + 41.0 * eta) * _XE_PI2)) + e2 * (-4435200.0 * chiS2 * (-593801.0 - 5790757.0 * eta + 1987972.0 * eta2) - 4435200.0 * chiA2 * (-593801.0 + 851591.0 * eta + 4835712.0 * eta2) - 745113600.0 * chiS * (-14684.0 + 13789.0 * eta) * pi + chiA * (8870400.0 * chiS * delta * (593801.0 + 3657185.0 * eta) + 10941248102400.0 * delta * pi) + 32.0 * (-3181561351866.0 - 9693202750.0 * eta + 636423999750.0 * eta2 + 75206939500.0 * eta3 + 387246182400.0 * _XE_EULERGAMMA + 72893399040.0 * _XE_LOG2 + 700566783840.0 * _XE_LOG3 - 1819125.0 * (69632.0 + 9717.0 * eta) * _XE_PI2)) + e4 * (-1108800.0 * chiS2 * (-8736761.0 - 8677333.0 * eta + 5503540.0 * eta2) - 1108800.0 * chiA2 * (-8736761.0 + 31182647.0 * eta + 8270976.0 * eta2) - 93139200.0 * chiS * (4957.0 + 101614.0 * eta) * pi + chiA * (2217600.0 * chiS * delta * (8736761.0 + 6220865.0 * eta) - 461691014400.0 * delta * pi) + 8.0 * (-28536072962442.0 + 1156154672750.0 * eta + 7217960245650.0 * eta2 + 913931903500.0 * eta3 + 2944779425280.0 * _XE_EULERGAMMA + 65437201589760.0 * _XE_LOG2 - 31525505272800.0 * _XE_LOG3 + 363825.0 * (-2647552.0 + 591261.0 * eta) * _XE_PI2)) + e6 * (-1663200.0 * chiA2 * (2460015.0 - 10429817.0 * eta + 589232.0 * eta2) - 1663200.0 * chiS2 * (2460015.0 - 1205477.0 * eta + 797748.0 * eta2) - 3880800.0 * chiS * (1359484.0 + 214925.0 * eta) * pi + chiA * (3326400.0 * chiS * delta * (-2460015.0 + 897617.0 * eta) - 5275885507200.0 * delta * pi) + 12.0 * (664772613120.0 * _XE_EULERGAMMA + 2.0 * (-3371321595729.0 + 2173303676775.0 * eta + 1510891888275.0 * eta2 + 207243883000.0 * eta3 - 220951471910400.0 * _XE_LOG2 + 53938777308570.0 * _XE_LOG3 + 60344238281250.0 * _XE_LOG5) + 121275.0 * (-1793024.0 + 783633.0 * eta) * _XE_PI2)) + e2m1sqrt * (-35765452800.0 * chiA * chiS * delta * (-11.0 + 9.0 * eta) + 17882726400.0 * chiS2 * (11.0 - 16.0 * eta + 4.0 * eta2) + 17882726400.0 * chiA2 * (11.0 - 46.0 * eta + 6.0 * eta2) + 3326400.0 * e8 * (-428445.0 + 70634.0 * eta + 50176.0 * eta2) + 177408.0 * (19954466.0 - 1990800.0 * eta2 + 75.0 * eta * (-19748.0 + 861.0 * _XE_PI2)) + e4 * (-109531699200.0 * chiA * chiS * delta * (-11.0 + 9.0 * eta) + 54765849600.0 * chiS2 * (11.0 - 16.0 * eta + 4.0 * eta2) + 54765849600.0 * chiA2 * (11.0 - 46.0 * eta + 6.0 * eta2) + 33264.0 * (652068196.0 + 39620000.0 * eta2 + 175.0 * eta * (-960472.0 + 6027.0 * _XE_PI2))) + e6 * (-121453516800.0 * chiA * chiS * delta * (-11.0 + 9.0 * eta) + 60726758400.0 * chiS2 * (11.0 - 16.0 * eta + 4.0 * eta2) + 60726758400.0 * chiA2 * (11.0 - 46.0 * eta + 6.0 * eta2) + 11088.0 * (-558162656.0 + 239618400.0 * eta2 + 25.0 * eta * (-18306272.0 + 140343.0 * _XE_PI2))) + e2 * (266750668800.0 * chiA * chiS * delta * (-11.0 + 9.0 * eta) - 133375334400.0 * chiS2 * (11.0 - 16.0 * eta + 4.0 * eta2) - 133375334400.0 * chiA2 * (11.0 - 46.0 * eta + 6.0 * eta2) - 22176.0 * (-1237347456.0 + 170839200.0 * eta2 + 25.0 * eta * (-19288240.0 + 154119.0 * _XE_PI2))))) * x3) / e2m13
    xdot3PN_e0 = -1.8639961596310804e-11 * ((-35765452800.0 * chiA * chiS * delta * (-11.0 + 9.0 * eta) + 17882726400.0 * chiS2 * (11.0 - 16.0 * eta + 4.0 * eta2) + 17882726400.0 * chiA2 * (11.0 - 46.0 * eta + 6.0 * eta2) - 79833600.0 * chiS2 * (55817.0 - 117201.0 * eta + 33460.0 * eta2) - 79833600.0 * chiA2 * (55817.0 - 243029.0 * eta + 59136.0 * eta2) - 874720788480.0 * log1xsqrt - 8941363200.0 * chiS * (-225.0 + 148.0 * eta) * pi + chiA * (159667200.0 * chiS * delta * (-55817.0 + 68481.0 * eta) + 2011806720000.0 * delta * pi) + 640.0 * (-15399771333.0 + 22047056185.0 * eta + 501236505.0 * eta2 + 181265700.0 * eta3 + 1366751232.0 * _XE_EULERGAMMA + 2733502464.0 * _XE_LOG2 - 436590.0 * (1024.0 + 1845.0 * eta) * _XE_PI2) + 177408.0 * (19954466.0 - 1990800.0 * eta2 + 75.0 * eta * (-19748.0 + 861.0 * _XE_PI2))) * x3)
    xdot_e0 = facte0 * (xdot0PN_e0 + xdot1PN_e0 + xdot15PN_e0 + xdot2PN_e0 + xdot25PN_e0 + xdot3PN_e0)
    xdot = fact0 * (xdot0PN + xdot1PN + xdot15PN + xdot2PN + xdot25PN + xdot3PN)
    xdot_diff = xdot - xdot_e0
    fact0 = -304.0 * e * eta * x4 / (15.0 * e2m152)
    edot0PN = 1.0 + 121.0 * e2 / 304.0
    edot1PN = -1.9580200501253132e-05 * ((164376.0 + 464376.0 * e2 + 94887.0 * e4 + 196448.0 * eta + 257124.0 * e2 * eta + 19768.0 * e4 * eta) * x) / e2m1
    edot15PN = (-16.0 * chiA * delta * (16232.0 + 248.0 * e2 - 1869.0 * e4) + 16.0 * chiS * (-16232.0 + 1869.0 * e4 + 13312.0 * eta - 276.0 * e4 * eta + 4.0 * e2 * (-62.0 + 1993.0 * eta)) + (189120.0 + 286512.0 * e2 + 24217.0 * e4 - 98.0 * e6) * pi) * x32 / (29184.0 * e2m132)
    edot2PN = (5423040.0 * chiA * chiS * delta - 5040.0 * (304.0 + 121.0 * e2) * e2m132 * (-5.0 + 2.0 * eta) - 10080.0 * chiS2 * (-269.0 + 36.0 * eta) - 10080.0 * chiA2 * (-269.0 + 1040.0 * eta) - 12.0 * e2 * (949135.0 + 148008.0 * chiA * chiS * delta - 4137363.0 * eta + 84.0 * chiS2 * (881.0 + 540.0 * eta) - 84.0 * chiA2 * (-881.0 + 4064.0 * eta) - 987567.0 * eta2) + 3.0 * e6 * (1056441.0 + 339608.0 * eta + 60256.0 * eta2) + 16.0 * (-765197.0 + 772695.0 * eta + 225792.0 * eta2) + e4 * (22189718.0 - 485352.0 * chiA * chiS * delta + 27599862.0 * eta - 2268.0 * chiS2 * (107.0 + 20.0 * eta) + 2268.0 * chiA2 * (-107.0 + 448.0 * eta) + 5073936.0 * eta2)) * x2 / (612864.0 * e2m12)
    edot25PN = -4.079208437761069e-07 * ((2.0 * chiA * delta * (-3691584.0 + 29552328.0 * e4 + 10324017.0 * e6 + 24.0 * e2 * (-4476734.0 + 204960.0 * e2m1sqrt * (-3.0 + eta) - 1745023.0 * eta) + 26880.0 * e2m1sqrt * (-304.0 + 121.0 * e4) * (-3.0 + eta) - 60983552.0 * eta + 19163060.0 * e4 * eta + 2307732.0 * e6 * eta) + 2.0 * chiS * (-3691584.0 + 29552328.0 * e4 + 10324017.0 * e6 - 63079040.0 * eta + 67302476.0 * e4 * eta + 3203016.0 * e6 * eta + 13440.0 * e2m1sqrt * (-304.0 + 121.0 * e4) * (-6.0 + 8.0 * eta - eta2) + 42913024.0 * eta2 + 9741032.0 * e4 * eta2 - 817656.0 * e6 * eta2 + 24.0 * e2 * (-4476734.0 + 3800571.0 * eta + 102480.0 * e2m1sqrt * (-6.0 + 8.0 * eta - eta2) + 3015082.0 * eta2)) + (50657472.0 + 8626445.0 * e6 - 12348.0 * e8 + 110954496.0 * eta + 1850510.0 * e6 * eta + 1372.0 * e8 * eta + 25.0 * e4 * (9908007.0 + 4499656.0 * eta) + 12.0 * e2 * (25180851.0 + 26615044.0 * eta)) * pi) * x52) / e2m152
    edot3PN = -5.886303661992885e-12 * ((4435200.0 * chiA2 * (-2069461.0 + 9635809.0 * eta - 4395216.0 * eta2) + 4435200.0 * chiS2 * (-2069461.0 + 6612781.0 * eta - 2014516.0 * eta2) + 175.0 * e8 * (46932564429.0 - 4051745280.0 * eta - 2937623040.0 * eta2 + 178984960.0 * eta3) + 142369920.0 * (1.0 + e2m1sqrt) * (24608.0 + 89024.0 * e2 + 42884.0 * e4 + 1719.0 * e6) * (2.0 * LOG2xe2m1 - loge2m1sqrt2) - 1490227200.0 * chiS * (-6148.0 + 4937.0 * eta) * pi + 8870400.0 * chiA * delta * (chiS * (-2069461.0 + 3985373.0 * eta) + 1032864.0 * pi) + 32.0 * e2 * (-6068102350278.0 - 1938077318100.0 * eta + 69300.0 * chiS2 * (9106018.0 + 4372739.0 * eta - 3135692.0 * eta2) + 1564141592475.0 * eta2 - 69300.0 * chiA2 * (-9106018.0 + 34227913.0 * eta + 5762232.0 * eta2) + 186429656875.0 * eta3 + 792146234880.0 * _XE_EULERGAMMA + 7937977259520.0 * _XE_LOG2 - 2860647700680.0 * _XE_LOG3 - 5821200.0 * chiS * (-39223.0 + 74954.0 * eta) * pi + 138600.0 * chiA * delta * (chiS * (9106018.0 + 3284449.0 * eta) + 1647366.0 * pi) + 5821200.0 * (-44512.0 + 8077.0 * eta) * _XE_PI2) + 4.0 * e6 * (4502397710583.0 + 11968410346200.0 * eta + 103950.0 * chiA2 * (-6397597.0 + 26767760.0 * eta - 3295488.0 * eta2) + 103950.0 * chiS2 * (-6397597.0 + 6830264.0 * eta - 1593536.0 * eta2) + 1488261297600.0 * eta2 + 392650720000.0 * eta3 + 122366946240.0 * _XE_EULERGAMMA - 584404081430400.0 * _XE_LOG2 + 159670846150200.0 * _XE_LOG3 + 144826171875000.0 * _XE_LOG5 + 207900.0 * chiA * delta * (chiS * (-6397597.0 + 4003818.0 * eta) - 618044.0 * pi) + 970200.0 * chiS * (-132438.0 + 42287.0 * eta) * pi + 363825.0 * (-110016.0 + 9061.0 * eta) * _XE_PI2) + 64.0 * (-641828882523.0 + 380870611275.0 * eta + 30985883775.0 * eta2 + 12481353500.0 * eta3 + 109482468480.0 * _XE_EULERGAMMA + 102648712320.0 * _XE_LOG2 + 116761130640.0 * _XE_LOG3 - 727650.0 * (49216.0 + 17917.0 * eta) * _XE_PI2) + 12.0 * e4 * (138600.0 * chiA2 * (-3131923.0 + 13837483.0 * eta - 1693552.0 * eta2) + 138600.0 * chiS2 * (-3131923.0 + 3055871.0 * eta - 1466108.0 * eta2) + 92400.0 * chiA * delta * (chiS * (-9395769.0 + 6548493.0 * eta) - 7879319.0 * pi) - 646800.0 * chiS * (1125617.0 + 95561.0 * eta) * pi + 6.0 * (-2864930559373.0 + 864070982775.0 * eta2 + 99510719000.0 * eta3 + 169594212480.0 * _XE_EULERGAMMA - 37468288063360.0 * _XE_LOG2 + 7085130274530.0 * _XE_LOG3 + 12068847656250.0 * _XE_LOG5 - 55474742400.0 * _XE_PI2) + 1925.0 * eta * (3400573970.0 + 49782159.0 * _XE_PI2)) + e2m1sqrt * (4435200.0 * chiA2 * (-2069461.0 + 9635809.0 * eta - 4395216.0 * eta2) + 4435200.0 * chiS2 * (-2069461.0 + 6612781.0 * eta - 2014516.0 * eta2) + 175.0 * e8 * (10539419949.0 + 3285893952.0 * eta + 1302605568.0 * eta2 + 178984960.0 * eta3) - 1490227200.0 * chiS * (-6148.0 + 4937.0 * eta) * pi + 8870400.0 * chiA * delta * (chiS * (-2069461.0 + 3985373.0 * eta) + 1032864.0 * pi) + 4.0 * e6 * (566190348111.0 + 9153699093000.0 * eta + 935550.0 * chiA2 * (-379573.0 + 1588880.0 * eta - 185472.0 * eta2) + 311850.0 * chiS2 * (-1138719.0 + 831208.0 * eta - 169792.0 * eta2) + 3117359044800.0 * eta2 + 392650720000.0 * eta3 + 122366946240.0 * _XE_EULERGAMMA - 584404081430400.0 * _XE_LOG2 + 159670846150200.0 * _XE_LOG3 + 144826171875000.0 * _XE_LOG5 + 207900.0 * chiA * delta * (3.0 * chiS * (-1138719.0 + 521486.0 * eta) - 618044.0 * pi) + 970200.0 * chiS * (-132438.0 + 42287.0 * eta) * pi + 3274425.0 * (-12224.0 + 6519.0 * eta) * _XE_PI2) + 64.0 * (-637780831131.0 + 380870611275.0 * eta + 30985883775.0 * eta2 + 12481353500.0 * eta3 + 109482468480.0 * _XE_EULERGAMMA + 102648712320.0 * _XE_LOG2 + 116761130640.0 * _XE_LOG3 - 727650.0 * (49216.0 + 17917.0 * eta) * _XE_PI2) + 160.0 * e2 * (-1094927702262.0 - 239784537300.0 * eta + 13860.0 * chiS2 * (7701538.0 + 6415619.0 * eta - 3646412.0 * eta2) + 275391017055.0 * eta2 - 13860.0 * chiA2 * (-7701538.0 + 28354633.0 * eta + 6528312.0 * eta2) + 37285931375.0 * eta3 + 158429246976.0 * _XE_EULERGAMMA + 1587595451904.0 * _XE_LOG2 - 572129540136.0 * _XE_LOG3 - 1164240.0 * chiS * (-39223.0 + 74954.0 * eta) * pi + 27720.0 * chiA * delta * (chiS * (7701538.0 + 4433569.0 * eta) + 1647366.0 * pi) + 291060.0 * (-178048.0 + 28413.0 * eta) * _XE_PI2) + 12.0 * e4 * (138600.0 * chiA2 * (-2004643.0 + 9123403.0 * eta - 1078672.0 * eta2) + 138600.0 * chiS2 * (-2004643.0 + 1416191.0 * eta - 1056188.0 * eta2) + 92400.0 * chiA * delta * (chiS * (-6013929.0 + 3781533.0 * eta) - 7879319.0 * pi) - 646800.0 * chiS * (1125617.0 + 95561.0 * eta) * pi + 6.0 * (-2199381667981.0 + 846453444375.0 * eta2 + 99510719000.0 * eta3 + 169594212480.0 * _XE_EULERGAMMA - 37468288063360.0 * _XE_LOG2 + 7085130274530.0 * _XE_LOG3 + 12068847656250.0 * _XE_LOG5 - 55474742400.0 * _XE_PI2) + 1925.0 * eta * (2808444530.0 + 54509049.0 * _XE_PI2)))) * x3) / (e2m13 * (1.0 + e2m1sqrt))
    edot = fact0 * (edot0PN + edot1PN + edot15PN + edot2PN + edot25PN + edot3PN)
    delta2 = delta * delta
    ldot0pn = 1.0
    ldot1pn = -3.0 / e2m1
    ldot2pn = (-18.0 + 28.0 * eta + e2 * (-51.0 + 26.0 * eta)) / (4.0 * e2m12)
    ldot15pn = (4.0 * chiS + 4.0 * chiA * delta - 2.0 * chiS * eta) / e2m132
    ldot2pn_ss = -3.0 * (chiS2 + 2.0 * chiA * chiS * delta + chiA2 * delta2) / (2.0 * e2m12)
    ldot25pn = (chiA * delta * (20.0 + e2 * (60.0 - 28.0 * eta) - 17.0 * eta) / 2.0 + chiS * (10.0 - 57.0 * eta / 2.0 + 2.0 * eta2 + e2 * (30.0 - 29.0 * eta + 7.0 * eta2))) / e2m152
    ldot3pn_ss = (2.0 * chiA * chiS * delta * (17.0 + e2 * (51.0 - 25.0 * eta) - 39.0 * eta) + chiS2 * (17.0 + e2 * (51.0 - 34.0 * eta) - 67.0 * eta + 20.0 * eta2) + chiA2 * (17.0 - 79.0 * eta + 20.0 * eta2 + e2 * (51.0 - 220.0 * eta + 58.0 * eta2))) / (2.0 * e2m13)
    ldot3pn = (16.0 * e4 * (156.0 + 5.0 * eta * (-22.0 + 13.0 * eta)) + e2 * (96.0 * (89.0 + 40.0 * e2m1sqrt) + eta * (-64.0 * (279.0 - 80.0 * eta + 24.0 * e2m1sqrt) + 123.0 * _XE_PI2)) + 4.0 * (-48.0 + 480.0 * e2m1sqrt + eta * (-8.0 * (457.0 - 28.0 * eta + 24.0 * e2m1sqrt) + 123.0 * _XE_PI2))) * x3 / (128.0 * e2m13)
    ldot = x32 * (ldot0pn + ldot1pn * x + x2 * (ldot2pn + ldot2pn_ss) + x32 * ldot15pn + x52 * ldot25pn + x3 * (ldot3pn + ldot3pn_ss))
    lamdot = x32
    return (xdot_diff, edot, ldot, lamdot)
end

function _et_from_eEOB(e, x, eta, delta, chiS, chiA, pn_order = 6.0)
    etasq = eta * eta
    etacu = etasq * eta
    chiSsq = chiS * chiS
    chiSchiA = chiS * chiA
    chiAsq = chiA * chiA
    esq = e * e
    omesq = 1.0 - esq
    sqrt_omesq = omesq^0.5
    ee_coeffs_0 = 1.0 * 1.0
    ee_coeffs_2 = eta - 3.0
    ee_coeffs_3 = 2.0 * (delta * chiA + chiS) / sqrt_omesq
    ee_coeffs_4 = (30.0 - 45.0 * eta + 4.0 * etasq + (-18.0 + 3.0 * eta - 4.0 * etasq) * esq + 9.0 * (2.0 * eta - 5.0) * sqrt_omesq) / (6.0 * omesq) - ((3.0 - 12.0 * eta) * chiAsq + 6.0 * delta * chiSchiA + 3.0 * chiSsq) / (2.0 * omesq)
    ee_coeffs_5 = -(delta * (96.0 + (-47.0 + 123.0 * esq) * eta) * chiA + (96.0 + (475.0 + 201.0 * esq) * eta + 2.0 * (5.0 - 9.0 * esq) * etasq) * chiS + (192.0 * delta * (eta - 3.0) * chiA - 96.0 * (6.0 - 8.0 * eta + etasq) * chiS) * sqrt_omesq) / (48.0 * omesq * sqrt_omesq)
    ee_coeffs_6 = (-26880.0 + 47144.0 * eta + 24640.0 * etasq + 2240.0 * etacu + 560.0 * (66.0 + 3.0 * eta + 120.0 * etasq - 8.0 * etacu) * esq + 280.0 * (-21.0 + 8.0 * (eta + etasq)) * eta * esq * esq - 35.0 * (1440.0 - (9824.0 - 123.0 * _XE_PI2 + 672.0 * esq) * eta + 96.0 * (8.0 + 7.0 * esq) * etasq) * sqrt_omesq) / (6720.0 * omesq * omesq) + ((-3.0 + 30.0 * eta - 6.0 * etasq + 3.0 * (9.0 - 7.0 * eta + 2.0 * etasq) * esq) * chiSsq - 6.0 * delta * (1.0 - 2.0 * eta + 3.0 * (eta - 3.0) * esq) * chiSchiA + 3.0 * (4.0 * eta - 1.0) * (1.0 + 6.0 * eta - (9.0 + eta) * esq) * chiAsq + ((-66.0 + 96.0 * eta - 24.0 * etasq) * chiSsq + 12.0 * delta * (-11.0 + 9.0 * eta) * chiSchiA - 6.0 * (11.0 - 46.0 * eta + 6.0 * etasq) * chiAsq) * sqrt_omesq) / (6.0 * omesq * omesq)
    x2 = x * x
    x_sqrt = sqrt(x)
    x15 = x * x_sqrt
    x25 = x2 * x_sqrt
    x3 = x2 * x
    if pn_order == 0
        et_sum = ee_coeffs_0
    elseif pn_order == 1
        et_sum = ee_coeffs_0
    elseif pn_order == 2
        et_sum = ee_coeffs_0 + x * ee_coeffs_2
    elseif pn_order == 3
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3
    elseif pn_order == 4
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4
    elseif pn_order == 5
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5
    elseif pn_order == 6
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5 + x3 * ee_coeffs_6
    else
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5 + x3 * ee_coeffs_6
    end
    et_HM = et_sum * e
    return et_HM
end

function _eEOB_from_et(e, x, eta, delta, chiS, chiA, pn_order = 6.0)
    etasq = eta * eta
    etacu = etasq * eta
    chiSsq = chiS * chiS
    chiSchiA = chiS * chiA
    chiAsq = chiA * chiA
    esq = e * e
    e4 = esq * esq
    delta2 = delta * delta
    omesq = 1.0 - esq
    sqrt_omesq = omesq^0.5
    omesq32 = omesq * sqrt_omesq
    omesq52 = omesq * omesq * sqrt_omesq
    ee_coeffs_0 = 1.0
    ee_coeffs_2 = 4.0 - eta
    ee_coeffs_3 = (-2.0 * chiS - 2.0 * chiA * delta + sqrt_omesq) / sqrt_omesq
    ee_coeffs_4 = (48.0 + 9.0 * chiSsq + 18.0 * chiA * chiS * delta - 60.0 * esq + chiAsq * (9.0 - 36.0 * eta) + 3.0 * eta + 39.0 * esq * eta + 2.0 * etasq - 2.0 * esq * etasq + 45.0 * sqrt_omesq - 18.0 * eta * sqrt_omesq) / (6.0 - 6.0 * esq)
    ee_coeffs_5 = 0.020833333333333332 * (-48.0 * (-4.0 + eta) * omesq32 + chiA * delta * (9.0 * esq * (32.0 + 3.0 * eta) - 576.0 * (1.0 + sqrt_omesq) + eta * (145.0 + 192.0 * sqrt_omesq)) + chiS * (-3.0 * esq * (-96.0 - 35.0 * eta + 6.0 * etasq) + etasq * (10.0 - 96.0 * sqrt_omesq) - 576.0 * (1.0 + sqrt_omesq) + eta * (667.0 + 768.0 * sqrt_omesq))) / omesq32
    ee_coeffs_6 = (403200.0 + 151200.0 * e4 + 13440.0 * chiS * (-1.0 + esq) - 554400.0 * esq - 585760.0 * eta - 87360.0 * e4 * eta + 673120.0 * esq * eta + 67200.0 * etasq - 3360.0 * e4 * etasq - 63840.0 * esq * etasq + 4305.0 * eta * _XE_PI2 - 4305.0 * esq * eta * _XE_PI2 + 60480.0 * sqrt_omesq + 127680.0 * e4 * sqrt_omesq - 305760.0 * esq * sqrt_omesq + 144376.0 * eta * sqrt_omesq - 158760.0 * e4 * eta * sqrt_omesq + 374640.0 * esq * eta * sqrt_omesq - 89600.0 * etasq * sqrt_omesq + 26880.0 * e4 * etasq * sqrt_omesq - 125440.0 * esq * etasq * sqrt_omesq + 3360.0 * chiAsq * (22.0 - 92.0 * eta + 12.0 * etasq + 22.0 * sqrt_omesq + 8.0 * delta2 * sqrt_omesq - 88.0 * eta * sqrt_omesq + esq * (-22.0 + 4.0 * etasq * (-3.0 + sqrt_omesq) - 6.0 * sqrt_omesq + 23.0 * eta * (4.0 + sqrt_omesq))) - 3360.0 * chiSsq * (-2.0 * (11.0 + 15.0 * sqrt_omesq - 8.0 * eta * (2.0 + sqrt_omesq) + etasq * (4.0 + sqrt_omesq)) + esq * (22.0 + 6.0 * sqrt_omesq + 2.0 * etasq * (4.0 + sqrt_omesq) - eta * (32.0 + 7.0 * sqrt_omesq))) + 6720.0 * chiA * delta * (2.0 * (-1.0 + esq) + chiS * (22.0 + 30.0 * sqrt_omesq - 2.0 * eta * (9.0 + 4.0 * sqrt_omesq) + esq * (-22.0 - 6.0 * sqrt_omesq + 3.0 * eta * (6.0 + sqrt_omesq))))) / (6720.0 * omesq52)
    x2 = x * x
    x_sqrt = sqrt(x)
    x15 = x * x_sqrt
    x25 = x2 * x_sqrt
    x3 = x2 * x
    if pn_order == 0
        et_sum = ee_coeffs_0
    elseif pn_order == 1
        et_sum = ee_coeffs_0
    elseif pn_order == 2
        et_sum = ee_coeffs_0 + x * ee_coeffs_2
    elseif pn_order == 3
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3
    elseif pn_order == 4
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4
    elseif pn_order == 5
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5
    elseif pn_order == 6
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5 + x3 * ee_coeffs_6
    else
        et_sum = ee_coeffs_0 + x * ee_coeffs_2 + x15 * ee_coeffs_3 + x2 * ee_coeffs_4 + x25 * ee_coeffs_5 + x3 * ee_coeffs_6
    end
    et_HM = et_sum * e
    return et_HM
end
function _compute_dldt_3pn(x, e, eta, delta, chiS, chiA)
    delta2 = delta * delta
    chiA2 = chiA * chiA
    chiS2 = chiS * chiS
    eta2 = eta * eta
    e2 = e * e
    e4 = e2 * e2
    _XE_PI2 = pi * pi
    e2m1 = 1.0 - e2
    e2m1_sqrt = e2m1^0.5
    e2m1_2 = e2m1 * e2m1
    e2m1_3 = e2m1_2 * e2m1
    e2m1_32 = e2m1 * e2m1_sqrt
    e2m1_52 = e2m1_2 * e2m1_sqrt
    x2 = x * x
    x3 = x2 * x
    x12 = x^0.5
    x32 = x * x12
    x52 = x2 * x12
    ldot0pn = 1.0
    ldot1pn = -3.0 / e2m1
    ldot15pn = 2.0 * (2.0 * chiA * delta - chiS * (-2.0 + eta)) / e2m1_32
    ldot2pn_ns = -18.0 + 28.0 * eta + e2 * (-51.0 + 26.0 * eta)
    ldot2pn_ss = -3.0 * (chiS2 + 2.0 * chiA * chiS * delta + chiA2 * delta2)
    ldot2pn = (ldot2pn_ns + ldot2pn_ss) / (4.0 * e2m1_2)
    ldot25pn = (chiA * delta * (20.0 + e2 * (60.0 - 28.0 * eta) - 17.0 * eta) + chiS * (20.0 - 57.0 * eta + 4.0 * eta2 + 2.0 * e2 * (30.0 - 29.0 * eta + 7.0 * eta2))) / (2.0 * e2m1_52)
    ldot3pn_ss = -0.5 * (2.0 * chiA * chiS * delta * (17.0 + e2 * (51.0 - 25.0 * eta) - 39.0 * eta) + chiS2 * (17.0 + e2 * (51.0 - 34.0 * eta) - 67.0 * eta + 20.0 * eta2) + chiA2 * (17.0 - 79.0 * eta + 20.0 * eta2 + e2 * (51.0 - 220.0 * eta + 58.0 * eta2))) / e2m1_3
    ldot3pn = (16.0 * e4 * (156.0 + 5.0 * eta * (-22.0 + 13.0 * eta)) + e2 * (96.0 * (89.0 + 40.0 * e2m1_sqrt) + eta * (-64.0 * (279.0 - 80.0 * eta + 24.0 * e2m1_sqrt) + 123.0 * _XE_PI2)) + 4.0 * (-48.0 + 480.0 * e2m1_sqrt + eta * (-8.0 * (457.0 - 28.0 * eta + 24.0 * e2m1_sqrt) + 123.0 * _XE_PI2))) * x3 / (128.0 * e2m1_3)
    dldt = x32 * (ldot0pn + ldot1pn * x + x2 * ldot2pn + x32 * ldot15pn + x52 * ldot25pn + x3 * (ldot3pn + ldot3pn_ss))
    return dldt
end

function _compute_dldtnorm_3pn(x, e, eta, delta, chiS, chiA)
    delta2 = delta * delta
    chiA2 = chiA * chiA
    chiS2 = chiS * chiS
    eta2 = eta * eta
    e2 = e * e
    e4 = e2 * e2
    _XE_PI2 = pi * pi
    e2m1 = 1.0 - e2
    e2m1_sqrt = e2m1^0.5
    e2m1_2 = e2m1 * e2m1
    e2m1_3 = e2m1_2 * e2m1
    e2m1_32 = e2m1 * e2m1_sqrt
    e2m1_52 = e2m1_2 * e2m1_sqrt
    x2 = x * x
    x3 = x2 * x
    x12 = x^0.5
    x32 = x * x12
    x52 = x2 * x12
    ldot0pn = 1.0
    ldot1pn = -3.0 / e2m1
    ldot15pn = 2.0 * (2.0 * chiA * delta - chiS * (-2.0 + eta)) / e2m1_32
    ldot2pn_ns = -18.0 + 28.0 * eta + e2 * (-51.0 + 26.0 * eta)
    ldot2pn_ss = -3.0 * (chiS2 + 2.0 * chiA * chiS * delta + chiA2 * delta2)
    ldot2pn = (ldot2pn_ns + ldot2pn_ss) / (4.0 * e2m1_2)
    ldot25pn = (chiA * delta * (20.0 + e2 * (60.0 - 28.0 * eta) - 17.0 * eta) + chiS * (20.0 - 57.0 * eta + 4.0 * eta2 + 2.0 * e2 * (30.0 - 29.0 * eta + 7.0 * eta2))) / (2.0 * e2m1_52)
    ldot3pn_ss = -0.5 * (2.0 * chiA * chiS * delta * (17.0 + e2 * (51.0 - 25.0 * eta) - 39.0 * eta) + chiS2 * (17.0 + e2 * (51.0 - 34.0 * eta) - 67.0 * eta + 20.0 * eta2) + chiA2 * (17.0 - 79.0 * eta + 20.0 * eta2 + e2 * (51.0 - 220.0 * eta + 58.0 * eta2))) / e2m1_3
    ldot3pn = (16.0 * e4 * (156.0 + 5.0 * eta * (-22.0 + 13.0 * eta)) + e2 * (96.0 * (89.0 + 40.0 * e2m1_sqrt) + eta * (-64.0 * (279.0 - 80.0 * eta + 24.0 * e2m1_sqrt) + 123.0 * _XE_PI2)) + 4.0 * (-48.0 + 480.0 * e2m1_sqrt + eta * (-8.0 * (457.0 - 28.0 * eta + 24.0 * e2m1_sqrt) + 123.0 * _XE_PI2))) * x3 / (128.0 * e2m1_3)
    dldt = ldot0pn + ldot1pn * x + x2 * ldot2pn + x32 * ldot15pn + x52 * ldot25pn + x3 * (ldot3pn + ldot3pn_ss)
    return dldt
end

##############################################################################
#   Generic helpers
##############################################################################

# Value of a (possibly dual) number; identity for ordinary reals
_val(x::Real) = ForwardDiff.value(x)

# sqrt(1 - 4 eta) protected against eta slightly above 1/4
_xe_delta(eta) = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), zero(eta))

# linspace that works with dual endpoints
function _linspace(a, b, n::Integer)
    n <= 1 && return [a]
    return [a + (b - a) * (k / (n - 1)) for k in 0:(n - 1)]
end

# Numpy-like gradient (second order interior, first order edges) on a non-uniform grid
function _np_gradient(f::AbstractVector, x::AbstractVector)
    n = length(f)
    g = similar(f, promote_type(eltype(f), eltype(x)))
    n < 2 && return g
    g[1] = (f[2] - f[1]) / (x[2] - x[1])
    g[n] = (f[n] - f[n - 1]) / (x[n] - x[n - 1])
    for i in 2:(n - 1)
        hs = x[i] - x[i - 1]
        hd = x[i + 1] - x[i]
        g[i] = (hs * hs * f[i + 1] + (hd * hd - hs * hs) * f[i] - hd * hd * f[i - 1]) / (hs * hd * (hd + hs))
    end
    return g
end

"""
Brent's root finder (as scipy.optimize.brentq) on the values of `g`, followed by one Newton correction performed
with the full (possibly dual) arithmetic. For Float64 inputs the correction is negligible; for ForwardDiff
dual numbers it propagates the derivative of the root with respect to the parameters (implicit function theorem).
"""
function _root_brent(g, a, b; xtol = 2e-12, rtol = 8.881784197001252e-16, maxiter = 100)
    gv = t -> _val(g(t))
    av = _val(a); bv = _val(b)
    fa = gv(av); fb = gv(bv)
    if fa * fb > 0
        error("PhenomXE: root is not bracketed in ($av, $bv)")
    end
    c = av; fc = fa
    d = bv - av; e = d
    xa = av; xb = bv
    for _ in 1:maxiter
        if fb * fc > 0
            c = xa; fc = fa; d = xb - xa; e = d
        end
        if abs(fc) < abs(fb)
            xa = xb; xb = c; c = xa
            fa = fb; fb = fc; fc = fa
        end
        tol = 0.5 * (xtol + rtol * abs(xb))
        xm = 0.5 * (c - xb)
        if abs(xm) <= tol || fb == 0.0
            break
        end
        if abs(e) >= tol && abs(fa) > abs(fb)
            s = fb / fa
            if xa == c
                p = 2.0 * xm * s; q = 1.0 - s
            else
                q = fa / fc; r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (xb - xa) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            end
            if p > 0
                q = -q
            end
            p = abs(p)
            if 2.0 * p < min(3.0 * xm * q - abs(tol * q), abs(e * q))
                e = d; d = p / q
            else
                d = xm; e = d
            end
        else
            d = xm; e = d
        end
        xa = xb; fa = fb
        if abs(d) > tol
            xb += d
        else
            xb += copysign(tol, xm)
        end
        fb = gv(xb)
    end
    t0 = xb
    # Newton correction with dual arithmetic (exact derivative propagation of the root)
    h = 1e-6 * max(1.0, abs(t0))
    dg = (gv(t0 + h) - gv(t0 - h)) / (2.0 * h)
    return t0 - g(t0) / dg
end

##############################################################################
#   Cubic splines (SciPy layout: coefficients c[4, n-1], segment polynomial
#   c1 u^3 + c2 u^2 + c3 u + c4 with u = x - x_i)
##############################################################################

"""
Not-a-knot cubic spline (identical to scipy.interpolate.CubicSpline default) for strictly increasing knots `x`.
Returns the coefficient matrix (4 x (n-1)) and the derivatives of the spline at the knots.
"""
function _notaknot_spline(x::AbstractVector, y::AbstractVector)
    n = length(x)
    (n < 3 || length(y) != n) && error("PhenomXE: not-a-knot spline inputs must have equal length >= 3")
    Tx = eltype(x)
    Ty = promote_type(eltype(x), eltype(y))
    dx = Vector{Tx}(undef, n - 1)
    slopes = Vector{Ty}(undef, n - 1)
    for i in 1:(n - 1)
        h = x[i + 1] - x[i]
        h <= 0 && error("PhenomXE: spline knots must be strictly increasing")
        dx[i] = h
        slopes[i] = (y[i + 1] - y[i]) / h
    end
    diagonal = Vector{Ty}(undef, n)
    rhs = Vector{Ty}(undef, n)
    if n == 3
        diagonal[1] = 1.0
        diagonal[2] = 2.0 * (dx[1] + dx[2])
        diagonal[3] = 1.0
        rhs[1] = 2.0 * slopes[1]
        rhs[2] = 3.0 * (dx[1] * slopes[2] + dx[2] * slopes[1])
        rhs[3] = 2.0 * slopes[2]
        upper_first = one(Tx)
        lower_last = one(Tx)
    else
        for i in 2:(n - 1)
            diagonal[i] = 2.0 * (dx[i - 1] + dx[i])
            rhs[i] = 3.0 * (dx[i] * slopes[i - 1] + dx[i - 1] * slopes[i])
        end
        d_start = x[3] - x[1]
        diagonal[1] = dx[2]
        rhs[1] = ((dx[1] + 2.0 * d_start) * dx[2] * slopes[1] + dx[1]^2 * slopes[2]) / d_start
        upper_first = d_start
        d_end = x[n] - x[n - 2]
        diagonal[n] = dx[n - 2]
        rhs[n] = (dx[n - 1]^2 * slopes[n - 2] + (2.0 * d_end + dx[n - 1]) * dx[n - 2] * slopes[n - 1]) / d_end
        lower_last = d_end
    end
    # Thomas algorithm: upper[i] is upper_first at the first row and dx[i-1] afterwards; lower[i-1] is dx[i] except on the last row
    for i in 2:n
        lower = i == n ? lower_last : dx[i]
        upper_previous = i == 2 ? upper_first : dx[i - 2]
        factor = lower / diagonal[i - 1]
        diagonal[i] -= factor * upper_previous
        rhs[i] -= factor * rhs[i - 1]
    end
    rhs[n] = rhs[n] / diagonal[n]
    for i in (n - 1):-1:1
        upper = i == 1 ? upper_first : dx[i - 1]
        rhs[i] = (rhs[i] - upper * rhs[i + 1]) / diagonal[i]
    end
    coeffs = Matrix{Ty}(undef, 4, n - 1)
    for i in 1:(n - 1)
        transition = (rhs[i] + rhs[i + 1] - 2.0 * slopes[i]) / dx[i]
        coeffs[1, i] = transition / dx[i]
        coeffs[2, i] = (slopes[i] - rhs[i]) / dx[i] - transition
        coeffs[3, i] = rhs[i]
        coeffs[4, i] = y[i]
    end
    return coeffs, rhs
end

"""
Shape-preserving PCHIP interpolant (scipy.interpolate.PchipInterpolator), same coefficient layout as `_notaknot_spline`.
"""
function _pchip_spline(x::AbstractVector, y::AbstractVector)
    n = length(x)
    Ty = promote_type(eltype(x), eltype(y))
    hk = [x[i + 1] - x[i] for i in 1:(n - 1)]
    mk = [(y[i + 1] - y[i]) / hk[i] for i in 1:(n - 1)]
    dk = Vector{Ty}(undef, n)
    if n == 2
        dk[1] = mk[1]; dk[2] = mk[1]
    else
        for i in 2:(n - 1)
            m0 = mk[i - 1]; m1 = mk[i]
            if sign(m0) != sign(m1) || m0 == 0 || m1 == 0
                dk[i] = zero(Ty)
            else
                w1 = 2.0 * hk[i] + hk[i - 1]
                w2 = hk[i] + 2.0 * hk[i - 1]
                whmean = (w1 / m0 + w2 / m1) / (w1 + w2)
                dk[i] = 1.0 / whmean
            end
        end
        dk[1] = _pchip_edge(hk[1], hk[2], mk[1], mk[2])
        dk[n] = _pchip_edge(hk[n - 1], hk[n - 2], mk[n - 1], mk[n - 2])
    end
    coeffs = Matrix{Ty}(undef, 4, n - 1)
    for i in 1:(n - 1)
        h = hk[i]
        coeffs[1, i] = (dk[i] + dk[i + 1] - 2.0 * mk[i]) / (h * h)
        coeffs[2, i] = (3.0 * mk[i] - 2.0 * dk[i] - dk[i + 1]) / h
        coeffs[3, i] = dk[i]
        coeffs[4, i] = y[i]
    end
    return coeffs
end

function _pchip_edge(h0, h1, m0, m1)
    d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if sign(d) != sign(m0)
        return zero(d)
    elseif sign(m0) != sign(m1) && abs(d) > 3.0 * abs(m0)
        return 3.0 * m0
    end
    return d
end

# Segment index (1-based) of a scalar query, as scipy's eval_cubic_spline: searchsorted(left) - 1, clamped
@inline function _spline_segment(x::AbstractVector, xq)
    i = searchsortedfirst(x, xq) - 1
    n = length(x)
    i < 1 && (i = 1)
    i > n - 1 && (i = n - 1)
    return i
end

@inline function _spline_eval_segment(c::AbstractMatrix, i::Integer, u)
    return ((c[1, i] * u + c[2, i]) * u + c[3, i]) * u + c[4, i]
end

@inline function _spline_deriv_segment(c::AbstractMatrix, i::Integer, u)
    return (3.0 * c[1, i] * u + 2.0 * c[2, i]) * u + c[3, i]
end

function _spline_eval(c::AbstractMatrix, x::AbstractVector, xq)
    i = _spline_segment(x, xq)
    return _spline_eval_segment(c, i, xq - x[i])
end

function _spline_deriv(c::AbstractMatrix, x::AbstractVector, xq)
    i = _spline_segment(x, xq)
    return _spline_deriv_segment(c, i, xq - x[i])
end

# Antiderivative of a piecewise cubic (zero at the first knot), evaluated at the knots and at an arbitrary point
function _spline_antiderivative_knots(c::AbstractMatrix, x::AbstractVector)
    n = length(x)
    F = Vector{promote_type(eltype(c), eltype(x))}(undef, n)
    F[1] = 0.0
    for i in 1:(n - 1)
        h = x[i + 1] - x[i]
        F[i + 1] = F[i] + ((c[1, i] * h / 4.0 + c[2, i] / 3.0) * h + c[3, i] / 2.0) * h * h + c[4, i] * h
    end
    return F
end

function _spline_antiderivative(c::AbstractMatrix, x::AbstractVector, F::AbstractVector, xq)
    i = _spline_segment(x, xq)
    u = xq - x[i]
    return F[i] + ((c[1, i] * u / 4.0 + c[2, i] / 3.0) * u + c[3, i] / 2.0) * u * u + c[4, i] * u
end

"""
Whether a SciPy-layout cubic spline is strictly positive on every segment (checked exactly at the end points and
at the roots of the derivative), see phenomte.utils_tehm.cubic_spline_is_positive.
"""
function _cubic_spline_is_positive(c::AbstractMatrix, breaks::AbstractVector)
    for i in 1:size(c, 2)
        c0 = _val(c[1, i]); c1 = _val(c[2, i]); c2 = _val(c[3, i]); c3 = _val(c[4, i])
        h = _val(breaks[i + 1]) - _val(breaks[i])
        if c3 <= 0.0 || c0 * h^3 + c1 * h^2 + c2 * h + c3 <= 0.0
            return false
        end
        a = 3.0 * c0; b = 2.0 * c1
        if a == 0.0
            b == 0.0 && continue
            roots = (-c2 / b, -c2 / b)
        else
            disc = b * b - 4.0 * a * c2
            disc < 0.0 && continue
            sq = sqrt(disc)
            roots = ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a))
        end
        for r in roots
            if 0.0 < r < h && c0 * r^3 + c1 * r^2 + c2 * r + c3 <= 0.0
                return false
            end
        end
    end
    return true
end

"""
Replace dx/dt on [x_transition - width, x_transition + width] by a quintic Hermite polynomial matching value,
first and second derivative at the window edges (phenomte.utils_tehm.smooth_dxdt_local). The patch is discarded
when it leaves the range of the data it replaces.
"""
function _smooth_dxdt_local(x::AbstractVector, dxdt::AbstractVector, x_transition, width)
    n = length(x)
    xa = x_transition - width
    xb = x_transition + width
    ia = clamp(searchsortedfirst(x, xa), 1, n)
    ib = clamp(searchsortedfirst(x, xb), 1, n)
    if ib <= ia + 2
        error("PhenomXE: dx/dt smoothing interval contains too few points. Increase width or refine sampling.")
    end
    d1 = _np_gradient(dxdt, x)
    d2 = _np_gradient(d1, x)
    h = x[ib] - x[ia]
    y0 = dxdt[ia]; y1 = dxdt[ib]
    p0 = d1[ia] * h; p1 = d1[ib] * h
    q0 = d2[ia] * h * h; q1 = d2[ib] * h * h
    patched = similar(dxdt, ib - ia + 1)
    for (k, idx) in enumerate(ia:ib)
        s = (x[idx] - x[ia]) / h
        s2 = s * s; s3 = s2 * s; s4 = s3 * s; s5 = s4 * s
        H0 = 1.0 - 10.0 * s3 + 15.0 * s4 - 6.0 * s5
        H1 = s - 6.0 * s3 + 8.0 * s4 - 3.0 * s5
        H2 = 0.5 * (s2 - 3.0 * s3 + 3.0 * s4 - s5)
        H3 = 10.0 * s3 - 15.0 * s4 + 6.0 * s5
        H4 = -4.0 * s3 + 7.0 * s4 - 3.0 * s5
        H5 = 0.5 * (s3 - 2.0 * s4 + s5)
        patched[k] = y0 * H0 + p0 * H1 + q0 * H2 + y1 * H3 + p1 * H4 + q1 * H5
    end
    lower = minimum(_val, view(dxdt, ia:ib))
    upper = maximum(_val, view(dxdt, ia:ib))
    pmin = minimum(_val, patched)
    pmax = maximum(_val, patched)
    out = copy(dxdt)
    if pmin <= 0.0 || pmin < lower || pmax > upper
        return out
    end
    out[ia:ib] .= patched
    return out
end

# Mask selecting a strictly increasing subsequence (running maximum), see phenomte.utils_tehm.strictly_increasing_mask
function _strictly_increasing_mask(x::AbstractVector)
    n = length(x)
    keep = falses(n)
    n == 0 && return keep
    keep[1] = true
    running_max = x[1]
    for i in 2:n
        if x[i] > running_max
            keep[i] = true
            running_max = x[i]
        end
    end
    return keep
end

# Clip dx/dt from below at `relative_floor` times its running maximum (phenomte.utils_tehm.enforce_positive_dxdt)
function _enforce_positive_dxdt(dxdt::AbstractVector; relative_floor = 2.0e-2)
    n = length(dxdt)
    n == 0 && return dxdt
    running_max = dxdt[1]
    needs_floor = false
    for i in 1:n
        if dxdt[i] > running_max
            running_max = dxdt[i]
        elseif dxdt[i] < relative_floor * running_max
            needs_floor = true
            break
        end
    end
    needs_floor || return dxdt
    floored = copy(dxdt)
    running_max = dxdt[1]
    for i in 1:n
        if dxdt[i] > running_max
            running_max = dxdt[i]
        end
        mn = relative_floor * running_max
        if floored[i] < mn
            floored[i] = mn
        end
    end
    return floored
end

##############################################################################
#   DOP853 integrator for the quasi-Keplerian dynamics
#   (port of phenomxpy.eccentricity.compiled_dop853; same step controller as SciPy)
##############################################################################

const _DOP853_MAX_FACTOR = 10.0
const _DOP853_MIN_FACTOR = 0.2
const _DOP853_SAFETY = 0.9
const _DOP853_ERROR_EXPONENT = -1.0 / 8.0

# Right-hand side: 3PN EOB quasi-Keplerian equations (without quasi-circular x-dot) + PhenomT dx/dt(x) interpolant
@inline function _xe_rhs(y::AbstractVector, eta, delta, chiS, chiA, sx::AbstractVector, sc::AbstractMatrix)
    x = y[1]; e = y[2]
    xdot_ecc, edot, ldot, lamdot = _rhs_3pn_eob_noqc(x, e, eta, delta, chiS, chiA)
    xdot = _spline_eval(sc, sx, x) + xdot_ecc
    return [xdot, edot, ldot, lamdot]
end

_rms_norm(v) = sqrt(sum(abs2, v) / length(v))

function _dop853_select_initial_step(y, t_bound, f, direction, rtol, atol, eta, delta, chiS, chiA, sx, sc)
    interval_length = abs(t_bound)
    scale = atol .+ abs.(y) .* rtol
    d0 = _rms_norm(y ./ scale)
    d1 = _rms_norm(f ./ scale)
    h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1
    h0 = min(h0, interval_length)
    y1 = y .+ h0 * direction .* f
    f1 = _xe_rhs(y1, eta, delta, chiS, chiA, sx, sc)
    d2 = _rms_norm((f1 .- f) ./ scale) / h0
    if d1 <= 1e-15 && d2 <= 1e-15
        h1 = max(1e-6, h0 * 1e-3)
    else
        h1 = (0.01 / max(d1, d2))^(1.0 / 8.0)
    end
    return min(100.0 * h0, h1, interval_length)
end

# One DOP853 step (12 stages) without error control; returns y_new, f_new and the stage matrix K (13 x 4)
function _dop853_step(y, f, h, eta, delta, chiS, chiA, sx, sc, K::AbstractMatrix)
    K[1, :] .= f
    for s in 2:12
        stage_y = copy(y)
        for j in 1:(s - 1)
            a = DOP853_A[s, j]
            a == 0.0 && continue
            for m in 1:4
                stage_y[m] += h * a * K[j, m]
            end
        end
        K[s, :] .= _xe_rhs(stage_y, eta, delta, chiS, chiA, sx, sc)
    end
    y_new = copy(y)
    for j in 1:12
        b = DOP853_B[j]
        for m in 1:4
            y_new[m] += h * b * K[j, m]
        end
    end
    f_new = _xe_rhs(y_new, eta, delta, chiS, chiA, sx, sc)
    K[13, :] .= f_new
    return y_new, f_new
end

function _dop853_error_norm(K::AbstractMatrix, y, y_new, h, rtol, atol)
    err5n2 = zero(eltype(y)); err3n2 = zero(eltype(y))
    for m in 1:4
        scale = atol + max(abs(y[m]), abs(y_new[m])) * rtol
        e5 = zero(eltype(y)); e3 = zero(eltype(y))
        for j in 1:13
            e5 += K[j, m] * DOP853_E5[j]
            e3 += K[j, m] * DOP853_E3[j]
        end
        err5n2 += (e5 / scale)^2
        err3n2 += (e3 / scale)^2
    end
    if err5n2 == 0.0 && err3n2 == 0.0
        return zero(eltype(y))
    end
    return abs(h) * err5n2 / sqrt((err5n2 + 0.01 * err3n2) * 4.0)
end

"""
Forward DOP853 integration from the reference point until x >= xf (the overshooting point is not stored) or until
e < 0, sqrt(x) >= 1, 1/x < 1, x decreases or e increases (phenomxpy integrate_forward_dop853).
"""
function _integrate_forward_dop853(y0::AbstractVector, xf, eta, delta, chiS, chiA, sx, sc, rtol, atol; max_points = 10000)
    T = eltype(y0)
    times = Vector{T}(undef, max_points)
    states = Matrix{T}(undef, 4, max_points)
    times[1] = 0.0
    states[:, 1] .= y0
    count = 1
    t = zero(T)
    y = copy(y0)
    f = _xe_rhs(y, eta, delta, chiS, chiA, sx, sc)
    h_abs = _dop853_select_initial_step(y, 1e25, f, 1.0, rtol, atol, eta, delta, chiS, chiA, sx, sc)
    K = Matrix{T}(undef, 13, 4)
    last_x = y[1]; last_e = y[2]
    while count < max_points
        step_rejected = false
        local y_new, f_new, h
        while true
            h = h_abs
            y_new, f_new = _dop853_step(y, f, h, eta, delta, chiS, chiA, sx, sc, K)
            error_norm = _dop853_error_norm(K, y, y_new, h, rtol, atol)
            if error_norm < 1.0
                factor = error_norm == 0.0 ? _DOP853_MAX_FACTOR : min(_DOP853_MAX_FACTOR, _DOP853_SAFETY * error_norm^_DOP853_ERROR_EXPONENT)
                if step_rejected
                    factor = min(1.0, factor)
                end
                h_abs *= factor
                break
            end
            h_abs *= max(_DOP853_MIN_FACTOR, _DOP853_SAFETY * error_norm^_DOP853_ERROR_EXPONENT)
            step_rejected = true
        end
        t += h
        y = y_new
        f = f_new
        x = y[1]; ecc = y[2]
        xdiff = x - last_x; ediff = ecc - last_e
        last_x = x; last_e = ecc
        if x >= xf
            break
        end
        count += 1
        times[count] = t
        states[:, count] .= y
        if ecc < 0.0 || sqrt(x) >= 1.0 || 1.0 / x < 1.0 || xdiff < 0.0 || ediff > 0.0
            break
        end
    end
    return times[1:count], states[:, 1:count]
end

"""
Backward DOP853 integration from the reference point to x = x_min. The last step is refined (by bisection on the
step fraction) so that the returned end point lies on x = x_min (phenomxpy integrate_backward_dop853).
"""
function _integrate_backward_dop853(y0::AbstractVector, x_min, e_max, eta, delta, chiS, chiA, sx, sc, rtol, atol; max_points = 10000)
    T = eltype(y0)
    times = Vector{T}(undef, max_points)
    states = Matrix{T}(undef, 4, max_points)
    times[1] = 0.0
    states[:, 1] .= y0
    count = 1
    t = zero(T)
    y = copy(y0)
    f = _xe_rhs(y, eta, delta, chiS, chiA, sx, sc)
    h_abs = _dop853_select_initial_step(y, -2e25, f, -1.0, rtol, atol, eta, delta, chiS, chiA, sx, sc)
    K = Matrix{T}(undef, 13, 4)
    Ktmp = Matrix{T}(undef, 13, 4)
    while count < max_points
        step_rejected = false
        local y_new, f_new, h
        while true
            h = -h_abs
            y_new, f_new = _dop853_step(y, f, h, eta, delta, chiS, chiA, sx, sc, K)
            error_norm = _dop853_error_norm(K, y, y_new, h, rtol, atol)
            if error_norm < 1.0
                factor = error_norm == 0.0 ? _DOP853_MAX_FACTOR : min(_DOP853_MAX_FACTOR, _DOP853_SAFETY * error_norm^_DOP853_ERROR_EXPONENT)
                if step_rejected
                    factor = min(1.0, factor)
                end
                h_abs *= factor
                break
            end
            h_abs *= max(_DOP853_MIN_FACTOR, _DOP853_SAFETY * error_norm^_DOP853_ERROR_EXPONENT)
            step_rejected = true
        end
        if y_new[1] <= x_min
            lower = 0.0; upper = 1.0
            terminal_y = y_new
            for _ in 1:12
                fraction = lower + (upper - lower) * _val((y[1] - x_min) / (y[1] - terminal_y[1]))
                if fraction <= lower || fraction >= upper
                    fraction = 0.5 * (lower + upper)
                end
                candidate_y, _ = _dop853_step(y, f, h * fraction, eta, delta, chiS, chiA, sx, sc, Ktmp)
                if candidate_y[1] > x_min
                    lower = fraction
                else
                    upper = fraction
                    terminal_y = candidate_y
                end
            end
            terminal_y, _ = _dop853_step(y, f, h * upper, eta, delta, chiS, chiA, sx, sc, Ktmp)
            count += 1
            times[count] = t + h * upper
            states[:, count] .= terminal_y
            break
        end
        t += h
        y = y_new
        f = f_new
        count += 1
        times[count] = t
        states[:, count] .= y
        if y[2] > e_max
            break
        end
    end
    return times[1:count], states[:, 1:count]
end

"""
Numerical evolution of the orbit-averaged quasi-Keplerian equations (3PN EOB + PhenomT quasi-circular x-dot)
from the reference point (x_ref, e_ref, l_ref, lam_ref) forward to x_peak and backward to x_min
(phenomxpy.eccentricity.eccentric_dynamics.eccentric_pn_evolve_PhenT with rhs_eqs = "eob_phenomT").
Returns the time array (starting at zero), the 4 x N state matrix (x, e, l, lambda) and the time of the reference point.
"""
function _eccentric_evolve_phenT(x_min, x_ref, e_ref, l_ref, lam_ref, eta, delta, chiA, chiS, sx, sc, x_peak; rtol = 1e-12, atol = 1e-12, e0_max = 0.9)
    if x_peak <= x_ref
        error("PhenomXE: reference frequency (x = $(_val(x_ref))) is higher than the maximum frequency (x = $(_val(x_peak)))")
    end
    T = promote_type(typeof(x_ref), typeof(e_ref), typeof(l_ref), typeof(lam_ref), typeof(eta), typeof(chiS), typeof(chiA), eltype(sc))
    z0 = T[x_ref, e_ref, l_ref, lam_ref]
    ft, fy = _integrate_forward_dop853(z0, x_peak, eta, delta, chiS, chiA, sx, sc, rtol, atol)
    if x_min < x_ref
        bt, by = _integrate_backward_dop853(z0, x_min, e0_max, eta, delta, chiS, chiA, sx, sc, rtol, atol)
        t_ref = -bt[end]
        nb = length(bt)
        combined_t = vcat(reverse(bt)[1:(nb - 1)], ft)
        combined_y = hcat(by[:, nb:-1:2], fy)
    else
        t_ref = ft[1]
        combined_t = copy(ft)
        combined_y = copy(fy)
    end
    combined_t .-= combined_t[1]
    # avoid repeated first / last points (checked on lambda, the last state component)
    if abs(combined_t[1] - combined_t[2]) < 1e-10 || combined_y[4, 2] - combined_y[4, 1] < 0
        combined_t = combined_t[2:end]
        combined_y = combined_y[:, 2:end]
    end
    if abs(combined_t[end] - combined_t[end - 1]) < 1e-10 || combined_y[4, end] - combined_y[4, end - 1] < 0
        combined_t = combined_t[1:(end - 1)]
        combined_y = combined_y[:, 1:(end - 1)]
    end
    return combined_t, combined_y, t_ref
end

##############################################################################
#   Quasi-Keplerian dynamics driven by the PhenomT (2,2) frequency
#   (phenomxpy.phenomxe.spa_xe_num.compute_qkp_dynamics_opt / compute_num_spa_j0_phase_opt)
##############################################################################

"""
Build the PhenomT inspiral-merger frequency grid, the dx/dt(x) interpolant and evolve the eccentric dynamics.
Returns a NamedTuple with the dynamics (tt, xt, et, lt, lamt, et_eob, omega, ndot, omega_dot), the reference time,
the PhenomT structure and the dx/dt(x) spline.
"""
function _compute_qkp_dynamics(eta, delta, chi1, chi2, chiA, chiS, x_start, x_ref, e0, l0, Mfmin;
                               rtol = 1e-12, atol = 1e-12, dt_fine = 0.5, dtheta_min = 0.002, nn_max = 10000, e_merger_max = 0.2)
    phT = _phenomT22_setup(eta, chi1, chi2, Mfmin; rtol = rtol, atol = atol)
    times_T, om22_T = _phT_omega22_grid(phT, dtheta_min, dt_fine, nn_max)

    x_inspiral_cut = (phT.omegaCut / 2.0)^(2.0 / 3.0)

    # dx/dt from the quasi-circular PhenomT frequency (compute_dxdt_no_max_check)
    omOrb = om22_T ./ 2.0
    x_T = omOrb .^ (2.0 / 3.0)
    _, domOrbdt = _notaknot_spline(times_T, omOrb)
    dxdt_T = 2.0 .* domOrbdt .* (2.0^(1.0 / 3.0)) ./ 3.0 .* om22_T .^ (-1.0 / 3.0)
    keep = _strictly_increasing_mask(x_T)
    if !all(keep)
        x_T = x_T[keep]
        dxdt_T = dxdt_T[keep]
    end
    dxdt_T = _enforce_positive_dxdt(dxdt_T)
    if times_T[1] < phT.tCut
        dxdt_T = _smooth_dxdt_local(x_T, dxdt_T, x_inspiral_cut, 0.01)
    end
    if length(x_T) < 3
        error("PhenomXE: the PhenomT frequency evolution left fewer than three usable samples of x(t); the dxdt(x) interpolant cannot be built.")
    end
    sc, _ = _notaknot_spline(x_T, dxdt_T)
    if !_cubic_spline_is_positive(sc, x_T)
        # the not-a-knot cubic undershoots to a negative dx/dt: use the monotone PCHIP interpolant instead
        sc = _pchip_spline(x_T, dxdt_T)
    end
    sx = x_T

    x_peak = x_T[end]
    x_peak_T = (phT.omegaPeak / 2.0)^(2.0 / 3.0)

    # evolve the quasi-Keplerian dynamics (lambda_ref = 0 by convention)
    tt, dyn, t_ref = _eccentric_evolve_phenT(x_start, x_ref, e0, l0, zero(e0), eta, delta, chiA, chiS, sx, sc, 1.2 * x_peak; rtol = rtol, atol = atol)
    xt = dyn[1, :]; et = dyn[2, :]; lt = dyn[3, :]; lamt = dyn[4, :]

    # keep only the part where x is non-decreasing
    idx_neg = findfirst(<(0), diff(xt))
    if idx_neg !== nothing
        tt = tt[1:idx_neg]; xt = xt[1:idx_neg]; et = et[1:idx_neg]; lt = lt[1:idx_neg]; lamt = lamt[1:idx_neg]
    end

    if et[end] > e_merger_max
        error("PhenomXE: eccentricity at merger e = $(_val(et[end])) is too high. Maximum allowed is $e_merger_max")
    end

    # go from EOB variables to PN ones (x is gauge invariant)
    et_eob = et
    et = [_et_from_eEOB(et_eob[i], xt[i], eta, delta, chiS, chiA, 6) for i in eachindex(et_eob)]
    if et[end] < 0
        for i in eachindex(et)
            if et[i] < 0
                et[i] = zero(et[i])
            end
        end
    end

    dldt_pn = [_compute_dldt_3pn(xt[i], et[i], eta, delta, chiS, chiA) for i in eachindex(xt)]
    omega = xt .^ 1.5
    c_dldt, ndot = _notaknot_spline(tt, dldt_pn)
    _, omega_dot = _notaknot_spline(tt, omega)

    # mean anomaly from the integrated PN mean motion, anchored at the reference point
    F = _spline_antiderivative_knots(c_dldt, tt)
    lt_ref = l0 - _spline_antiderivative(c_dldt, tt, F, t_ref)
    lt = lt_ref .+ F

    return (tt = tt, xt = xt, et = et, lt = lt, lamt = lamt, et_eob = et_eob, omega = omega, ndot = ndot, omega_dot = omega_dot,
            t_ref = t_ref, phT = phT, sx = sx, sc = sc, x_peak = x_peak, x_peak_T = x_peak_T)
end

"""
Numerical SPA of the j = 0 harmonic (compute_num_spa_j0_phase_opt): cut the dynamics at the PhenomT frequency peak and
return the SPA frequencies and phase together with the (cut) dynamics.
"""
function _compute_num_spa_j0_phase(dynT; mm = 2)
    tt = dynT.tt .- dynT.tt[end]
    xt = dynT.xt; et = dynT.et; lt = dynT.lt; lamt = dynT.lamt; et_eob = dynT.et_eob
    omega = dynT.omega; ndot = dynT.ndot; omega_dot = dynT.omega_dot
    if length(tt) < 4
        error("PhenomXE: the dynamics has only $(length(tt)) points, at least 4 are needed for interpolation")
    end
    x_peak_T = dynT.x_peak_T
    x_final = xt[end]
    if x_peak_T < x_final
        cx, _ = _notaknot_spline(tt, xt)
        t_peak = _root_brent(t -> _spline_eval(cx, tt, t) - x_peak_T, tt[1], tt[end])
    else
        t_peak = tt[end]
    end
    if x_final > x_peak_T
        # cut the dynamics at the PhenomT peak and append the interpolated peak point (cut_dynamics_PhT_peak)
        idx = findall(xt .<= x_peak_T)
        rows = (et, lt, lamt, et_eob, omega, ndot, omega_dot)
        vals = map(r -> (c = _notaknot_spline(tt, r)[1]; _spline_eval(c, tt, t_peak)), rows)
        tt = vcat(tt[idx], t_peak); xt = vcat(xt[idx], x_peak_T)
        et = vcat(et[idx], vals[1]); lt = vcat(lt[idx], vals[2]); lamt = vcat(lamt[idx], vals[3]); et_eob = vcat(et_eob[idx], vals[4])
        omega = vcat(omega[idx], vals[5]); ndot = vcat(ndot[idx], vals[6]); omega_dot = vcat(omega_dot[idx], vals[7])
    end
    v_spa_qc = sqrt.(xt)
    dxdt = (2.0 / 3.0) .* omega_dot ./ v_spa_qc
    tc = tt[end]
    Mfs_spa = xt .^ 1.5 ./ pi
    phase_j_spa = 2.0 .* pi .* Mfs_spa .* (tt .- tc) .- (mm .* lamt .+ pi / 4.0)
    return (Mfs_spa = Mfs_spa, phase_j_spa = phase_j_spa, tt = tt, xt = xt, et = et, lt = lt, lamt = lamt, omega = omega,
            ndot = ndot, omega_dot = omega_dot, dxdt = dxdt, v_spa_qc = v_spa_qc, tc = tc, t_peak = t_peak)
end
