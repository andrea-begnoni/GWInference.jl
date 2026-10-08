#! format: off
#######################################################
# IMRPhenomT (2,2) mode: parameter-space fits, frequency and amplitude ansaetze
#######################################################
#
# Port of phenomxpy/phenomt (fits.py, internals.py, numba_ansaetze.py) restricted to the (2,2) mode.
# Used by IMRPhenomXE (see PhenomXE.jl) to obtain the quasi-circular dx/dt(x) that drives the eccentric dynamics
# and the NR-calibrated inspiral-merger amplitude of the j = 0 harmonic.
#
# The parameter-space fits below (functions _IMRPhenomT_*) were transliterated automatically from
# phenomxpy/phenomt/fits.py (mode 22 branches) and checked numerically against the Python implementation.

function _IMRPhenomT_Inspiral_TaylorT3_t0(eta, S, dchi, delta)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    S5 = S * S4
    dchi2 = dchi * dchi
    fit = 1.0 / eta * ((-20.74399646637014 - 106.27711276502542 * eta) / (1.0 + 0.6516016033332481 * eta) + 0.0012450290074562259 * dchi * delta * (1.0 - 4701633.367918768 * eta) * eta2 - 111.5049997379579 * dchi * delta * (1.0 + 19.95458485773613 * eta) * S * eta2 + 1204.6829118499857 * (1.0 - 4.025474056585855 * eta) * dchi2 * eta3 + S * (338.7318821277009 - 1553.5891860091408 * eta + 19614.263378999745 * eta2 - 156449.78737303324 * eta3 + 577363.3090369126 * eta4 - 802867.433363341 * eta5) + (-55.75053935847546 - 290.36341163610575 * eta + 7873.7667183299345 * eta2 - 43585.59040070178 * eta3 + 87229.84668746481 * eta4 - 32469.263449695136 * eta5) * S2 + (-102.8269343111326 + 5121.845705262981 * eta - 93026.46878769135 * eta2 + 650989.6793529999 * eta3 - 1884606.1037110784 * eta4 + 1861602.620702142 * eta5) * S3 + (-7.294950933078567 + 314.24955197427136 * eta - 3751.8509582195657 * eta2 + 21205.339564205595 * eta3 - 46448.94771114493 * eta4 + 20310.512558558552 * eta5) * S4 + (97.22312282683716 - 4556.60375328623 * eta + 76308.73046927384 * eta2 - 468784.4188333802 * eta3 + 998692.0246600509 * eta4 - 322905.9042578296 * eta5) * S5)
    return fit
end

function _IMRPhenomT_Inspiral_Freq_CP(eta, S, dchi, delta, idx)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    S5 = S * S4
    S6 = S * S5
    dchi2 = dchi * dchi
    if idx == 1
        fit = -0.014968864336704284 * dchi * delta * (1.0 - 1.942061808318584 * eta) * eta2 + 0.0017312772309375462 * dchi * delta * (1.0 - 0.07106994121956058 * eta) * S * eta2 + S * (0.0019208448318368731 - 0.0013579968243452476 * eta - 0.0033501404728414627 * eta2 + 0.008914420175326192 * eta3) + 6.687615165457298e-06 * dchi2 * eta3 + (0.02104073275966069 + 717.1534194224539 * eta + 85.37320237350282 * eta2 + 12.789214868358362 * eta3 - 16.00243777208413 * eta4) / (1.0 + 32934.586638893634 * eta) + (-8.306810248117731e-06 + 9.918593182087119e-05 * eta - 0.003805916669791129 * eta2 + 0.009854209286892323 * eta3) * S2 + (-5.578836442449699e-06 - 0.0030378960591856616 * eta + 0.03746366675135751 * eta2 - 0.10298471015315146 * eta3) * S3 + (4.425141111368952e-05 - 0.0008702073302258368 * eta + 0.006538604805919268 * eta2 - 0.01578597166324495 * eta3) * S4 + (-1.9469656288570753e-05 + 0.002969863931498354 * eta - 0.03643271052162611 * eta2 + 0.09959495981802587 * eta3) * S5 + (-4.2037164406446896e-05 + 0.0007336074135429041 * eta - 0.005603356997202016 * eta2 + 0.013439843000090702 * eta3) * S6
    elseif idx == 2
        fit = -0.04486391236129559 * dchi * delta * (1.0 - 1.8997912248414794 * eta) * eta2 - 0.003531802135161727 * dchi * delta * (1.0 - 8.001211450141325 * eta) * S * eta2 + S * (0.0061664395419698285 - 0.0040934633081508905 * eta - 0.009180337242551828 * eta2 + 0.020338583755834694 * eta3) + 6.524644306613066e-05 * dchi2 * eta3 + 1.0 / (1.0 - 3.2125452791404148 * eta) * (0.03711511661217631 - 0.10663782888636487 * eta - 0.09963406984414182 * eta2 + 0.6597367702009397 * eta3 - 2.777344875144891 * eta4 + 4.220674345359693 * eta5) + (0.00044302547647888445 + 0.000424246501303979 * eta - 0.01394093576260671 * eta2 + 0.02634851560709597 * eta3) * S2 + (0.00011582043047950321 - 0.008282652950117982 * eta + 0.08965067576998058 * eta2 - 0.23963885130463913 * eta3) * S3 + (0.0006123158975881322 - 0.007809160444435783 * eta + 0.028517174579539676 * eta2 - 0.03717957419042746 * eta3) * S4 + (-8.85530893214531e-05 + 0.005939789043536808 * eta - 0.07106551435109858 * eta2 + 0.1891131957235774 * eta3) * S5 + (-0.0005110853374341054 + 0.0038762476596420855 * eta + 0.005094077179675256 * eta2 - 0.047971766995287136 * eta3) * S6
    elseif idx == 3
        fit = -0.10196878573773932 * dchi * delta * (1.0 - 1.8918584778973513 * eta) * eta2 - 0.018820536453940443 * dchi * delta * (1.0 - 3.7307154599131183 * eta) * S * eta2 - 0.00013162098437956188 * dchi2 * eta3 + S * (0.0145572994468378 - 0.0017482433991394227 * eta - 0.10299007619034371 * eta2 + 0.4581039376357615 * eta3 - 0.7123678787549022 * eta4) + (0.05489007025458171 + 5.852073438961151 * eta + 2.74597705533403 * eta2 + 4.834336623113389 * eta3 - 26.931994454691022 * eta4 + 57.67035368809743 * eta5) / (1.0 + 105.52132834236778 * eta) + (0.003001211395915229 + 0.0017929418998452987 * eta - 0.13776590125456148 * eta2 + 0.7471133710854526 * eta3 - 1.3620323111858437 * eta4) * S2 + (0.001143282743686261 - 0.05793457776296727 * eta + 0.7841331051705482 * eta2 - 3.4936244160305323 * eta3 + 4.802357041496856 * eta4) * S3 + (0.0009168588840889624 - 0.03261437094899735 * eta + 0.3472881896838799 * eta2 - 1.3634383958859384 * eta3 + 1.7313939586675267 * eta4) * S4 + (-0.0002794014744432316 + 0.055911057147527664 * eta - 0.8686311380514122 * eta2 + 4.096191294930781 * eta3 - 6.009676060669872 * eta4) * S5 + (-0.0005046018052528331 + 0.029804593053788925 * eta - 0.3792653361049425 * eta2 + 1.6366976231421981 * eta3 - 2.26904099961476 * eta4) * S6
    elseif idx == 4
        fit = -0.1831889759662071 * dchi * delta * (1.0 - 1.8484261527766557 * eta) * eta2 - 0.07586202965525136 * dchi * delta * (1.0 - 3.2918162656371983 * eta) * S * eta2 + 0.0019259052728265817 * dchi2 * eta3 + S * (0.02685637375751212 + 0.013341664908359861 * eta - 0.3057217933283597 * eta2 + 1.395763446325911 * eta3 - 2.2559396974665376 * eta4) + (0.0725639467287476 + 12.39400068457852 * eta + 12.907450928972402 * eta2 - 7.422660061864399 * eta3 + 66.32985901506036 * eta4 - 117.85875779454518 * eta5) / (1.0 + 168.63492460136445 * eta) + (0.0087781653701194 + 0.006944161553839352 * eta - 0.3301149078235105 * eta2 + 1.6835714783903248 * eta3 - 2.950404929598742 * eta4) * S2 + (0.0037229746496019625 - 0.17155338099487646 * eta + 2.5881802140836774 * eta2 - 13.14710199375518 * eta3 + 21.366803256010915 * eta4) * S3 + (0.00278507305662002 - 0.12475855143364532 * eta + 1.8640209516178643 * eta2 - 10.117078727717564 * eta3 + 17.94244821676711 * eta4) * S4 + (0.0010273954584773936 + 0.1713357629442166 * eta - 3.017249223460983 * eta2 + 15.855096360798678 * eta3 - 26.444621592311933 * eta4) * S5 + (-0.00012207946532225968 + 0.11709700788855186 * eta - 2.0950821618097026 * eta2 + 11.925324501640054 * eta3 - 21.683978511818076 * eta4) * S6
    elseif idx == 5
        fit = -0.2508206617297265 * dchi * delta * (1.0 - 1.861010982421798 * eta) * eta2 - 0.1392163711259171 * dchi * delta * (1.0 - 3.2669366465555796 * eta) * S * eta2 + 0.0023126403170013045 * dchi2 * eta3 + S * (0.036750064163293766 + 0.036904343404333906 * eta - 0.5238739410356437 * eta2 + 2.3292117112945223 * eta3 - 3.654184701923543 * eta4) + (0.08373610487663233 + 6.301736487754372 * eta + 9.03911386193751 * eta2 + 4.91153188278086 * eta3) / (1.0 + 72.64820846804257 * eta) + (0.014963449678540705 + 0.008354571522567225 * eta - 0.41723078020683 * eta2 + 2.2007932082378785 * eta3 - 4.245354787320365 * eta4) * S2 + (0.005706180633326235 - 0.15748500622007494 * eta + 2.3477109912232845 * eta2 - 11.413877195221694 * eta3 + 17.033120593116756 * eta4) * S3 + (0.003890296981717687 - 0.15985471334551038 * eta + 2.560312006077997 * eta2 - 14.400920672743332 * eta3 + 26.10406142567958 * eta4) * S4 + (0.005305988847210204 + 0.10869207132210629 * eta - 2.4201307115268875 * eta2 + 12.544899744864924 * eta3 - 19.550600837316903 * eta4) * S5 + (0.002917248769788225 + 0.11851143848720952 * eta - 2.6640023622893416 * eta2 + 15.993378498844761 * eta3 - 29.752144941054446 * eta4) * S6
    end
    return fit
end

function _IMRPhenomT_Intermediate_Freq_CP1(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    S5 = S * S4
    dchi2 = dchi * dchi
    fit = -0.3926039690467202 * dchi * delta * (1.0 - 2.359180951434749 * eta) * eta2 - 0.28551098014898896 * dchi * delta * (1.0 - 3.414696100901444 * eta) * S * eta2 + 0.003414004344822246 * dchi2 * eta3 + S * (0.05697014130854102 + 0.07170430925984912 * eta - 0.9606499306623374 * eta2 + 5.440955307244598 * eta3 - 10.594319036394571 * eta4) + (0.10030959768350425 + 44.56725135920024 * eta + 163.96290948585087 * eta2 - 143.05635831020462 * eta3 + 393.8084861740473 * eta4) / (1.0 + 436.6494065618 * eta) + (0.021213606590798472 + 0.2148355967310081 * eta - 2.7747405367196265 * eta2 + 13.771088220299802 * eta3 - 25.128755397215368 * eta4) * S2 + (-0.003645992092251503 + 0.2137524962844931 * eta - 0.644979226062801 * eta2 - 1.7314849842209137 * eta3 + 5.573297392347478 * eta4) * S3 + (0.029352214609533665 - 0.6020287633594307 * eta + 7.014738679280164 * eta2 - 36.027159248248296 * eta3 + 63.42605850359639 * eta4) * S4 + (0.0356519646654399 - 0.5569780178251297 * eta + 4.017784725334053 * eta2 - 15.05881246593488 * eta3 + 22.94821359434365 * eta4) * S5
    return fit
end

function _IMRPhenomT_PeakFrequency(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    S5 = S * S4
    dchi2 = dchi * dchi
    fit = 0.27212130745330404 + 0.40972689759932074 * eta - 0.0018392172960247433 * eta * dchi2 + S * (0.09558832959428547 - 0.04834585264918328 * eta - 0.15275173823699056 * eta2) - 3.4232387074402153 * eta2 + 32.853772442252605 * eta3 - 1.4976829186605336 * dchi * delta * (1.0 - 4.775645585721007 * eta) * eta3 - 0.9981117852179613 * dchi * delta * (1.0 - 5.260098925354571 * eta) * S * eta3 - 125.22505746137587 * eta4 + 179.3797198714914 * eta5 + (0.054391696704622204 - 0.1482682698299456 * eta + 0.08938162810617255 * eta2) * S2 + (-0.020719540055375383 + 0.5090144456500953 * eta - 1.5809441589349338 * eta2) * S3 + (0.024240736699062685 - 0.09490089674418004 * eta + 0.09518501714836035 * eta2) * S4 + (0.09759303647532228 - 1.105520690228567 * eta + 2.921271981239294 * eta2) * S5
    return fit
end

function _IMRPhenomT_RD_Freq_C2(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    S2 = S * S
    S3 = S * S2
    dchi2 = dchi * dchi
    fit = 0.1598180460429256 + 0.19120040104567676 * eta + (-0.012853620630980167 - 0.006532392920798404 * eta) * S - 0.7733759581766899 * eta2 + 0.18151402648790957 * dchi * delta * (1.0 - 9.041198282315879 * eta) * eta2 + 0.27147713896183995 * dchi * delta * (1.0 - 5.653323210961101 * eta) * S * eta2 - 0.01603489049446065 * dchi2 * eta3 + (-0.046785083372074494 + 0.102759380109996 * eta) * S2 + (0.0009883572415502464 - 0.050384608002279486 * eta) * S3
    return fit
end

function _IMRPhenomT_RD_Freq_C3(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    S2 = S * S
    S3 = S * S2
    dchi2 = dchi * dchi
    fit = 2.6456463496860927 - 28.079375863863458 * eta + 323.1691069138812 * eta2 - 0.5040057675360762 * dchi * delta * (1.0 + 21.786482297795278 * eta) * eta2 + 1.561247215701216 * dchi * delta * (1.0 - 1.7508069810164308 * eta) * S * eta2 + S * (3.091917073632116 - 17.345283345692266 * eta + 33.40735388809028 * eta2) - 1490.8128941604907 * eta3 + 0.1619056474567525 * dchi2 * eta3 + 2376.3257196613886 * eta4 + (0.734022429223849 - 0.029342234233198747 * eta - 9.281610698291932 * eta2) * S2
    return fit
end

function _IMRPhenomT_Inspiral_Amp_CP(eta, S, dchi, delta, mode, idx)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    eta6 = eta * eta5
    eta7 = eta * eta6
    eta8 = eta * eta7
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    dchi2 = dchi * dchi
    if idx == 1
        fit = 6.480771730217768e-05 * eta * dchi2 - 0.3543965558027252 * dchi * delta * (1.0 - 2.463526130684083 * eta) * eta3 + 0.01879295038873938 * dchi * delta * (1.0 - 5.236796607517272 * eta) * S * eta3 + S * (0.1472653807120573 * eta - 1.9636752493349356 * eta2 + 14.177521724634461 * eta3 - 48.94620901701877 * eta4 + 63.83730899015984 * eta5) + eta * (0.8493442097893826 - 13.211067914003836 * eta + 311.99021467938235 * eta2 - 4731.025904601601 * eta3 + 44821.93042533854 * eta4 - 264474.1374080295 * eta5 + 943246.2317701122 * eta6 - 1858813.5904328802 * eta7 + 1552477.8581809246 * eta8) + (0.04902976057622393 * eta - 1.0152511131279736 * eta2 + 8.286289152216145 * eta3 - 30.19775956110767 * eta4 + 40.670065442751955 * eta5) * S2 + (0.04780630695082567 * eta - 1.2177827888317065 * eta2 + 11.505675146308567 * eta3 - 46.733420749352135 * eta4 + 68.40821782168776 * eta5) * S3
    elseif idx == 2
        fit = 0.000100027278976821 * eta * dchi2 - 0.7578403155712378 * dchi * delta * (1.0 - 2.056456271350877 * eta) * eta3 - 0.14126282637778914 * dchi * delta * (1.0 - 2.5840771007494916 * eta) * S * eta3 + S * (0.2331970217833686 * eta - 1.5473968380422929 * eta2 + 5.973401506474942 * eta3 - 9.110484789161045 * eta4) + eta * (0.9904613241626621 - 6.708006572605403 * eta + 127.40270095439482 * eta2 - 1723.355339710798 * eta3 + 15430.10086310527 * eta4 - 88744.26044058547 * eta5 + 313650.01696201024 * eta6 - 617887.8122937253 * eta7 + 518220.9267888211 * eta8) + (0.08934817374146888 * eta - 0.8887847358339216 * eta2 + 3.7233864099350784 * eta3 - 5.814765403882651 * eta4) * S2 + (0.04471990627820145 * eta - 0.642458648615624 * eta2 + 3.393481171493086 * eta3 - 6.092083983738554 * eta4) * S3
    elseif idx == 3
        fit = 0.0002459376633671657 * eta * dchi2 - 0.8794763631110696 * dchi * delta * (1.0 - 2.0751630535350096 * eta) * eta3 - 0.3319387797134261 * dchi * delta * (1.0 - 3.1838055629892184 * eta) * S * eta3 + S * (0.23505507416274007 * eta - 1.2449030421324767 * eta2 + 4.315803728759738 * eta3 - 6.384257606413192 * eta4) + eta * (1.0208762064809185 - 3.3799457394243957 * eta + 16.242639717123314 * eta2 + 299.2297416582362 * eta3 - 5913.920743907752 * eta4 + 46388.231537995445 * eta5 - 192261.0498470111 * eta6 + 413750.14250475995 * eta7 - 364403.84935539874 * eta8) + (0.09630827896641526 * eta - 0.7915321134872877 * eta2 + 2.86907420250287 * eta3 - 4.038995403653199 * eta4) * S2 + (0.07395420485618898 * eta - 1.0289224187583748 * eta2 + 5.275845823734598 * eta3 - 9.206158044409037 * eta4) * S3
    end
    return fit
end

function _IMRPhenomT_Intermediate_Amp_CP1(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    eta6 = eta * eta5
    eta7 = eta * eta6
    eta8 = eta * eta7
    S2 = S * S
    S3 = S * S2
    dchi2 = dchi * dchi
    fit = 0.0004059354652663733 * eta * dchi2 - 0.9382383412276684 * dchi * delta * (1.0 - 2.509151362054917 * eta) * eta3 - 0.6560748977864668 * dchi * delta * (1.0 - 3.426294113321932 * eta) * S * eta3 + S * (0.23465398091766254 * eta - 1.3398914201113978 * eta2 + 5.9073801933446495 * eta3 - 10.84221896204708 * eta4) + eta * (1.2946032382158479 - 3.3343035556341816 * eta + 91.6430240976277 * eta2 - 1687.6195123629968 * eta3 + 19726.50907350641 * eta4 - 140798.18973779568 * eta5 + 594095.3303894227 * eta6 - 1358657.562562124 * eta7 + 1295891.2179017465 * eta8) + (0.03174875260265387 * eta + 0.23082150180902375 * eta2 - 1.9901867982613048 * eta3 + 4.009389679757772 * eta4) * S2 + (-0.04033221614773138 * eta + 0.8426888041517518 * eta2 - 4.742283264846479 * eta3 + 9.059923021547936 * eta4) * S3
    return fit
end

function _IMRPhenomT_PeakAmp(eta, S, dchi, delta, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    eta6 = eta * eta5
    S2 = S * S
    S3 = S * S2
    S4 = S * S3
    dchi2 = dchi * dchi
    fit = 0.0017885007700308166 * eta * dchi2 - 0.5846280668038513 * dchi * delta * (1.0 - 4.879882766464646 * eta) * eta3 - 0.874161608112943 * dchi * delta * (1.0 - 1.690095043235707 * eta) * S * eta3 + S * (0.203557188205307 * eta - 2.4368458739010563 * eta2 + 12.206344183078137 * eta3 - 23.417979354674692 * eta4) + eta * (1.4701266133411792 - 1.387711607537906 * eta + 25.641251409467607 * eta2 - 186.013359336165 * eta3 + 801.3039484150348 * eta4 - 1893.8181854645718 * eta5 + 1946.531703997353 * eta6) + (-0.0018659293826992745 * eta - 0.1888206507658455 * eta2 + 1.4677324802664107 * eta3 - 1.4019283350536489 * eta4) * S2 + (-0.14699838946027494 * eta + 2.6186847787143837 * eta2 - 15.574381075605208 * eta3 + 31.239292792717016 * eta4) * S3
    return fit
end

function _IMRPhenomT_Ringdown_Amp_C3(eta, S, dchi, mode)
    eta2 = eta * eta
    eta3 = eta * eta2
    S2 = S * S
    S3 = S * S2
    fit = -0.48053994718185694 + 0.7023672141561462 * eta + S * (-0.3597773028596323 + 1.4330280386796503 * eta - 3.239121799338561 * eta2) - 0.1993836305574211 * eta2 + (-0.2651107472061685 + 1.6433443489711386 * eta - 2.757772023954491 * eta2) * S2 + (-0.01973537883495192 - 0.2410762147438714 * eta + 2.7315015976869756 * eta2) * S3
    return fit
end

function _IMRPhenomT_fringfit(finalDimlessSpin, mode)
    x = finalDimlessSpin
    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2
    x5 = x3 * x2
    x6 = x3 * x3
    x7 = x4 * x3
    x8 = x6 * x2
    return_val = (0.05947169566573468 - 0.14989771215394762 * finalDimlessSpin + 0.09535606290986028 * x2 + 0.02260924869042963 * x3 - 0.02501704155363241 * x4 - 0.005852438240997211 * x5 + 0.0027489038393367993 * x6 + 0.0005821983163192694 * x7) / (1.0 - 2.8570126619966296 * finalDimlessSpin + 2.373335413978394 * x2 - 0.6036964688511505 * x4 + 0.0873798215084077 * x6)
    return return_val
end

function _IMRPhenomT_fdampfit(finalDimlessSpin, mode)
    x = finalDimlessSpin
    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2
    x5 = x3 * x2
    x6 = x3 * x3
    x7 = x4 * x3
    x8 = x4 * x4
    return_val = (0.014158792290965177 - 0.036989395871554566 * finalDimlessSpin + 0.026822526296575368 * x2 + 0.0008490933750566702 * x3 - 0.004843996907020524 * x4 - 0.00014745235759327472 * x5 + 0.0001504546201236794 * x6) / (1.0 - 2.5900842798681376 * finalDimlessSpin + 1.8952576220623967 * x2 - 0.31416610693042507 * x4 + 0.009002719412204133 * x6)
    return return_val
end

function _IMRPhenomT_fdampn2fit(finalDimlessSpin, mode)
    x = finalDimlessSpin
    x2 = finalDimlessSpin * finalDimlessSpin
    x3 = x2 * finalDimlessSpin
    x4 = x2 * x2
    x5 = x3 * x2
    x6 = x3 * x3
    x7 = x6 * x
    x8 = x7 * x
    return_val = 0.043611742588188715 + 1.0 / (2.0 - 1.9477781396815619 * x) * (-0.004016191313442792 * x - 0.0027646155943395426 * x2 + 0.001141927763953028 * x3 + 0.007938320030300492 * x4 - 0.0008263166671238823 * x5 - 0.014025760257115768 * x6 + 0.001792158578158245 * x7 + 0.008824138122361842 * x8)
    return return_val
end

##############################################################################
#   Final state fits (IMRPhenomX_FinalMass2017 / IMRPhenomX_FinalSpin2017, as used by PhenomT)
##############################################################################

function _IMRPhenomX_FinalMass2017(eta, s1, s2)
    delta = _xe_delta(eta)
    eta2 = eta * eta; eta3 = eta2 * eta; eta4 = eta3 * eta
    m1 = 0.5 * (1.0 + delta); m2 = 0.5 * (1.0 - delta)
    S = (m1 * m1 * s1 + m2 * m2 * s2) / (m1 * m1 + m2 * m2)
    S2 = S * S; S3 = S2 * S
    dchi = s1 - s2; dchi2 = dchi * dchi
    eqSpin = ((0.057190958417936644 * eta + 0.5609904135313374 * eta2 - 0.84667563764404 * eta3 + 3.145145224278187 * eta4) *
        (1.0 + (-0.13084389181783257 - 1.1387311580238488 * eta + 5.49074464410971 * eta2) * S +
         (-0.17762802148331427 + 2.176667900182948 * eta2) * S2 +
         (-0.6320191645391563 + 4.952698546796005 * eta - 10.023747993978121 * eta2) * S3)) /
        (1.0 + (-0.9919475346968611 + 0.367620218664352 * eta + 4.274567337924067 * eta2) * S)
    uneqSpin = -0.09803730445895877 * dchi * delta * (1.0 - 3.2283713377939134 * eta) * eta2 +
        0.01118530335431078 * dchi2 * eta3 -
        0.01978238971523653 * dchi * delta * (1.0 - 4.91667749015812 * eta) * eta * S
    return 1.0 - (eqSpin + uneqSpin)
end

function _IMRPhenomX_FinalSpin2017(eta, s1, s2)
    delta = _xe_delta(eta)
    m1 = 0.5 * (1.0 + delta); m2 = 0.5 * (1.0 - delta)
    m1Sq = m1 * m1; m2Sq = m2 * m2
    eta2 = eta * eta; eta3 = eta * eta2
    S = (m1Sq * s1 + m2Sq * s2) / (m1Sq + m2Sq)
    S2 = S * S; S3 = S * S2
    dchi = s1 - s2; dchi2 = dchi * dchi
    noSpin = (3.4641016151377544 * eta + 20.0830030082033 * eta2 - 12.333573402277912 * eta3) / (1.0 + 7.2388440419467335 * eta)
    eqSpin = (m1Sq + m2Sq) * S + (
        (-0.8561951310209386 * eta - 0.09939065676370885 * eta2 + 1.668810429851045 * eta3) * S +
        (0.5881660363307388 * eta - 2.149269067519131 * eta2 + 3.4768263932898678 * eta3) * S2 +
        (0.142443244743048 * eta - 0.9598353840147513 * eta2 + 1.9595643107593743 * eta3) * S3) /
        (1.0 + (-0.9142232693081653 + 2.3191363426522633 * eta - 9.710576749140989 * eta3) * S)
    uneqSpin = 0.3223660562764661 * dchi * delta * (1.0 + 9.332575956437443 * eta) * eta2 -
        0.059808322561702126 * dchi2 * eta3 +
        2.3170397514509933 * dchi * delta * (1.0 - 3.2624649875884852 * eta) * eta3 * S
    return noSpin + eqSpin + uneqSpin
end

##############################################################################
#   PhenomT (2,2) frequency ansaetze (phenomt/numba_ansaetze.py)
##############################################################################

# TaylorT3 omega as a function of theta = (-eta t / 5)^(-1/8); pn = [omega1PN, omega1halfPN, omega2PN, omega2halfPN, omega3PN, omega3halfPN]
function _phT_pn_omega(theta, pn)
    theta2 = theta * theta; theta3 = theta2 * theta; theta4 = theta2 * theta2
    theta5 = theta3 * theta2; theta6 = theta3 * theta3; theta7 = theta4 * theta3
    logterm = 107.0 * log(theta) / 280.0
    return theta3 / 4.0 * (1.0 + pn[1] * theta2 + pn[2] * theta3 + pn[3] * theta4 + pn[4] * theta5 + pn[5] * theta6 + logterm * theta6 + pn[6] * theta7)
end

# Inspiral omega: TaylorT3 augmented with 6 pseudo-PN terms (ps)
function _phT_inspiral_omega(t, eta, pn, ps)
    theta = (-eta * t / 5.0)^(-1.0 / 8.0)
    theta2 = theta * theta; theta3 = theta * theta2; theta4 = theta * theta3; theta5 = theta * theta4
    theta6 = theta * theta5; theta7 = theta * theta6; theta8 = theta4 * theta4; theta9 = theta8 * theta
    theta10 = theta9 * theta; theta11 = theta10 * theta; theta12 = theta11 * theta; theta13 = theta12 * theta
    logterm = 107.0 * log(theta) / 280.0
    fac = theta3 / 4.0
    taylort3 = 1.0 + pn[1] * theta2 + pn[2] * theta3 + pn[3] * theta4 + pn[4] * theta5 + pn[5] * theta6 + logterm * theta6 + pn[6] * theta7
    out = ps[1] * theta8 + ps[2] * theta9 + ps[3] * theta10 + ps[4] * theta11 + ps[5] * theta12 + ps[6] * theta13
    return fac * (taylort3 + out)
end

# Time derivative of the inspiral omega
function _phT_inspiral_domega(t, eta, pn, ps)
    theta = (-eta * t / 5.0)^(-1.0 / 8.0)
    theta2 = theta * theta; theta3 = theta * theta2; theta4 = theta * theta3; theta5 = theta * theta4
    theta6 = theta * theta5; theta7 = theta * theta6; theta8 = theta4 * theta4; theta9 = theta8 * theta
    theta10 = theta9 * theta; theta11 = theta10 * theta; theta12 = theta11 * theta; theta13 = theta12 * theta
    logterm = log(theta)
    der_omega = 0.25 * theta2 * (3.0 + 5.0 * pn[1] * theta2 + 6.0 * pn[2] * theta3 + 7.0 * pn[3] * theta4 + 8.0 * pn[4] * theta5 +
        107.0 / 280.0 * theta6 + 9.0 * pn[5] * theta6 + 10.0 * pn[6] * theta7 + 11.0 * ps[1] * theta8 + 12.0 * ps[2] * theta9 +
        13.0 * ps[3] * theta10 + 14.0 * ps[4] * theta11 + 15.0 * ps[5] * theta12 + 16.0 * ps[6] * theta13) +
        963.0 * theta8 * logterm / 1120.0
    der_theta = 0.125 * (5.0 / eta)^(1.0 / 8.0) * (-t)^(-9.0 / 8.0)
    return der_theta * der_omega
end

# Intermediate (merger) omega ansatz
function _phT_intermediate_omega(t, alpha1RD, omegaPeak, omegaRING, domegaPeak, C1, C2, C3)
    x = asinh(alpha1RD * t)
    w = 1.0 - omegaPeak / omegaRING + x * (domegaPeak / alpha1RD + x * (C1 + x * (C2 + x * C3)))
    return omegaRING * (1.0 - w)
end

# Time derivative of the intermediate omega ansatz
function _phT_intermediate_domega(t, alpha1RD, omegaRING, domegaPeak, C1, C2, C3)
    as = asinh(alpha1RD * t)
    return -omegaRING / sqrt(1.0 + (alpha1RD * t)^2) * (domegaPeak + alpha1RD * (2.0 * C1 * as + 3.0 * C2 * as * as + 4.0 * C3 * as^3))
end

# Ringdown omega ansatz and its derivative
function _phT_ringdown_omega(t, c1, c2, c3, c4, omegaRING)
    expC = exp(-c2 * t); expC2 = expC * expC
    num = -c1 * c2 * (2.0 * c4 * expC2 + c3 * expC)
    den = 1.0 + c4 * expC2 + c3 * expC
    return num / den + omegaRING
end

function _phT_ringdown_domega(t, c1, c2, c3, c4)
    expC = exp(c2 * t); expC2 = expC * expC
    num = c1 * c2 * c2 * expC * (4.0 * c4 * expC + c3 * (c4 + expC2))
    den = c4 + expC * (c3 + expC)
    return num / (den * den)
end

##############################################################################
#   PhenomT (2,2) amplitude ansaetze
##############################################################################

# Inspiral amplitude without the overall fac0 * x factor (amp_inspiral_ansatz_PhT in phenomxpy)
function _phT_inspiral_amp_bare(x, pr, pim, ps)
    xhalf = sqrt(x); x1half = x * xhalf; x2 = x * x; x2half = x2 * xhalf; x3 = x2 * x
    x3half = x3 * xhalf; x4 = x2 * x2; x4half = x4 * xhalf; x5 = x3 * x2
    ampreal = pr[1] + pr[2] * xhalf + pr[3] * x + pr[4] * x1half + pr[5] * x2 + pr[6] * x2half + pr[7] * x3 + pr[8] * x3half + pr[9] * log(16.0 * x) * x3
    ampimag = pim[1] * xhalf + pim[2] * x + pim[3] * x1half + pim[4] * x2 + pim[5] * x2half + pim[6] * x3 + pim[7] * x3half
    ampreal += ps[1] * x4 + ps[2] * x4half + ps[3] * x5
    return complex(ampreal, ampimag)
end

_phT_inspiral_amp(x, fac0, pr, pim, ps) = fac0 * x * _phT_inspiral_amp_bare(x, pr, pim, ps)

function _phT_intermediate_amp(t, alpha1RD, C1, C2, C3, C4, tshift)
    sech1 = 1.0 / cosh(alpha1RD * (t - tshift))
    sech2 = 1.0 / cosh(2.0 * alpha1RD * (t - tshift))
    return C1 + C2 * sech1 + C3 * sech2^(1.0 / 7.0) + C4 * (t - tshift) * (t - tshift)
end

##############################################################################
#   PhenomT (2,2) structure: coefficients of the omega and amplitude ansaetze
##############################################################################

# Piecewise IMR omega of the PhenomT 22 mode
function _phT_imr_omega(T, t)
    if t < T.tCut
        return _phT_inspiral_omega(t, T.eta, T.omega_pn, T.omega_ps)
    elseif t >= 0.0
        return _phT_ringdown_omega(t, T.c1, T.c2, T.c3, T.c4, T.omegaRING)
    else
        return _phT_intermediate_omega(t, T.alpha1RD, T.omegaPeak, T.omegaRING, T.domegaPeak, T.omegaMergerC1, T.omegaMergerC2, T.omegaMergerC3)
    end
end

# Derivative of the absolute value of the complex inspiral amplitude (pAmp._der_complex_amp_orientation with return_phase=False)
function _phT_der_abs_amp(t, eta, tCut, omega, pr, pim, ps, fac0, alpha1RD, omegaRING, domegaPeak, C1, C2, C3, omega_pn, omega_ps)
    x = (omega * 0.5)^(2.0 / 3.0)
    xhalf = sqrt(x); x1half = x * xhalf; x2 = x * x; x2half = x2 * xhalf; x3 = x2 * x
    x3half = x3 * xhalf; x4 = x2 * x2; x4half = x4 * xhalf; x5 = x3 * x2
    ampreal = pr[1] + pr[2] * xhalf + pr[3] * x + pr[4] * x1half + pr[5] * x2 + pr[6] * x2half + pr[7] * x3 + pr[8] * x3half +
        pr[9] * log(16.0 * x) * x3 + ps[1] * x4 + ps[2] * x4half + ps[3] * x5
    ampimag = pim[1] * xhalf + pim[2] * x + pim[3] * x1half + pim[4] * x2 + pim[5] * x2half + pim[6] * x3 + pim[7] * x3half
    dampreal = 0.5 * pr[2] / xhalf + pr[3] + 1.5 * pr[4] * xhalf + 2.0 * pr[5] * x + 2.5 * pr[6] * x1half + 3.0 * pr[7] * x2 +
        3.5 * pr[8] * x2half + pr[9] * x2 * (1.0 + 3.0 * log(16.0 * x)) + 4.0 * ps[1] * x3 + 4.5 * ps[2] * x3half + 5.0 * ps[3] * x4
    dampimag = 0.5 * pim[1] / xhalf + pim[2] + 1.5 * pim[3] * xhalf + 2.0 * pim[4] * x + 2.5 * pim[5] * x1half + 3.0 * pim[6] * x2 + 3.5 * pim[7] * x2half
    der_x_per_omega = cbrt(2.0 / omega) / 3.0
    if t < tCut
        der_omega_per_t = _phT_inspiral_domega(t, eta, omega_pn, omega_ps)
    else
        der_omega_per_t = _phT_intermediate_domega(t, alpha1RD, omegaRING, domegaPeak, C1, C2, C3)
    end
    amp = abs(complex(ampreal, ampimag))
    return fac0 * (ampreal * (dampreal * x + ampreal) + ampimag * (dampimag * x + ampimag)) / amp * der_x_per_omega * der_omega_per_t
end

"""
Set up the IMRPhenomT (2,2)-mode frequency and amplitude coefficients (phenomxpy pWF/pPhase/pAmp, mode 22),
for the aligned-spin binary with symmetric mass ratio `eta`, spins `chi1`, `chi2`, and starting geometric frequency `Mfmin`.
Returns a NamedTuple with every quantity needed by IMRPhenomXE.
"""
function _phenomT22_setup(eta, chi1, chi2, Mfmin; rtol = 1e-12, atol = 1e-12)
    T = promote_type(typeof(eta), typeof(chi1), typeof(chi2), typeof(Mfmin))
    delta = _xe_delta(eta)
    m1 = 0.5 * (1.0 + delta); m2 = 0.5 * (1.0 - delta)
    chiS = 0.5 * (chi1 + chi2); chiA = 0.5 * (chi1 - chi2)
    S = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)
    dchi = chi1 - chi2
    eta2 = eta * eta; eta3 = eta * eta2
    chi12 = chi1 * chi1; chi22 = chi2 * chi2; chi23 = chi2 * chi22

    afinal = _IMRPhenomX_FinalSpin2017(eta, chi1, chi2)
    Mfinal = _IMRPhenomX_FinalMass2017(eta, chi1, chi2)
    fring = _IMRPhenomT_fringfit(afinal, 22) / Mfinal
    fdamp = _IMRPhenomT_fdampfit(afinal, 22) / Mfinal
    fdampn2 = _IMRPhenomT_fdampn2fit(afinal, 22) / Mfinal

    # ---------------------------------------------------------------- pPhase (omega)
    tCut = -26.982976386771437 / eta   # inspiral-intermediate transition time of the 22 mode
    omega1PN = 743.0 / 2688.0 + (11.0 * eta) / 32.0
    omega1halfPN = (-19.0 * (chi1 + chi2) * eta) / 80.0 + (-113.0 * (-2.0 * chi1 * m1 - 2.0 * chi2 * m2) - 96.0 * pi) / 320.0
    omega2PN = ((56975.0 + 61236.0 * chi12 - 119448.0 * chi1 * chi2 + 61236.0 * chi22) * eta) / 258048.0 + (371.0 * eta2) / 2048.0 +
        (1855099.0 - 3429216.0 * chi12 * m1 - 3429216.0 * chi22 * m2) / 14450688.0
    omega2halfPN = (-17.0 * (chi1 + chi2) * eta2) / 128.0 + (-146597.0 * (-2.0 * chi1 * m1 - 2.0 * chi2 * m2) - 46374.0 * pi) / 129024.0 +
        (eta * (-2.0 * (chi1 * (1213.0 - 63.0 * delta) + chi2 * (1213.0 + 63.0 * delta)) + 117.0 * pi)) / 2304.0
    omega3PN = -720817631400877.0 / 288412611379200.0 - (16928263.0 * chi12) / 137625600.0 - (16928263.0 * chi22) / 137625600.0 -
        (16928263.0 * chi12 * delta) / 137625600.0 + (16928263.0 * chi22 * delta) / 137625600.0 +
        ((-2318475.0 + 18767224.0 * chi12 - 54663952.0 * chi1 * chi2 + 18767224.0 * chi22) * eta2) / 137625600.0 +
        (235925.0 * eta3) / 1769472.0 + (107.0 * MathConstants.eulergamma) / 280.0 - (6127.0 * chi1 * pi) / 12800.0 -
        (6127.0 * chi2 * pi) / 12800.0 - (6127.0 * chi1 * delta * pi) / 12800.0 + (6127.0 * chi2 * delta * pi) / 12800.0 +
        (eta * (632550449425.0 + 35200873512.0 * chi12 - 28527282000.0 * chi1 * chi2 + 9605339856.0 * chi12 * delta -
                1512.0 * chi22 * (-23281001.0 + 6352738.0 * delta) + 34172264448.0 * (chi1 + chi2) * pi - 22912243200.0 * pi^2)) / 104044953600.0 +
        (53.0 * pi^2) / 200.0 + (107.0 * log(2.0)) / 280.0
    omega3halfPN = (-12029.0 * (chi1 + chi2) * eta3) / 92160.0 +
        (eta2 * (507654.0 * chi1 * chi22 - 838782.0 * chi23 + chi2 * (-840149.0 + 507654.0 * chi12 - 870576.0 * delta) +
                 chi1 * (-840149.0 - 838782.0 * chi12 + 870576.0 * delta) + 1701228.0 * pi)) / 15482880.0 +
        (eta * (-1134.0 * chi23 * (-206917.0 + 71931.0 * delta) + chi1 * (-1496368361.0 - 429508815.0 * delta + 1134.0 * chi12 * (206917.0 + 71931.0 * delta)) -
                chi2 * (1496368361.0 - 429508815.0 * delta + 437064012.0 * chi12 * m1) - 437064012.0 * chi1 * chi22 * m2 -
                144.0 * (488825.0 + 923076.0 * chi12 - 1782648.0 * chi1 * chi2 + 923076.0 * chi22) * pi)) / 185794560.0 +
        (-2.0 * chi1 * (-6579635551.0 + 535759434.0 * chi12) * m1 + 13159271102.0 * chi2 * m2 - 1071518868.0 * chi23 * m2 +
         (-565550067.0 + 930460608.0 * chi12 * m1 + 930460608.0 * chi22 * m2) * pi) / 1300561920.0
    omega_pn = [omega1PN, omega1halfPN, omega2PN, omega2halfPN, omega3PN, omega3halfPN]

    # Inspiral collocation points in theta and solution of the pseudo-PN coefficients
    thetapoints = [0.33, 0.45, 0.55, 0.65, 0.75, 0.82]
    tt0 = _IMRPhenomT_Inspiral_TaylorT3_t0(eta, S, dchi, delta)
    tEarly = -5.0 / (eta * thetapoints[1]^8)
    thetaini = (eta * (tt0 - tEarly) / 5.0)^(-0.125)
    omegapoints = Vector{T}(undef, 6)
    omegapoints[1] = _phT_pn_omega(thetaini, omega_pn)
    for idx in 2:6
        omegapoints[idx] = _IMRPhenomT_Inspiral_Freq_CP(eta, S, dchi, delta, idx - 1)
    end
    Amat = Matrix{T}(undef, 6, 6)
    Bvec = Vector{T}(undef, 6)
    for i in 1:6
        th = thetapoints[i]
        Bvec[i] = 4.0 / (th * th * th) * (omegapoints[i] - _phT_pn_omega(th, omega_pn))
        thp = th^8
        for j in 1:6
            Amat[i, j] = thp
            thp *= th
        end
    end
    omega_ps = Amat \ Bvec

    # Ringdown coefficients
    omegaRING = 2.0 * pi * fring
    alpha1RD = 2.0 * pi * fdamp
    omegaPeak = _IMRPhenomT_PeakFrequency(eta, S, dchi, delta, 22)
    c2 = _IMRPhenomT_RD_Freq_C2(eta, S, dchi, delta, 22)
    c3 = _IMRPhenomT_RD_Freq_C3(eta, S, dchi, delta, 22)
    c4 = zero(T)
    c1 = (1.0 + c3 + c4) * (omegaRING - omegaPeak) / c2 / (c3 + 2.0 * c4)

    # Intermediate coefficients
    omegaCut = _phT_inspiral_omega(tCut, eta, omega_pn, omega_ps)
    tcpMerger = -5.0 / (eta * 0.95^8)
    omegaMergerCP = 1.0 - _IMRPhenomT_Intermediate_Freq_CP1(eta, S, dchi, delta, 22) / omegaRING
    omegaCutBar = 1.0 - omegaCut / omegaRING
    domegaCut = -_phT_inspiral_domega(tCut, eta, omega_pn, omega_ps) / omegaRING
    domegaPeak = -_phT_ringdown_domega(zero(T), c1, c2, c3, c4) / omegaRING
    ascut = asinh(alpha1RD * tCut); ascut2 = ascut * ascut; ascut3 = ascut * ascut2; ascut4 = ascut * ascut3
    bascut = asinh(alpha1RD * tcpMerger); bascut2 = bascut * bascut; bascut3 = bascut * bascut2; bascut4 = bascut * bascut3
    dencut = sqrt(1.0 + tCut * tCut * alpha1RD * alpha1RD)
    M3 = Matrix{T}(undef, 3, 3)
    B3 = Vector{T}(undef, 3)
    B3[1] = omegaCutBar - (1.0 - omegaPeak / omegaRING) - (domegaPeak / alpha1RD) * ascut
    M3[1, 1] = ascut2; M3[1, 2] = ascut3; M3[1, 3] = ascut4
    B3[2] = omegaMergerCP - (1.0 - omegaPeak / omegaRING) - (domegaPeak / alpha1RD) * bascut
    M3[2, 1] = bascut2; M3[2, 2] = bascut3; M3[2, 3] = bascut4
    B3[3] = domegaCut - domegaPeak / dencut
    M3[3, 1] = 2.0 * alpha1RD * ascut / dencut; M3[3, 2] = 3.0 * alpha1RD * ascut2 / dencut; M3[3, 3] = 4.0 * alpha1RD * ascut3 / dencut
    sol3 = M3 \ B3
    omegaMergerC1, omegaMergerC2, omegaMergerC3 = sol3[1], sol3[2], sol3[3]

    Tph = (eta = eta, tCut = tCut, omega_pn = omega_pn, omega_ps = omega_ps, omegaRING = omegaRING, alpha1RD = alpha1RD,
           omegaPeak = omegaPeak, domegaPeak = domegaPeak, omegaMergerC1 = omegaMergerC1, omegaMergerC2 = omegaMergerC2,
           omegaMergerC3 = omegaMergerC3, c1 = c1, c2 = c2, c3 = c3, c4 = c4)

    # Starting time of the waveform: time at which the 22 frequency equals 2 pi Mfmin (pPhase._get_time_of_freq)
    v = cbrt(pi * Mfmin)
    tlow_pn = -5.0 / (256.0 * eta * v^8) * 1.2
    tlow_test = ifelse(tlow_pn < -1e9, tlow_pn, -1e9)
    tEnd = 500.0
    time_of_freq = t -> 2.0 * pi * Mfmin - (t < tEarly ? _phT_imr_omega(Tph, t - tt0) : _phT_imr_omega(Tph, t))
    tmin = _root_brent(time_of_freq, tlow_test, tEnd; xtol = atol, rtol = rtol, maxiter = 1000)
    if _val(tmin) >= 0.0
        error("PhenomXE: PhenomT tmin = $(_val(tmin)) > 0 is invalid. Try lowering f_min.")
    end

    # ---------------------------------------------------------------- pAmp (22 mode)
    amp_inspiral_cut = -150.0
    tshift = zero(T)   # IMRPhenomT_tshift is zero for the 22 mode
    fac0 = 2.0 * eta * sqrt(16.0 * pi / 5.0)
    S0 = m1 * chi1 + m2 * chi2
    ampN = 1.0
    amp1PNreal = -107.0 / 42.0 + (55.0 * eta) / 42.0
    amp1halfPNreal = (-4.0 * chiS) / 3.0 - (4.0 * chiA * delta) / 3.0 + (4.0 * chiS * eta) / 3.0 + 2.0 * pi
    amp2PNreal = -2173.0 / 1512.0 - (1069.0 * eta) / 216.0 + (2047.0 * eta2) / 1512.0 + S0^2
    amp2halfPNreal = (-107.0 * pi) / 21.0 + (34.0 * eta * pi) / 21.0
    amp2halfPNimag = -24.0 * eta
    amp3PNreal = 27027409.0 / 646800.0 - (278185.0 * eta) / 33264.0 - (20261.0 * eta2) / 2772.0 + (114635.0 * eta3) / 99792.0 -
        (856.0 * MathConstants.eulergamma) / 105.0 + (2.0 * pi^2) / 3.0 + (41.0 * eta * pi^2) / 96.0
    amp3PNimag = (428.0 * pi) / 105.0
    amp3halfPNreal = (-2173.0 * pi) / 756.0 - (2495.0 * eta * pi) / 378.0 + (40.0 * eta2 * pi) / 27.0
    amp3halfPNimag = (14333.0 * eta) / 162.0 - (4066.0 * eta2) / 945.0
    amplog = -428.0 / 105.0
    pn_real = T[ampN, 0.0, amp1PNreal, amp1halfPNreal, amp2PNreal, amp2halfPNreal, amp3PNreal, amp3halfPNreal, amplog]
    pn_imag = T[0.0, 0.0, 0.0, 0.0, amp2halfPNimag, amp3PNimag, amp3halfPNimag]

    # Inspiral collocation points (t = -2000, -250, -150 M) and pseudo-PN amplitude coefficients
    tinsp = [-2000.0, -250.0, -150.0]
    zeros3 = zeros(T, 3)
    Ma = Matrix{T}(undef, 3, 3)
    Ba = Vector{T}(undef, 3)
    for i in 1:3
        om = _phT_imr_omega(Tph, tinsp[i])
        xx = (0.5 * om)^(2.0 / 3.0)
        xxhalf = sqrt(xx)
        ampoffset = real(_phT_inspiral_amp(xx, fac0, pn_real, pn_imag, zeros3))
        Ba[i] = (1.0 / fac0 / xx) * (_IMRPhenomT_Inspiral_Amp_CP(eta, S, dchi, delta, 22, i) - ampoffset)
        xp = xx * xx * xx * xx
        for j in 1:3
            Ma[i, j] = xp
            xp *= xxhalf
        end
    end
    pseudo = Ma \ Ba

    # Ringdown amplitude coefficients
    alpha2RD = 2.0 * pi * fdampn2
    alpha21RD = 0.5 * (alpha2RD - alpha1RD)
    ampPeak = _IMRPhenomT_PeakAmp(eta, S, dchi, delta, 22)
    c3a = _IMRPhenomT_Ringdown_Amp_C3(eta, S, dchi, 22)
    c2a = alpha21RD
    coshc3 = cosh(c3a); tanhc3 = tanh(c3a)
    if c2a > abs(0.5 * alpha1RD / tanhc3)
        c2a = -0.5 * alpha1RD / tanhc3
    end
    c1a = ampPeak * alpha1RD * coshc3 * coshc3 / c2a
    c4a = ampPeak - c1a * tanhc3

    # Intermediate amplitude coefficients (quasi-circular branch)
    om_cut = _phT_imr_omega(Tph, amp_inspiral_cut)
    x_cut = (om_cut * 0.5)^(2.0 / 3.0)
    ampinsp_cplx = _phT_inspiral_amp(x_cut, fac0, pn_real, pn_imag, pseudo)
    sign_ampinsp = _val(real(ampinsp_cplx)) >= 0.0 ? 1.0 : -1.0   # sign of the real part (value-based, ForwardDiff friendly)
    ampinsp = sign_ampinsp * abs(ampinsp_cplx)
    phi = alpha1RD * (amp_inspiral_cut - tshift); phi2 = 2.0 * phi
    sech1 = 1.0 / cosh(phi); sech2 = 1.0 / cosh(phi2)
    ampMergerCP1 = _IMRPhenomT_Intermediate_Amp_CP1(eta, S, dchi, delta, 22)
    tcpMergerA = -25.0
    phib = alpha1RD * (tcpMergerA - tshift)
    sech1b = 1.0 / cosh(phib); sech2b = 1.0 / cosh(2.0 * phib)
    dampMECO = sign_ampinsp * _phT_der_abs_amp(amp_inspiral_cut, eta, tCut, om_cut, pn_real, pn_imag, pseudo, fac0,
        alpha1RD, omegaRING, domegaPeak, omegaMergerC1, omegaMergerC2, omegaMergerC3, omega_pn, omega_ps)
    Mi = Matrix{T}(undef, 4, 4)
    Bi = Vector{T}(undef, 4)
    Mi[1, 1] = 1.0; Mi[1, 2] = sech1; Mi[1, 3] = sech2^(1.0 / 7.0); Mi[1, 4] = (amp_inspiral_cut - tshift)^2; Bi[1] = ampinsp
    Mi[2, 1] = 1.0; Mi[2, 2] = sech1b; Mi[2, 3] = sech2b^(1.0 / 7.0); Mi[2, 4] = (tcpMergerA - tshift)^2; Bi[2] = ampMergerCP1
    Mi[3, 1] = 1.0; Mi[3, 2] = 1.0; Mi[3, 3] = 1.0; Mi[3, 4] = 0.0; Bi[3] = ampPeak
    Mi[4, 1] = 0.0; Mi[4, 2] = -alpha1RD * sech1 * tanh(phi); Mi[4, 3] = (-2.0 / 7.0) * alpha1RD * sinh(phi2) * sech2^(8.0 / 7.0)
    Mi[4, 4] = 2.0 * (amp_inspiral_cut - tshift); Bi[4] = dampMECO
    soli = Mi \ Bi

    return (eta = eta, delta = delta, chiS = chiS, chiA = chiA, afinal = afinal, Mfinal = Mfinal, fring = fring, fdamp = fdamp,
            tCut = tCut, tt0 = tt0, tEarly = tEarly, omega_pn = omega_pn, omega_ps = omega_ps, omegaRING = omegaRING,
            alpha1RD = alpha1RD, omegaPeak = omegaPeak, domegaPeak = domegaPeak, omegaMergerC1 = omegaMergerC1,
            omegaMergerC2 = omegaMergerC2, omegaMergerC3 = omegaMergerC3, c1 = c1, c2 = c2, c3 = c3, c4 = c4, omegaCut = omegaCut,
            tmin = tmin, amp_inspiral_cut = amp_inspiral_cut, tshift = tshift, fac0 = fac0, pn_real = pn_real, pn_imag = pn_imag,
            pseudo = pseudo, ampPeak = ampPeak, mergerC1 = soli[1], mergerC2 = soli[2], mergerC3 = soli[3], mergerC4 = soli[4])
end

"""
Time grid and PhenomT (2,2) frequency in the inspiral and intermediate regions
(phenomte.utils_tehm.compute_omega22_inspint_PhenomT_adaptive_grids_v1): uniform grid in the TaylorT3 variable theta
during the inspiral and a fine grid with step `dt_fine` in the intermediate region up to the peak (t = 0).
"""
function _phT_omega22_grid(T, dtheta, dt_fine, nn_max)
    eta = T.eta
    tmin = T.tmin
    tCut = T.tCut
    if tmin < tCut
        theta_min = (-eta * tmin / 5.0)^(-1.0 / 8.0)
        theta_max = (-eta * tCut / 5.0)^(-1.0 / 8.0)
        npts = Int(floor(abs(_val(theta_max) - _val(theta_min)) / dtheta))
        npts = min(npts, nn_max)
        theta_grid = _linspace(theta_min, theta_max, npts)
        t_insp = (-5.0 / eta) .* theta_grid .^ (-8.0)
        # fine grid from the peak (t = 0) down to just above tCut (create_fine_time_grid)
        nfine = ceil(Int, (_val(tCut) + dt_fine) / (-dt_fine))
        t_int = [-dt_fine * k for k in (nfine - 1):-1:0]
        while tCut >= t_int[1]
            popfirst!(t_int)
        end
        om_insp = [_phT_inspiral_omega(t, eta, T.omega_pn, T.omega_ps) for t in t_insp]
        om_int = [_phT_intermediate_omega(t, T.alpha1RD, T.omegaPeak, T.omegaRING, T.domegaPeak, T.omegaMergerC1, T.omegaMergerC2, T.omegaMergerC3) for t in t_int]
        return vcat(t_insp, t_int), vcat(om_insp, om_int)
    else
        n_pts = abs(Int(floor(_val(tCut) / dt_fine)))
        t = reverse(_linspace(zero(tmin), tmin, n_pts))
        om = [_phT_intermediate_omega(ti, T.alpha1RD, T.omegaPeak, T.omegaRING, T.domegaPeak, T.omegaMergerC1, T.omegaMergerC2, T.omegaMergerC3) for ti in t]
        return t, om
    end
end
