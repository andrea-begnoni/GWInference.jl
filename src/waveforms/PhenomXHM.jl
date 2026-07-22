#! format: off

"""
ToDo: Need documentation
"""
function PolAbs(model::PhenomXHM,
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
function Pol(model::PhenomXHM,
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
ToDo: Need documentation
"""
function Phi(model::PhenomXHM,
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

# *************************************************
#                                                  
#          Amplitude Cutting Frequencies           
#                                                  
# *************************************************

#= Inspiral cutting frequency for the Amplitude =#
function _fcutInsp(model::PhenomXHM, eta, chi1, emm, fMECOlm)

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
  #Return the end frequency of the inspiral region and the beginning of the intermediate for the amplitude of one mode.



  
    if q < 20.
        fcut = fMECOlm;
    else
        #fcutEMR is the cutting frequency for extreme mass ratios that is given by a fit to the frequency of a particular geometrical structure of the amplitude
        fcutEMR = 1.25*emm*((0.011671068725758493 - 0.0000858396080377194*chi1 + 0.000316707064291237*chi1^2)*(0.8447212540381764 + 6.2873167352395125*eta))/(1.2857082764038923 - 0.9977728883419751*chi1);

        transition_eta = 0.0192234; # q=50
        sharpness = 0.004;
        funcs = 0.5 + 0.5 * tanh((eta - transition_eta)/sharpness);
        fcut = funcs * fMECOlm + (1 - funcs) * fcutEMR;
    end
    return fcut

end

#= Ringdown cutting frequency for the amplitude =#
function _fcutRD(model::PhenomXHM, fring, fdamp, ModeMixingOn, fring_22, fdamp_22)

  #Returns the end of the intermediate region and the beginning of the ringdown for the amplitude of one mode

    if ModeMixingOn == true
        fcut = fring_22 - 0.5 * fdamp_22; #v8
    else
        fcut = fring - fdamp; #v2
    end

    return fcut;
end




function _spinWeightedSphericalHarmonic(model::PhenomXHM,theta, modes, m)
    # Taken from arXiv:0709.0093v3 eq. (II.7), (II.8) and LALSimulation for the s=-2 case and up to l=4.
    # We assume already phi=0 and s=-2 to simplify the function

    if m >0
        Ylm = @. ifelse(
            modes == 21,
            sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 + cos(theta)),
            ifelse(
                modes == 22,
                sqrt(5.0 / (64.0 * pi)) * (1.0 + cos(theta)) * (1.0 + cos(theta)),
                ifelse(
                    modes == 32,
                    sqrt(7.0 / pi) * (cos(theta * 0.5)^4) * (-2.0 + 3.0 * cos(theta)) * 0.5,
                    ifelse(
                        modes == 33,
                        -sqrt(21.0 / (2.0 * pi)) * (cos(theta * 0.5)^5) * sin(theta * 0.5),
                        ifelse(
                            modes == 43,
                            -3.0 *
                            sqrt(7.0 / (2.0 * pi)) *
                            (cos(theta * 0.5)^5) *
                            (-1.0 + 2.0 * cos(theta)) *
                            sin(theta * 0.5),
                            3.0 * sqrt(7.0 / pi) * (cos(theta * 0.5)^6) * (sin(theta * 0.5)^2),
                        ),
                    ),
                ),
            ),
        )
    elseif m <0
        Ylm = @. ifelse(
            modes == 21,
            sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 - cos(theta)),
            ifelse(
                modes == 22,
                sqrt(5.0 / (64.0 * pi)) * (1.0 - cos(theta)) * (1.0 - cos(theta)),
                ifelse(
                    modes == 32,
                    sqrt(7.0 / (4.0 * pi)) *
                    (2.0 + 3.0 * cos(theta)) *
                    ((sin(theta * 0.5))^(4.0)),
                    ifelse(
                        modes == 33,
                        sqrt(21.0 / (2.0 * pi)) * cos(theta * 0.5) * ((sin(theta * 0.5))^(5.0)),
                        ifelse(
                            modes == 43,
                            3.0 *
                            sqrt(7.0 / (2.0 * pi)) *
                            cos(theta * 0.5) *
                            (1.0 + 2.0 * cos(theta)) *
                            ((sin(theta * 0.5))^5.0),
                            3.0 * sqrt(7.0 / pi) * (cos(theta * 0.5) * cos(theta * 0.5)) * ((sin(theta * 0.5))^6.0),
                        ),
                    ),
                ),
            ),
        )
    end

    Ylm = ifelse(m == 0, Ylm, Ylm * exp(1im * m * pi/2))
    return Ylm
end

#*************  ANSATZ INTEGRATED PHASE *************

function RD_Phase_AnsatzInt(model::PhenomXHM, ff, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)

    invf   = ff^(-1);
    frd=fring
    fda=fdamp

    phaseRD = 0.0
    if ModeMixingOn == false
        #   rescaling of the 22 ansatz -- used for (21),(33),(44) 
        # ansatz:
        #  alpha0 f - fRDlm^2*alpha2)/f  + alphaL*ArcTan[(f - fRDlm)/fdamplm]
        phaseRD = alpha0*ff -frd*frd *(alpha2)*invf +  (alphaL)* atan((ff-frd)/fda);
    else
        #  calibration of spheroidal ringdown waveform for (32) 
        invf3 = ff^(-3);

        # if IMRPhenomXHMRingdownPhaseVersion == 122019
        #  ansatz: f alpha0 - (alpha4)/(3 f^3) - (alpha2)/f + alphaL ArcTan[(f - fRDlm)/fdamplm]
        alpha0_S = ModeMixingCoeffs[1]
        alphaL_S = ModeMixingCoeffs[2]
        alpha2_S = ModeMixingCoeffs[3]
        alpha4_S = ModeMixingCoeffs[4]
        phi0_S = ModeMixingCoeffs[5]
        phaseRD = phi0_S+alpha0_S*ff -(alpha2_S)*invf -1. ./3. .*(alpha4_S)*invf3 +(alphaL_S)* atan((ff-frd)/fda);
        
        
    end

    return phaseRD
end



#*************  ANSATZ PHASE DERIVATIVE *************

function RD_Phase_Ansatz(model::PhenomXHM, ff, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)

    frd   = fring
    fda   = fdamp

    if ModeMixingOn == false

        #  rescaling of the 22 ansatz -- used for (21),(33),(44) 
        #ansatz:
        #alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)
        dphaseRD = ( alpha0 +  frd*frd*(alpha2)*ff^(-2) + ( (alphaL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );

    else
        #  calibration of spheroidal ringdown waveform for (32) 

        #if IMRPhenomXHMRingdownPhaseVersion == 122019
        # ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)
        alpha0_S = ModeMixingCoeffs[1]
        alphaL_S = ModeMixingCoeffs[2]
        alpha2_S = ModeMixingCoeffs[3]
        alpha4_S = ModeMixingCoeffs[4]
        phi0_S = ModeMixingCoeffs[5]
        dphaseRD = ( alpha0_S +  (alpha2_S)*ff^(-2) + (alpha4_S)*ff^(-4) +( (alphaL_S)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
        
    end
    return dphaseRD;
end

function hphc(model::PhenomXHM, f, mc, eta, chi1, chi2, dL, iota;
                debug=false, return_modes=false, final_spin_override=nothing,
                fRef=nothing,
                fcutPar = 0.3, GMsun_over_c3 = uc.GMsun_over_c3, GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc)


    M = mc / (eta^(0.6))
    eta1 = eta
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    eta3 = eta * eta2
    eta4 = eta * eta3
    eta5 = eta * eta4
    eta6 = eta * eta5
    eta7 = eta * eta6

    etaInv = 1 ./ eta
    sqroot = sqrt(eta)

    typeofFD = typeof(eta) # This is the type of this variable when ForwardDiff is used, 
    # this is used to define the type of all the variables inside this function, keeping track of the Tag of the variables

    if debug == true
        println("typeofFD = $typeofFD")
    end

    if eta < 0.05
        println("Extreme mass ratio limit")
        println("Model not implemented yet")
    end

    #distance = dL * uGpc        #distance is dL in meters
    Mtot = M 

    amp0    =  M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL
    ampNorm = sqrt(2. *eta/3.) * (pi^(-1. /6.))
    factor_22 = 2. * sqrt(5. /(64. *pi))

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)  

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)

    f_min = minimum(f)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    Mf = fgrid
    len = length(fgrid)
    offset = 0 # offset for the frequency series from SetupWFArrays 

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    m1ByMSq = m1ByM*m1ByM
    m2ByMSq = m2ByM*m2ByM

    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2 = chi1*chi2

    totchi = ((m1ByMSq * chi1 + m2ByMSq * chi2) / (m1ByMSq + m2ByMSq))
    STotR = totchi
    
    totchi2 = totchi*totchi
    dchi = chi1-chi2
    # Normalised PN reduced spin parameter
    chiS = 0.5 * (chi1 + chi2)
    chiA = 0.5 * (chi1 - chi2)
    chidiff1 = chiA
    chidiff2 = chiA^2
    chi_eff = (m1ByM*chi1 + m2ByM*chi2)
    chiPN   = (chi_eff - (38. /113.)*eta*(chi1 + chi2)) / (1. - (76. *eta/113.))
    chiPN1  = chiPN
    chiPN2  = chiPN*chiPN
    chiPN3  = chiPN*chiPN2

    S = STotR;
    S1 = STotR
    S2 = STotR*STotR
    S3 = STotR*S2
    S4 = STotR*S3
    S5 = STotR*S4

    q = 0.5*(1.0 + Seta - 2.0*eta)/eta

    # PN symmetry coefficient
    delta = Seta

    phiNorm = - (3. * (pi^(-5. /3.)))/ 128.
    dphase0     = 5.0 / (128.0 * pi^(5. ./3.))
    Ampzero = 0

    fRef      = isnothing(fRef) ? f_min : fRef
    MfRef = fRef * M * GMsun_over_c3 
    phi0 = 0.0

    etaEMR=0.05 # Extreme mass ratio limit
    if eta < etaEMR
        println("Extreme mass ratio limit")
        println("Model not implemented yet")
    end


    fMECO    = (((0.018744340279608845 + 0.0077903147004616865*eta + 0.003940354686136861*eta2 - 0.00006693930988501673*eta2*eta)/(1. - 0.10423384680638834*eta)) + ((chiPN*(0.00027180386951683135 - 0.00002585252361022052*chiPN + eta2*eta2*(-0.0006807631931297156 + 0.022386313074011715*chiPN - 0.0230825153005985*chiPN2) + eta2*(0.00036556167661117023 - 0.000010021140796150737*chiPN - 0.00038216081981505285*chiPN2) + eta*(0.00024422562796266645 - 0.00001049013062611254*chiPN - 0.00035182990586857726*chiPN2) + eta2*eta*(-0.0005418851224505745 + 0.000030679548774047616*chiPN + 4.038390455349854e-6*chiPN2) - 0.00007547517256664526*chiPN2))/(0.026666543809890402 + (-0.014590539285641243 - 0.012429476486138982*eta + 1.4861197211952053*eta2*eta2 + 0.025066696514373803*eta2 + 0.005146809717492324*eta2*eta)*chiPN + (-0.0058684526275074025 - 0.02876774751921441*eta - 2.551566872093786*eta2*eta2 - 0.019641378027236502*eta2 - 0.001956646166089053*eta2*eta)*chiPN2 + (0.003507640638496499 + 0.014176504653145768*eta + 1. *eta2*eta2 + 0.012622225233586283*eta2 - 0.00767768214056772*eta2*eta)*chiPN2*chiPN)) + (dchi*dchi*(0.00034375176678815234 + 0.000016343732281057392*eta)*eta2 + dchi*Seta*eta*(0.08064665214195679*eta2 + eta*(-0.028476219509487793 - 0.005746537021035632*chiPN) - 0.0011713735642446144*chiPN)))
    afinal   = (((3.4641016151377544*eta + 20.0830030082033*eta2 - 12.333573402277912*eta2*eta)/(1 + 7.2388440419467335*eta)) + ((m1ByMSq + m2ByMSq)*totchi + ((-0.8561951310209386*eta - 0.09939065676370885*eta2 + 1.668810429851045*eta2*eta)*totchi + (0.5881660363307388*eta - 2.149269067519131*eta2 + 3.4768263932898678*eta2*eta)*totchi2 + (0.142443244743048*eta - 0.9598353840147513*eta2 + 1.9595643107593743*eta2*eta)*totchi2*totchi) / (1 + (-0.9142232693081653 + 2.3191363426522633*eta - 9.710576749140989*eta2*eta)*totchi)) + (0.3223660562764661*dchi*Seta*(1 + 9.332575956437443*eta)*eta2 - 0.059808322561702126*dchi*dchi*eta2*eta + 2.3170397514509933*dchi*Seta*(1 - 3.2624649875884852*eta)*eta2*eta*totchi))
    afinal = isnothing(final_spin_override) ? afinal : final_spin_override
    Erad   = ((((0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta2*eta + 3.145145224278187*eta2*eta2)*(1. + (-0.13084389181783257 - 1.1387311580238488*eta + 5.49074464410971*eta2)*totchi + (-0.17762802148331427 + 2.176667900182948*eta2)*totchi2 + (-0.6320191645391563 + 4.952698546796005*eta - 10.023747993978121*eta2)*totchi*totchi2)) / (1. + (-0.9919475346968611 + 0.367620218664352*eta + 4.274567337924067*eta2)*totchi)) + (- 0.09803730445895877*dchi*Seta*(1. - 3.2283713377939134*eta)*eta2 + 0.01118530335431078*dchi*dchi*eta2*eta - 0.01978238971523653*dchi*Seta*(1. - 4.91667749015812*eta)*eta*totchi))
    Mfinal = 1. - Erad

    # Fitting function for hybrid minimum energy circular orbit (MECO) function and computation of ISCO frequency
    Z1tmp = 1. + cbrt((1. - afinal*afinal) ) * (cbrt(1. + afinal) + cbrt(1. - afinal))
    Z1tmp = ifelse(Z1tmp>3., 3., Z1tmp)
    Z2tmp = sqrt(3. *afinal*afinal + Z1tmp*Z1tmp)
    fISCO  = (1. / ((3. + Z2tmp - sign(afinal)*sqrt((3. - Z1tmp) * (3. + Z1tmp + 2. *Z2tmp)))^(3. /2.) + afinal))/pi


    ### fring and fdamp for each mode

    x2= afinal*afinal;
    x3= x2*afinal;
    x4= x2*x2;
    x5= x3*x2;
    x6= x3*x3;

    fring_21 = (0.059471695665734674 - 0.07585416297991414*afinal + 0.021967909664591865*x2 - 0.0018964744613388146*x3 + 0.001164879406179587*x4 - 0.0003387374454044957*x5)/(1 - 1.4437415542456158*afinal + 0.49246920313191234*x2)/Mfinal
    fring_22 = ((0.05947169566573468 - 0.14989771215394762*afinal + 0.09535606290986028*afinal*afinal + 0.02260924869042963*afinal*afinal*afinal - 0.02501704155363241*afinal*afinal*afinal*afinal - 0.005852438240997211*(afinal^5) + 0.0027489038393367993*(afinal^6) + 0.0005821983163192694*(afinal^7))/(1 - 2.8570126619966296*afinal + 2.373335413978394*afinal*afinal - 0.6036964688511505*afinal*afinal*afinal*afinal + 0.0873798215084077*(afinal^6)))/(1. - Erad)
    fring_33 = (0.09540436245212061 - 0.22799517865876945*afinal + 0.13402916709362475*x2 + 0.03343753057911253*x3 - 0.030848060170259615*x4 - 0.006756504382964637*x5 + 0.0027301732074159835*x6)/(1 - 2.7265947806178334*afinal + 2.144070539525238*x2 - 0.4706873667569393*x4 + 0.05321818246993958*x6)/Mfinal
    fring_32 = (0.09540436245212061 - 0.13628306966373951*afinal + 0.030099881830507727*x2 - 0.000673589757007597*x3 + 0.0118277880067919*x4 + 0.0020533816327907334*x5 - 0.0015206141948469621*x6)/(1 - 1.6531854335715193*afinal + 0.5634705514193629*x2 + 0.12256204148002939*x4 - 0.027297817699401976*x6)/Mfinal
    fring_44 = (0.1287821193485683 - 0.21224284094693793*afinal + 0.0710926778043916*x2 + 0.015487322972031054*x3 - 0.002795401084713644*x4 + 0.000045483523029172406*x5 + 0.00034775290179000503*x6)/(1 - 1.9931645124693607*afinal + 1.0593147376898773*x2 - 0.06378640753152783*x4)/Mfinal
    fring_array = [fring_21, fring_33, fring_32, fring_44]

    fdamp_21 = (2.0696914454467294 - 3.1358071947583093*afinal + 0.14456081596393977*x2 + 1.2194717985037946*x3 - 0.2947372598589144*x4 + 0.002943057145913646*x5)/(146.1779212636481 - 219.81790388304876*afinal + 17.7141194900164*x2 + 75.90115083917898*x3 - 18.975287709794745*x4)/Mfinal
    fdamp_22 = ((0.014158792290965177 - 0.036989395871554566*afinal + 0.026822526296575368*afinal*afinal + 0.0008490933750566702*afinal*afinal*afinal - 0.004843996907020524*afinal*afinal*afinal*afinal - 0.00014745235759327472*(afinal^5) + 0.0001504546201236794*(afinal^6))/(1 - 2.5900842798681376*afinal + 1.8952576220623967*afinal*afinal - 0.31416610693042507*afinal*afinal*afinal*afinal + 0.009002719412204133*(afinal^6)))/(1. - Erad)
    fdamp_33 = (0.014754148319335946 - 0.03124423610028678*afinal + 0.017192623913708124*x2 + 0.001034954865629645*x3 - 0.0015925124814622795*x4 - 0.0001414350555699256*x5)/(1 - 2.0963684630756894*afinal + 1.196809702382645*x2 - 0.09874113387889819*x4)/Mfinal
    fdamp_32 = (0.014754148319335946 - 0.03445752346074498*afinal + 0.02168855041940869*x2 + 0.0014945908223317514*x3 - 0.0034761714223258693*x4)/(1 - 2.320722660848874*afinal + 1.5096146036915865*x2 - 0.18791187563554512*x4)/Mfinal
    fdamp_44 = (0.014986847152355699 - 0.01722587715950451*afinal - 0.0016734788189065538*x2 + 0.0002837322846047305*x3 + 0.002510528746148588*x4 + 0.00031983835498725354*x5 + 0.000812185411753066*x6)/(1 - 1.1350205970682399*afinal - 0.0500827971270845*x2 + 0.13983808071522857*x4 + 0.051876225199833995*x6)/Mfinal
    fdamp_array = [fdamp_21, fdamp_33, fdamp_32, fdamp_44]



    ### 22 mode

    
    Phase22 = Phase_22_ConnectionCoefficients(mc, eta, chi1, chi2;
                                              final_spin_override=final_spin_override)
    Amp22 = Ampl_22_ConnectionCoefficients(mc, eta, chi1, chi2, dL;
                                            final_spin_override=final_spin_override)

    psi4tostrain = ((13.39320482758057 - 175.42481512989315*eta + 2097.425116152503*eta2 - 9862.84178637907*eta2*eta + 16026.897939722587*eta2*eta2) + ((4.7895602776763 - 163.04871764530466*eta + 609.5575850476959*eta2)*totchi + (1.3934428041390161 - 97.51812681228478*eta + 376.9200932531847*eta2)*totchi2 + (15.649521097877374 + 137.33317057388916*eta - 755.9566456906406*eta2)*totchi2*totchi + (13.097315867845788 + 149.30405703643288*eta - 764.5242164872267*eta2)*totchi2*totchi2) + (105.37711654943146*dchi*Seta*eta2))
    DeltaT =  -2. *pi*(500+psi4tostrain);

    lina = 0.0

    timeshift= TimeShift_22(model, eta, Seta, S, dchi, fring_22, fdamp_22, Phase22)
    phifRef = -(etaInv*_completePhase(model, MfRef, Phase22, fdamp_22, fring_22) + timeshift*MfRef + lina) + pi/4. + pi
    
    
    phase_22    =  etaInv .*_completePhase(model, fgrid, Phase22, fdamp_22, fring_22) .+ ifelse.(fgrid .<= fcutPar, timeshift .*fgrid .+ lina .+ phifRef, 0.)
    ampl_22 = Ampl(PhenomXAS(), f, mc, eta, chi1, chi2, dL;
                   fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3,
                   GMsun_over_c2_Gpc=GMsun_over_c2_Gpc,
                   final_spin_override=final_spin_override)

    if debug == true
        for j in 1:len
            println("data$j = [$(f[j]), $(ampl_22[j]), $(phase_22[j])]")
        end
    end

    Ym = _spinWeightedSphericalHarmonic(model, iota, 22, -2)    # since _spinWeightedSphericalHarmonic() is first called with -m
        
    # Equatorial symmetry: add in -m and m mode 
    Ystar = conj(_spinWeightedSphericalHarmonic(model, iota, 22, 2))
    factorp = 0.5 * (Ym +  Ystar);
    factorc = 1im * 0.5 * ( Ym -  Ystar);



    htildelm =  @. -1. *ampl_22 * exp(1im * phase_22) /factor_22 
    wf22 = - htildelm
    coprecessing_modes = return_modes ? Dict{Tuple{Int,Int},typeof(htildelm)}((2, -2) => copy(htildelm)) : nothing

    ### First calculate the 22 mode
    hp = Vector{Complex{typeofFD}}(undef, len) #initialize the vector
    hc = Vector{Complex{typeofFD}}(undef, len)
    if debug == true
        println("htildelm1 =", htildelm[1])
        println("hp1 =", hp[1])
    end

    for j in 1:len
        hlm = htildelm[j];
        hp[j] = (factorp * hlm);
        hc[j] = (factorc * hlm);
    end




    if debug == true
        println()
        println("ell_emm = 22")
        println("Ym = $Ym")
        println("Ystar = $Ystar")
        println("factorp = $factorp")
        println("factorc = $factorc\n")
    end




    psi4tostrain = ((13.39320482758057 - 175.42481512989315*eta + 2097.425116152503*eta2 - 9862.84178637907*eta2*eta + 16026.897939722587*eta2*eta2) + ((4.7895602776763 - 163.04871764530466*eta + 609.5575850476959*eta2)*totchi + (1.3934428041390161 - 97.51812681228478*eta + 376.9200932531847*eta2)*totchi2 + (15.649521097877374 + 137.33317057388916*eta - 755.9566456906406*eta2)*totchi2*totchi + (13.097315867845788 + 149.30405703643288*eta - 764.5242164872267*eta2)*totchi2*totchi2) + (105.37711654943146*dchi*Seta*eta2))
    DeltaT =  -2. *pi*(500+psi4tostrain);

    ### coefficients for mode mixing, see IMRPhenomXHM_Initialize_MixingCoeffs
    re_l2m2lp2 = (1 - 2.2956993576253635*afinal + 1.461988775298876*x2 + 0.0043296365593147035*x3 - 0.1695667458204109*x4 - 0.0006267849034466508*x5)/(1 - 2.2956977727459043*afinal + 1.4646339137818438*x2 - 0.16843226886562457*x4 - 0.00007150540890128118*x6);
    im_l2m2lp2 = (afinal*(0.3826673013161342 - 0.47531267226013896*afinal - 0.05898102880105067*x2 + 0.0724525431346487*x3 + 0.054714637311702986*x4 + 0.024544862718252784*x5))/(-38.70835035062785 + 69.82140084545878*afinal - 27.99036444363243*x2 - 4.152310472191899*x4 + 1. *x6);
    re_l3m2lp2 = (afinal*(0.47513455283841244 - 0.9016636384605536*afinal + 0.3844811236426182*x2 + 0.0855565148647794*x3 - 0.03620067426672167*x4 - 0.006557249133752502*x5))/(-6.76894063440646 + 15.170831931186493*afinal - 9.406169787571082*x2 + 1. *x4);
    im_l3m2lp2 = (afinal*(-2.8704762147145533 + 4.436434016918535*afinal - 1.0115343326360486*x2 - 0.08965314412106505*x3 - 0.4236810894599512*x4 - 0.041787576033810676*x5))/(-171.80908957903395 + 272.362882450877*afinal - 76.68544453077854*x2 - 25.14197656531123*x4 + 1. *x6);
    re_l2m2lp3 = (afinal*(18.522563276099167 - 37.978140351289014*afinal + 19.030390708998894*x2 + 3.0355668591803386*x3 - 2.210028290847915*x4 - 0.37117112862247975*x5))/(164.52480238697507 - 377.9093045285145*afinal + 243.3353695550844*x2 - 30.79738566181734*x4 + 1. *x6);
    im_l2m2lp3 = (afinal*(-49.7688437256778 + 120.43773704442333*afinal - 82.95323455645332*x2 + 1.721453011852496*x3 + 11.540237244397877*x4 - 0.9819458637589314*x5))/(2858.5790831181725 - 6305.619505422591*afinal + 3825.6742092829054*x2 - 377.7822297815406*x4 + 1. *x6);
    re_l3m2lp3 = (1 - 2.107852425643677*afinal + 1.1906393634562715*x2 + 0.02244848864087732*x3 - 0.09593447799423722*x4 - 0.0021343381708933025*x5 - 0.005319515989331159*x6)/(1 - 2.1078515887706324*afinal + 1.2043484690080966*x2 - 0.08910191596778137*x4 - 0.005471749827809503*x6);
    im_l3m2lp3 = (afinal*(12.45701482868677 - 29.398484595717147*afinal + 18.26221675782779*x2 + 1.9308599142669403*x3 - 3.159763242921214*x4 - 0.0910871567367674*x5))/(345.52914639836257 - 815.4349339779621*afinal + 538.3888932415709*x2 - 69.3840921447381*x4 + 1. *x6);
    


    mixingCoeffs1 = re_l2m2lp2 + 1im * im_l2m2lp2 # not used
    mixingCoeffs2 = re_l2m2lp3 + 1im * im_l2m2lp3 # not used
    mixingCoeffs3 = re_l3m2lp2 + 1im * im_l3m2lp2
    mixingCoeffs4 = re_l3m2lp3 + 1im * im_l3m2lp3

    mixingCoeffs = [mixingCoeffs1, mixingCoeffs2, mixingCoeffs3, mixingCoeffs4]
  
    # Adjust conventions so that they match the ones used for the hybrids
    mixingCoeffs[3]= -1. * mixingCoeffs[3];
    mixingCoeffs[4]= -1. * mixingCoeffs[4];

    if debug == true
        println("mixingCoeffs: ", mixingCoeffs)
        println("mixingCoeffs[3]: ", mixingCoeffs[3])
        println("mixingCoeffs[4]: ", mixingCoeffs[4])
    end

    IMRPhenomXHMInspiralPhaseVersion = 122019
    IMRPhenomXHMIntermediatePhaseVersion = 122019
    IMRPhenomXHMRingdownPhaseVersion = 122019
    IMRPhenomXHMReleaseVersion = 122022

    IMRPhenomXHMInspiralAmpFitsVersion     = 122022;
    IMRPhenomXHMIntermediateAmpFitsVersion = 122022;
    IMRPhenomXHMRingdownAmpFitsVersion     = 122022;


    IMRPhenomXHMInspiralPhaseFitsVersion = IMRPhenomXHMInspiralPhaseVersion;
    IMRPhenomXHMIntermediatePhaseFitsVersion = IMRPhenomXHMIntermediatePhaseVersion;
    IMRPhenomXHMRingdownPhaseFitsVersion = IMRPhenomXHMRingdownPhaseVersion;
  
    IMRPhenomXHMInspiralPhaseFreqsVersion = IMRPhenomXHMInspiralPhaseVersion;
    IMRPhenomXHMIntermediatePhaseFreqsVersion = IMRPhenomXHMIntermediatePhaseVersion;
    IMRPhenomXHMRingdownPhaseFreqsVersion = IMRPhenomXHMRingdownPhaseVersion;



    #= Take all the phenom coefficients accross the three regions (inspiral, intermeidate and ringdown) and all the needed parameters to reconstruct the amplitude (including mode-mixing). =#
    index = 0

    for ell_emm in [21,33,32,44]
        if debug == true
            println("ell_emm: ", ell_emm) 
        end


        htildelm__ = Vector{Complex{typeofFD}}(undef, len)
        
        ModeMixingCoeffs = zeros(5)
        alpha0, alpha2, alphaL = 0., 0., 0.
        index += 1

        ell = ell_emm .÷ 10
        emm = ell_emm .% 10

        # Multiply by (-1)^l to get the true h_l-m(f) 
        if ell%2 != 0
            minus1l = -1;
        else
            minus1l = 1;
        end


        Amp0 = minus1l * amp0; #//Transform NR units to Physical units

        # watch out for Amp0 or amp0, they differ by minus1l
    
        #  Loop over only positive m is intentional. The single mode function returns the negative mode h_l-m, and the positive is added automatically in IMRPhenomXHMFDAddMode. 
        #  First check if (l,m) mode is 'activated' in the ModeArray 
        #  if activated then generate the mode, else skip this mode. 
    
        fring = fring_array[index]
        fdamp = fdamp_array[index]

        fMECOlm = fMECO*emm*0.5


        if ell_emm == 21
            ModeMixingOn=false;
            if q == 1. && chi1 == chi2  # Odd mode for equal mass, equal spin is zero
                Ampzero = 1;
            end
        elseif ell_emm == 33
            ModeMixingOn=false;
            if q == 1. && chi1 == chi2 # Odd mode for equal mass, equal spin is zero
                Ampzero = 1;
            end
        elseif ell_emm == 32    # Mode with Mixing
            ModeMixingOn=true;
        
        elseif ell_emm == 44
            ModeMixingOn=false;
        end

        if Ampzero == 1 # Skip the mode
            continue
        end


        ###  from IMRPhenomXHM_FillAmpFitsArray

        ### Collocation points INSPIRAL AMPLITUDE
        if ell_emm == 21
            Insp_v1 = abs(chidiff1*eta5*(-3962.5020052272976 + 987.635855365408*chiPN1 - 134.98527058315528*chiPN2) + delta*(19.30531354642419 + 16.6640319856064*eta1 - 120.58166037019478*eta2 + 220.77233521626252*eta3)*sqroot + chidiff1*delta*(31.364509907424765*eta1 - 843.6414532232126*eta2 + 2638.3077554662905*eta3)*sqroot + chidiff1*delta*(32.374226994179054*eta1 - 202.86279451816662*eta2 + 347.1621871204769*eta3)*chiPN1*sqroot + delta*chiPN1*(-16.75726972301224*(1.1787350890261943 - 7.812073811917883*eta1 + 99.47071002831267*eta2 - 500.4821414428368*eta3 + 876.4704270866478*eta4) + 2.3439955698372663*(0.9373952326655807 + 7.176140122833879*eta1 - 279.6409723479635*eta2 + 2178.375177755584*eta3 - 4768.212511142035*eta4)*chiPN1)*sqroot)

            Insp_v2 = abs(chidiff1*eta5*(-2898.9172078672705 + 580.9465034962822*chiPN1 + 22.251142639924076*chiPN2) + delta*(chidiff2*(-18.541685007214625*eta1 + 166.7427445020744*eta2 - 417.5186332459383*eta3) + chidiff1*(41.61457952037761*eta1 - 779.9151607638761*eta2 + 2308.6520892707795*eta3))*sqroot + delta*(11.414934585404561 + 30.883118528233638*eta1 - 260.9979123967537*eta2 + 1046.3187137392433*eta3 - 1556.9475493549746*eta4)*sqroot + delta*chiPN1*(-10.809007068469844*(1.1408749895922659 - 18.140470190766937*eta1 + 368.25127088896744*eta2 - 3064.7291458207815*eta3 + 11501.848278358668*eta4 - 16075.676528787526*eta5) + 1.0088254664333147*(1.2322739396680107 - 192.2461213084741*eta1 + 4257.760834055382*eta2 - 35561.24587952242*eta3 + 130764.22485304279*eta4 - 177907.92440833704*eta5)*chiPN1)*sqroot + delta*(chidiff1*(36.88578491943111*eta1 - 321.2569602623214*eta2 + 748.6659668096737*eta3)*chiPN1 + chidiff1*(-95.42418611585117*eta1 + 1217.338674959742*eta2 - 3656.192371615541*eta3)*chiPN2)*sqroot)

            Insp_v3 = abs(chidiff1*eta5*(-2282.9983216879655 + 157.94791186394787*chiPN1 + 16.379731479465033*chiPN2) + chidiff1*delta*(21.935833431534224*eta1 - 460.7130131927895*eta2 + 1350.476411541137*eta3)*sqroot + delta*(5.390240326328237 + 69.01761987509603*eta1 - 568.0027716789259*eta2 + 2435.4098320959706*eta3 - 3914.3390484239667*eta4)*sqroot + chidiff1*delta*(29.731007410186827*eta1 - 372.09609843131386*eta2 + 1034.4897198648962*eta3)*chiPN1*sqroot + delta*chiPN1*(-7.1976397556450715*(0.7603360145475428 - 6.587249958654174*eta1 + 120.87934060776237*eta2 - 635.1835857158857*eta3 + 1109.0598539312573*eta4) - 0.0811847192323969*(7.951454648295709 + 517.4039644814231*eta1 - 9548.970156895082*eta2 + 52586.63520999897*eta3 - 93272.17990295641*eta4)*chiPN1 - 0.28384547935698246*(-0.8870770459576875 + 180.0378964169756*eta1 - 2707.9572896559484*eta2 + 14158.178124971111*eta3 - 24507.800226675925*eta4)*chiPN2)*sqroot)

        elseif ell_emm == 33
            Insp_v1 = chidiff1*eta5*(155.1434307076563 + 26.852777193715088*chiPN1 + 1.4157230717300835*chiPN2) + chidiff1*delta*(6.296698171560171*eta1 + 15.81328761563562*eta2 - 141.85538063933927*eta3)*sqroot + delta*(20.94372147101354 + 68.14577638017842*eta1 - 898.470298591732*eta2 + 4598.64854748635*eta3 - 8113.199260593833*eta4)*sqroot + chidiff1*delta*(29.221863857271703*eta1 - 348.1658322276406*eta2 + 965.4670353331536*eta3)*chiPN1*sqroot + delta*chiPN1*(-9.753610761811967*(1.7819678168496158 - 44.07982999150369*eta1 + 750.8933447725581*eta2 - 5652.44754829634*eta3 + 19794.855873435758*eta4 - 26407.40988450443*eta5) + 0.014210376114848208*(-196.97328616330392 + 7264.159472864562*eta1 - 125763.47850622259*eta2 + 1.1458022059130718e6*eta3 - 4.948175330328345e6*eta4 + 7.911048294733888e6*eta5)*chiPN1 - 0.26859293613553986*(-8.029069605349488 + 888.7768796633982*eta1 - 16664.276483466252*eta2 + 128973.72291098491*eta3 - 462437.2690007375*eta4 + 639989.1197424605*eta5)*chiPN2)*sqroot

            Insp_v2 = chidiff1*eta5*(161.62678370819597 + 37.141092711336846*chiPN1 - 0.16889712161410445*chiPN2) + chidiff1*delta*(3.4895829486899825*eta1 + 51.07954458810889*eta2 - 249.71072528701757*eta3)*sqroot + delta*(12.501397517602173 + 35.75290806646574*eta1 - 357.6437296928763*eta2 + 1773.8883882162215*eta3 - 3100.2396041211605*eta4)*sqroot + chidiff1*delta*(13.854211287141906*eta1 - 135.54916401086845*eta2 + 327.2467193417936*eta3)*chiPN1*sqroot + delta*chiPN1*(-5.2580116732827085*(1.7794900975289085 - 48.20753331991333*eta1 + 861.1650630146937*eta2 - 6879.681319382729*eta3 + 25678.53964955809*eta4 - 36383.824902258915*eta5) + 0.028627002336747746*(-50.57295946557892 + 734.7581857539398*eta1 - 2287.0465658878725*eta2 + 15062.821881048358*eta3 - 168311.2370167227*eta4 + 454655.37836367317*eta5)*chiPN1 - 0.15528289788512326*(-12.738184090548508 + 1129.44485109116*eta1 - 25091.14888164863*eta2 + 231384.03447562453*eta3 - 953010.5908118751*eta4 + 1.4516597366230418e6*eta5)*chiPN2)*sqroot

            Insp_v3 = chidiff1*delta*(-0.5869777957488564*eta1 + 32.65536124256588*eta2 - 110.10276573567405*eta3) + chidiff1*delta*(3.524800489907584*eta1 - 40.26479860265549*eta2 + 113.77466499598913*eta3)*chiPN1 + delta*chiPN1*(-1.2846335585108297*(0.09991079016763821 + 1.37856806162599*eta1 + 23.26434219690476*eta2 - 34.842921754693386*eta3 - 70.83896459998664*eta4) - 0.03496714763391888*(-0.230558571912664 + 188.38585449575902*eta1 - 3736.1574640444287*eta2 + 22714.70643022915*eta3 - 43221.0453556626*eta4)*chiPN1) + chidiff1*eta7*(2667.3441342894776 + 47.94869769580204*chidiff2 + 793.5988192446642*chiPN1 + 293.89657731755483*chiPN2) + delta*(5.148353856800232 + 148.98231189649468*eta1 - 2774.5868652930294*eta2 + 29052.156454239772*eta3 - 162498.31493332976*eta4 + 460912.76402476896*eta5 - 521279.50781871413*eta6)*sqroot

        elseif ell_emm == 32    
            Insp_v1 = (chidiff1*delta*(-0.739317114582042*eta1 - 47.473246070362634*eta2 + 278.9717709112207*eta3 - 566.6420939162068*eta4) + chidiff2*(-0.5873680378268906*eta1 + 6.692187014925888*eta2 - 24.37776782232888*eta3 + 23.783684827838247*eta4))*sqroot + (3.2940434453819694 + 4.94285331708559*eta1 - 343.3143244815765*eta2 + 3585.9269057886418*eta3 - 19279.186145681153*eta4 + 51904.91007211022*eta5 - 55436.68857586653*eta6)*sqroot + chidiff1*delta*(12.488240781993923*eta1 - 209.32038774208385*eta2 + 1160.9833883184604*eta3 - 2069.5349737049073*eta4)*chiPN1*sqroot + chiPN1*(0.6343034651912586*(-2.5844888818001737 + 78.98200041834092*eta1 - 1087.6241783616488*eta2 + 7616.234910399297*eta3 - 24776.529123239357*eta4 + 30602.210950069973*eta5) - 0.062088720220899465*(6.5586380356588565 + 36.01386705325694*eta1 - 3124.4712274775407*eta2 + 33822.437731298516*eta3 - 138572.93700180828*eta4 + 198366.10615196894*eta5)*chiPN1)*sqroot

            Insp_v2 = (chidiff2*(-0.03940151060321499*eta1 + 1.9034209537174116*eta2 - 8.78587250202154*eta3) + chidiff1*delta*(-1.704299788495861*eta1 - 4.923510922214181*eta2 + 0.36790005839460627*eta3))*sqroot + (2.2911849711339123 - 5.1846950040514335*eta1 + 60.10368251688146*eta2 - 1139.110227749627*eta3 + 7970.929280907627*eta4 - 25472.73682092519*eta5 + 30950.67053883646*eta6)*sqroot + chiPN1*(0.7718201508695763*(-1.3012906461000349 + 26.432880113146012*eta1 - 186.5001124789369*eta2 + 712.9101229418721*eta3 - 970.2126139442341*eta4) + 0.04832734931068797*(-5.9999628512498315 + 78.98681284391004*eta1 + 1.8360177574514709*eta2 - 2537.636347529708*eta3 + 6858.003573909322*eta4)*chiPN1)*sqroot

            Insp_v3 = (chidiff2*(-0.6358511175987503*eta1 + 5.555088747533164*eta2 - 14.078156877577733*eta3) + chidiff1*delta*(0.23205448591711159*eta1 - 19.46049432345157*eta2 + 36.20685853857613*eta3))*sqroot + (1.1525594672495008 + 7.380126197972549*eta1 - 17.51265776660515*eta2 - 976.9940395257111*eta3 + 8880.536804741967*eta4 - 30849.228936891763*eta5 + 38785.53683146884*eta6)*sqroot + chidiff1*delta*(1.904350804857431*eta1 - 25.565242391371093*eta2 + 80.67120303906654*eta3)*chiPN1*sqroot + chiPN1*(0.785171689871352*(-0.4634745514643032 + 18.70856733065619*eta1 - 167.9231114864569*eta2 + 744.7699462372949*eta3 - 1115.008825153004*eta4) + 0.13469300326662165*(-2.7311391326835133 + 72.17373498208947*eta1 - 483.7040402103785*eta2 + 1136.8367114738041*eta3 - 472.02962341590774*eta4)*chiPN1)*sqroot

        elseif ell_emm == 44
            Insp_v1 = (chidiff1*delta*(0.5697308729057493*eta1 + 8.895576813118867*eta2 - 34.98399465240273*eta3) + chidiff2*(1.6370346538130884*eta1 - 14.597095790380884*eta2 + 33.182723737396294*eta3))*sqroot + (5.2601381002242595 - 3.557926105832778*eta1 - 138.9749850448088*eta2 + 603.7453704122706*eta3 - 923.5495700703648*eta4)*sqroot + chiPN1*(-0.41839636169678796*(5.143510231379954 + 104.62892421207803*eta1 - 4232.508174045782*eta2 + 50694.024801783446*eta3 - 283097.33358214336*eta4 + 758333.2655404843*eta5 - 788783.0559069642*eta6) - 0.05653522061311774*(5.605483124564013 + 694.00652410087*eta1 - 17551.398321516353*eta2 + 165236.6480734229*eta3 - 761661.9645651339*eta4 + 1.7440315410044065e6*eta5 - 1.6010489769238676e6*eta6)*chiPN1 - 0.023693246676754775*(16.437107575918503 - 2911.2154288136217*eta1 + 89338.32554683842*eta2 - 1.0803340811860575e6*eta3 + 6.255666490084672e6*eta4 - 1.7434160932177313e7*eta5 + 1.883460394974573e7*eta6)*chiPN2)*sqroot

            Insp_v2 = (chidiff2*(-0.8318312659717388*eta1 + 7.6541168007977864*eta2 - 16.648660653220123*eta3) + chidiff1*delta*(2.214478316304753*eta1 - 7.028104574328955*eta2 + 5.56587823143958*eta3))*sqroot + (3.173191054680422 + 6.707695566702527*eta1 - 155.22519772642607*eta2 + 604.0067075996933*eta3 - 876.5048298377644*eta4)*sqroot + chidiff1*delta*(4.749663394334708*eta1 - 42.62996105525792*eta2 + 97.01712147349483*eta3)*chiPN1*sqroot + chiPN1*(-0.2627203100303006*(6.460396349297595 - 52.82425783851536*eta1 - 552.1725902144143*eta2 + 12546.255587592654*eta3 - 81525.50289542897*eta4 + 227254.37897941095*eta5 - 234487.3875219032*eta6) - 0.008424003742397579*(-109.26773035716548 + 15514.571912666677*eta1 - 408022.6805482195*eta2 + 4.620165968920881e6*eta3 - 2.6446950627957724e7*eta4 + 7.539643948937692e7*eta5 - 8.510662871580401e7*eta6)*chiPN1 - 0.008830881730801855*(-37.49992494976597 + 1359.7883958101172*eta1 - 23328.560285901796*eta2 + 260027.4121353132*eta3 - 1.723865744472182e6*eta4 + 5.858455766230802e6*eta5 - 7.756341721552802e6*eta6)*chiPN2 - 0.027167813927224657*(34.281932237450256 - 3312.7658728016568*eta1 + 84126.14531363266*eta2 - 956052.0170024392*eta3 + 5.570748509263883e6*eta4 - 1.6270212243584689e7*eta5 + 1.8855858173287075e7*eta6)*chiPN3)*sqroot

            Insp_v3 = (chidiff1*delta*(1.4739380748149558*eta1 + 0.06541707987699942*eta2 - 9.473290540936633*eta3) + chidiff2*(-0.3640838331639651*eta1 + 3.7369795937033756*eta2 - 8.709159662885131*eta3))*sqroot + (1.7335503724888923 + 12.656614578053683*eta1 - 139.6610487470118*eta2 + 456.78649322753824*eta3 - 599.2709938848282*eta4)*sqroot + chidiff1*delta*(2.3532739003216254*eta1 - 21.37216554136868*eta2 + 53.35003268489743*eta3)*chiPN1*sqroot + chiPN1*(-0.15782329022461472*(6.0309399412954345 - 229.16361598098678*eta1 + 3777.477006415653*eta2 - 31109.307191210424*eta3 + 139319.8239886073*eta4 - 324891.4001578353*eta5 + 307714.3954026392*eta6) - 0.03050157254864058*(4.232861441291087 + 1609.4251694451375*eta1 - 51213.27604422822*eta2 + 612317.1751155312*eta3 - 3.5589766538499263e6*eta4 + 1.0147654212772278e7*eta5 - 1.138861230369246e7*eta6)*chiPN1 - 0.026407497690308382*(-17.184685557542196 + 744.4743953122965*eta1 - 10494.512487701073*eta2 + 66150.52694069289*eta3 - 184787.79377504133*eta4 + 148102.4257785174*eta5 + 128167.89151782403*eta6)*chiPN2)*sqroot
        end
        # ok


        ### Collocation points INTERMEDIATE AMPLITUDE
        if ell_emm == 21
            int1 = abs(delta*eta1*(chidiff2*(5.159755997682368*eta1 - 30.293198248154948*eta2 + 63.70715919820867*eta3) + chidiff1*(8.262642080222694*eta1 - 415.88826990259116*eta2 + 1427.5951158851076*eta3)) + delta*eta1*(18.55363583212328 - 66.46950491124205*eta1 + 447.2214642597892*eta2 - 1614.178472020212*eta3 + 2199.614895727586*eta4) + chidiff1*eta5*(-1698.841763891122 - 195.27885562092342*S1 - 1.3098861736238572*S2) + delta*eta1*(chidiff1*(34.17829404207186*eta1 - 386.34587928670015*eta2 + 1022.8553774274128*eta3)*S1 + chidiff1*(56.76554600963724*eta1 - 491.4593694689354*eta2 + 1016.6019654342113*eta3)*S2) + delta*eta1*S1*(-8.276366844994188*(1.0677538075697492 - 24.12941323757896*eta1 + 516.7886322104276*eta2 - 4389.799658723288*eta3 + 16770.447637953577*eta4 - 23896.392706809565*eta5) - 1.6908277400304084*(3.4799140066657928 - 29.00026389706585*eta1 + 114.8330693231833*eta2 - 184.13091281984674*eta3 + 592.300353344717*eta4 - 2085.0821513466053*eta5)*S1 - 0.46006975902558517*(-2.1663474937625975 + 826.026625945615*eta1 - 17333.549622759732*eta2 + 142904.08962903373*eta3 - 528521.6231015554*eta4 + 731179.456702448*eta5)*S2))

            int2 = abs(delta*eta1*(13.757856231617446 - 12.783698329428516*eta1 + 12.048194546899204*eta2) + chidiff1*delta*eta1*(15.107530092096438*eta1 - 416.811753638553*eta2 + 1333.6181181686939*eta3) + chidiff1*eta5*(-1549.6199518612063 - 102.34716990474509*S1 - 3.3637011939285015*S2) + delta*eta1*(chidiff1*(36.358142200869295*eta1 - 384.2123173145321*eta2 + 984.6826660818275*eta3)*S1 + chidiff1*(4.159271594881928*eta1 + 105.10911749116399*eta2 - 639.190132707115*eta3)*S2) + delta*eta1*S1*(-8.097876227116853*(0.6569459700232806 + 9.861355377849485*eta1 - 116.88834714736281*eta2 + 593.8035334117192*eta3 - 1063.0692862578455*eta4) - 1.0546375154878165*(0.745557030602097 + 65.25215540635162*eta1 - 902.5751736558435*eta2 + 4350.442990924205*eta3 - 7141.611333893155*eta4)*S1 - 0.5006664599166409*(10.289020582277626 - 212.00728173197498*eta1 + 2334.0029399672358*eta2 - 11939.621138801092*eta3 + 21974.8201355744*eta4)*S2))
            
            int3 = abs(delta*eta1*(13.318990196097973 - 21.755549987331054*eta1 + 76.14884211156267*eta2 - 127.62161159798488*eta3) + chidiff1*delta*eta1*(17.704321326939414*eta1 - 434.4390350012534*eta2 + 1366.2408490833282*eta3) + chidiff1*delta*eta1*(11.877985158418596*eta1 - 131.04937626836355*eta2 + 343.79587860999874*eta3)*S1 + chidiff1*eta5*(-1522.8543551416456 - 16.639896279650678*S1 + 3.0053086651515843*S2) + delta*eta1*S1*(-8.665646058245033*(0.7862132291286934 + 8.293609541933655*eta1 - 111.70764910503321*eta2 + 576.7172598056907*eta3 - 1001.2370065269745*eta4) - 0.9459820574514348*(1.309016452198605 + 48.94077040282239*eta1 - 817.7854010574645*eta2 + 4331.56002883546*eta3 - 7518.309520232795*eta4)*S1 - 0.4308267743835775*(9.970654092010587 - 302.9708323417439*eta1 + 3662.099161055873*eta2 - 17712.883990278668*eta3 + 29480.158198408903*eta4)*S2))
            
            int4 = abs(delta*eta1*(13.094382343446163 - 22.831152256559523*eta1 + 83.20619262213437*eta2 - 139.25546924151664*eta3) + chidiff1*delta*eta1*(20.120192352555357*eta1 - 458.2592421214168*eta2 + 1430.3698681181*eta3) + chidiff1*delta*eta1*(12.925363020014743*eta1 - 126.87194512915104*eta2 + 280.6003655502327*eta3)*S1 + chidiff1*eta5*(-1528.956015503355 + 74.44462583487345*S1 - 2.2456928156392197*S2) + delta*eta1*S1*(-9.499741513411829*(0.912120958549489 + 2.400945118514037*eta1 - 33.651192908287236*eta2 + 166.04254881175257*eta3 - 248.5050377498615*eta4) - 0.7850652143322492*(1.534131218043425 + 60.81773903539479*eta1 - 1032.1319480683567*eta2 + 5381.481380750608*eta3 - 9077.037917192794*eta4)*S1 - 0.21540359093306097*(9.42805409480658 - 109.06544597367301*eta1 + 385.8345793110262*eta2 + 1889.9613367802453*eta3 - 9835.416414460055*eta4)*S2))
            
        elseif ell_emm == 33
            int1 = chidiff1*delta*eta1*(-0.3516244197696068*eta1 + 40.425151307421416*eta2 - 148.3162618111991*eta3) + delta*eta1*(26.998512565991778 - 146.29035440932105*eta1 + 914.5350366065115*eta2 - 3047.513201789169*eta3 + 3996.417635728702*eta4) + chidiff1*delta*eta1*(5.575274516197629*eta1 - 44.592719238427094*eta2 + 99.91399033058927*eta3)*S1 + delta*eta1*S1*(-0.5383304368673182*(-7.456619067234563 + 129.36947401891433*eta1 - 843.7897535238325*eta2 + 3507.3655567272644*eta3 - 9675.194644814854*eta4 + 11959.83533107835*eta5) - 0.28042799223829407*(-6.212827413930676 + 266.69059813274475*eta1 - 4241.537539226717*eta2 + 32634.43965039936*eta3 - 119209.70783201039*eta4 + 166056.27237509796*eta5)*S1) + chidiff1*eta5*(199.6863414922219 + 53.36849263931051*S1 + 7.650565415855383*S2)
        
            int2 = delta*eta1*(17.42562079069636 - 28.970875603981295*eta1 + 50.726220750178435*eta2) + chidiff1*delta*eta1*(-7.861956897615623*eta1 + 93.45476935080045*eta2 - 273.1170921735085*eta3) + chidiff1*delta*eta1*(-0.3265505633310564*eta1 - 9.861644053348053*eta2 + 60.38649425562178*eta3)*S1 + chidiff1*eta5*(234.13476431269862 + 51.2153901931183*S1 - 10.05114600643587*S2) + delta*eta1*S1*(0.3104472390387834*(6.073591341439855 + 169.85423386969634*eta1 - 4964.199967099143*eta2 + 42566.59565666228*eta3 - 154255.3408672655*eta4 + 205525.13910847943*eta5) + 0.2295327944679772*(19.236275867648594 - 354.7914372697625*eta1 + 1876.408148917458*eta2 + 2404.4151687877525*eta3 - 41567.07396803811*eta4 + 79210.33893514868*eta5)*S1 + 0.30983324991828787*(11.302200127272357 - 719.9854052004307*eta1 + 13278.047199998868*eta2 - 104863.50453518033*eta3 + 376409.2335857397*eta4 - 504089.07690692553*eta5)*S2)

            int3 = delta*eta1*(14.555522136327964 - 12.799844096694798*eta1 + 16.79500349318081*eta2) + chidiff1*delta*eta1*(-16.292654447108134*eta1 + 190.3516012682791*eta2 - 562.0936797781519*eta3) + chidiff1*delta*eta1*(-7.048898856045782*eta1 + 49.941617405768135*eta2 - 73.62033985436068*eta3)*S1 + chidiff1*eta5*(263.5151703818307 + 44.408527093031566*S1 + 10.457035444964653*S2) + delta*eta1*S1*(0.4590550434774332*(3.0594364612798635 + 207.74562213604057*eta1 - 5545.0086137386525*eta2 + 50003.94075934942*eta3 - 195187.55422847517*eta4 + 282064.174913521*eta5) + 0.657748992123043*(5.57939137343977 - 124.06189543062042*eta1 + 1276.6209573025596*eta2 - 6999.7659193505915*eta3 + 19714.675715229736*eta4 - 20879.999628681435*eta5)*S1 + 0.3695850566805098*(6.077183107132255 - 498.95526910874986*eta1 + 10426.348944657859*eta2 - 91096.64982858274*eta3 + 360950.6686625352*eta4 - 534437.8832860565*eta5)*S2)
        
            int4 = delta*eta1*(13.312095699772305 - 7.449975618083432*eta1 + 17.098576301150125*eta2) + delta*eta1*(chidiff1*(-31.171150896110156*eta1 + 371.1389274783572*eta2 - 1103.1917047361735*eta3) + chidiff2*(32.78644599730888*eta1 - 395.15713118955387*eta2 + 1164.9282236341376*eta3)) + chidiff1*delta*eta1*(-46.85669289852532*eta1 + 522.3965959942979*eta2 - 1485.5134187612182*eta3)*S1 + chidiff1*eta5*(287.90444670305715 - 21.102665129433042*chidiff2 + 7.635582066682054*S1 - 29.471275170013012*S2) + delta*eta1*S1*(0.6893003654021495*(3.1014226377197027 - 44.83989278653052*eta1 + 565.3767256471909*eta2 - 4797.429130246123*eta3 + 19514.812242035154*eta4 - 27679.226582207506*eta5) + 0.7068016563068026*(4.071212304920691 - 118.51094098279343*eta1 + 1788.1730303291356*eta2 - 13485.270489656365*eta3 + 48603.96661003743*eta4 - 65658.74746265226*eta5)*S1 + 0.2181399561677432*(-1.6754158383043574 + 303.9394443302189*eta1 - 6857.936471898544*eta2 + 59288.71069769708*eta3 - 216137.90827404748*eta4 + 277256.38289831823*eta5)*S2)
            
        elseif ell_emm == 32

            int1 = (chidiff2*(-0.2341404256829785*eta1 + 2.606326837996192*eta2 - 8.68296921440857*eta3) + chidiff1*delta*(0.5454562486736877*eta1 - 25.19759222940851*eta2 + 73.40268975811729*eta3))*sqroot + chidiff1*delta*(0.4422257616009941*eta1 - 8.490112284851655*eta2 + 32.22238925527844*eta3)*chiPN1*sqroot + chiPN1*(0.7067243321652764*(0.12885110296881636 + 9.608999847549535*eta1 - 85.46581740280585*eta2 + 325.71940024255775*eta3 + 175.4194342269804*eta4 - 1929.9084724384807*eta5) + 0.1540566313813899*(-0.3261041495083288 + 45.55785402900492*eta1 - 827.591235943271*eta2 + 7184.647314370326*eta3 - 28804.241518798244*eta4 + 43309.69769878964*eta5)*chiPN1)*sqroot + (480.0434256230109*eta1 + 25346.341240810478*eta2 - 99873.4707358776*eta3 + 106683.98302194536*eta4)*sqroot*(1 + 1082.6574834474493*eta1 + 10083.297670051445*eta2)^(-1)
        
            int2 = eta1*(chidiff2*(-4.175680729484314*eta1 + 47.54281549129226*eta2 - 128.88334273588077*eta3) + chidiff1*delta*(-0.18274358639599947*eta1 - 71.01128541687838*eta2 + 208.07105580635888*eta3)) + eta1*(4.760999387359598 - 38.57900689641654*eta1 + 456.2188780552874*eta2 - 4544.076411013166*eta3 + 24956.9592553473*eta4 - 69430.10468748478*eta5 + 77839.74180254337*eta6) + chidiff1*delta*eta1*(1.2198776533959694*eta1 - 26.816651899746475*eta2 + 68.72798751937934*eta3)*S1 + eta1*S1*(1.5098291294292217*(0.4844667556328104 + 9.848766999273414*eta1 - 143.66427232396376*eta2 + 856.9917885742416*eta3 - 1633.3295758142904*eta4) + 0.32413108737204144*(2.835358206961064 - 62.37317183581803*eta1 + 761.6103793011912*eta2 - 3811.5047139343505*eta3 + 6660.304740652403*eta4)*S1)
        
            int3 = 3.881450518842405*eta1 - 12.580316392558837*eta2 + 1.7262466525848588*eta3 + chidiff2*(-7.065118823041031*eta2 + 77.97950589523865*eta3 - 203.65975422378446*eta4) - 58.408542930248046*eta4 + chidiff1*delta*(1.924723094787216*eta2 - 90.92716917757797*eta3 + 387.00162600306226*eta4) + 403.5748987560612*eta5 + chidiff1*delta*(-0.2566958540737833*eta2 + 14.488550203412675*eta3 - 26.46699529970884*eta4)*chiPN1 + chiPN1*(0.3650871458400108*(71.57390929624825*eta2 - 994.5272351916166*eta3 + 6734.058809060536*eta4 - 18580.859291282686*eta5 + 16001.318492586077*eta6) + 0.0960146077440495*(451.74917589707513*eta2 - 9719.470997418284*eta3 + 83403.5743434538*eta4 - 318877.43061174755*eta5 + 451546.88775684836*eta6)*chiPN1 - 0.03985156529181297*(-304.92981902871617*eta2 + 3614.518459296278*eta3 - 7859.4784979916085*eta4 - 46454.57664737511*eta5 + 162398.81483375572*eta6)*chiPN2)
        
            int4 = eta1*(chidiff2*(-8.572797326909152*eta1 + 92.95723645687826*eta2 - 236.2438921965621*eta3) + chidiff1*delta*(6.674358856924571*eta1 - 171.4826985994883*eta2 + 645.2760206304703*eta3)) + eta1*(3.921660532875504 - 16.57299637423352*eta1 + 25.254017911686333*eta2 - 143.41033155133266*eta3 + 692.926425981414*eta4) + chidiff1*delta*eta1*(-3.582040878719185*eta1 + 57.75888914133383*eta2 - 144.21651114700492*eta3)*S1 + eta1*S1*(1.242750265695504*(-0.522172424518215 + 25.168480118950065*eta1 - 303.5223688400309*eta2 + 1858.1518762309654*eta3 - 3797.3561904195085*eta4) + 0.2927045241764365*(0.5056957789079993 - 15.488754837330958*eta1 + 471.64047356915603*eta2 - 3131.5783196211587*eta3 + 6097.887891566872*eta4)*S1)
            
            
            
        elseif ell_emm == 44
            
            int1 = eta1*(chidiff1*delta*(1.5378890240544967*eta1 - 3.4499418893734903*eta2 + 16.879953490422782*eta3) + chidiff2*(1.720226708214248*eta1 - 11.87925165364241*eta2 + 23.259283336239545*eta3)) + eta1*(8.790173464969538 - 64.95499142822892*eta1 + 324.1998823562892*eta2 - 1111.9864921907126*eta3 + 1575.602443847111*eta4) + eta1*S1*(-0.062333275821238224*(-21.630297087123807 + 137.4395894877131*eta1 + 64.92115530780129*eta2 - 1013.1110639471394*eta3) - 0.11014697070998722*(4.149721483857751 - 108.6912882442823*eta1 + 831.6073263887092*eta2 - 1828.2527520190122*eta3)*S1 - 0.07704777584463054*(4.581767671445529 - 50.35070009227704*eta1 + 344.9177692251726*eta2 - 858.9168637051405*eta3)*S2)
        
            int2 = eta1*(chidiff1*delta*(2.3123974306694057*eta1 - 12.237594841284904*eta2 + 44.78225529547671*eta3) + chidiff2*(2.9282931698944292*eta1 - 25.624210264341933*eta2 + 61.05270871360041*eta3)) + eta1*(6.98072197826729 - 46.81443520117986*eta1 + 236.76146303619544*eta2 - 920.358408667518*eta3 + 1478.050456337336*eta4) + eta1*S1*(-0.07801583359561987*(-28.29972282146242 + 752.1603553640072*eta1 - 10671.072606753183*eta2 + 83447.0461509547*eta3 - 350025.2112501252*eta4 + 760889.6919776166*eta5 - 702172.2934567826*eta6) + 0.013159545629626014*(91.1469833190294 - 3557.5003799977294*eta1 + 52391.684517955284*eta2 - 344254.9973814295*eta3 + 1.0141877915334814e6*eta4 - 1.1505186449682908e6*eta5 + 268756.85659532435*eta6)*S1)
            
            int3 = eta1*(chidiff1*delta*(-0.8765502142143329*eta1 + 22.806632458441996*eta2 - 43.675503209991184*eta3) + chidiff2*(0.48698617426180074*eta1 - 4.302527065360426*eta2 + 16.18571810759235*eta3)) + eta1*(6.379772583015967 - 44.10631039734796*eta1 + 269.44092930942793*eta2 - 1285.7635006711453*eta3 + 2379.538739132234*eta4) + eta1*chiPN1*(-0.23316184683282615*(-1.7279023138971559 - 23.606399143993716*eta1 + 409.3387618483284*eta2 - 1115.4147472977265*eta3) - 0.09653777612560172*(-5.310643306559746 - 2.1852511802701264*eta1 + 541.1248219096527*eta2 - 1815.7529908827103*eta3)*chiPN1 - 0.060477799540741804*(-14.578189130145661 + 175.6116682068523*eta1 - 569.4799973930861*eta2 + 426.0861915646515*eta3)*chiPN2)
        
            int4 = eta1*(chidiff1*delta*(-2.461738962276138*eta1 + 45.3240543970684*eta2 - 112.2714974622516*eta3) + chidiff2*(0.9158352037567031*eta1 - 8.724582331021695*eta2 + 28.44633544874233*eta3)) + eta1*(6.098676337298138 - 45.42463610529546*eta1 + 350.97192927929433*eta2 - 2002.2013283876834*eta3 + 4067.1685640401033*eta4) + eta1*chiPN1*(-0.36068516166901304*(-2.120354236840677 - 47.56175350408845*eta1 + 1618.4222330016048*eta2 - 14925.514654896673*eta3 + 60287.45399959349*eta4 - 91269.3745059139*eta5) - 0.09635801207669747*(-11.824692837267394 + 371.7551657959369*eta1 - 4176.398139238679*eta2 + 16655.87939259747*eta3 - 4102.218189945819*eta4 - 67024.98285179552*eta5)*chiPN1 - 0.06565232123453196*(-26.15227471380236 + 1869.0168486099005*eta1 - 33951.35186039629*eta2 + 253694.6032002248*eta3 - 845341.6001856657*eta4 + 1.0442282862506858e6*eta5)*chiPN2)
        
        end 

        ### Collocation points RINGDOWN AMPLITUDE
        # from IMRPhenomXHM_RD_Amp_Coefficients(pWF22, pWFHM, pAmp);


        if ell_emm == 21
            alambda = abs(delta*(0.24548180919287976 - 0.25565119457386487*eta1)*eta1 + chidiff1*delta*eta1*(0.5670798742968471*eta1 - 14.276514548218454*eta2 + 45.014547333879136*eta3) + chidiff1*delta*eta1*(0.4580805242442763*eta1 - 4.859294663135058*eta2 + 14.995447609839573*eta3)*chiPN1 + chidiff1*eta5*(-27.031582936285528 + 6.468164760468401*chiPN1 + 0.34222101136488015*chiPN2) + delta*eta1*chiPN1*(-0.2204878224611389*(1.0730799832007898 - 3.44643820338605*eta1 + 32.520429274459836*eta2 - 83.21097158567372*eta3) + 0.008901444811471891*(-5.876973170072921 + 120.70115519895002*eta1 - 916.5281661566283*eta2 + 2306.8425350489847*eta3)*chiPN1 + 0.015541783867953005*(2.4780170686140455 + 17.377013149762398*eta1 - 380.91157168170236*eta2 + 1227.5332509075172*eta3)*chiPN2))

            rdcp1 = abs(delta*eta1*(12.880905080761432 - 23.5291063016996*eta1 + 92.6090002736012*eta2 - 175.16681482428694*eta3) + chidiff1*delta*eta1*(26.89427230731867*eta1 - 710.8871223808559*eta2 + 2255.040486907459*eta3) + chidiff1*delta*eta1*(21.402708785047853*eta1 - 232.07306353130417*eta2 + 591.1097623278739*eta3)*chiPN1 + delta*eta1*chiPN1*(-10.090867481062709*(0.9580746052260011 + 5.388149112485179*eta1 - 107.22993216128548*eta2 + 801.3948756800821*eta3 - 2688.211889175019*eta4 + 3950.7894052628735*eta5 - 1992.9074348833092*eta6) - 0.42972412296628143*(1.9193131231064235 + 139.73149069609775*eta1 - 1616.9974609915555*eta2 - 3176.4950303461164*eta3 + 107980.65459735804*eta4 - 479649.75188253267*eta5 + 658866.0983367155*eta6)*chiPN1) + chidiff1*eta5*(-1512.439342647443 + 175.59081294852444*chiPN1 + 10.13490934572329*chiPN2))

            rdcp2 = abs(delta*(9.112452928978168 - 7.5304766811877455*eta1)*eta1 + chidiff1*delta*eta1*(16.236533863306132*eta1 - 500.11964987628926*eta2 + 1618.0818430353293*eta3) + chidiff1*delta*eta1*(2.7866868976718226*eta1 - 0.4210629980868266*eta2 - 20.274691328125606*eta3)*chiPN1 + chidiff1*eta5*(-1116.4039232324135 + 245.73200219767514*chiPN1 + 21.159179960295855*chiPN2) + delta*eta1*chiPN1*(-8.236485576091717*(0.8917610178208336 + 5.1501231412520285*eta1 - 87.05136337926156*eta2 + 519.0146702141192*eta3 - 997.6961311502365*eta4) + 0.2836840678615208*(-0.19281297100324718 - 57.65586769647737*eta1 + 586.7942442434971*eta2 - 1882.2040277496196*eta3 + 2330.3534917059906*eta4)*chiPN1 + 0.40226131643223145*(-3.834742668014861 + 190.42214703482531*eta1 - 2885.5110686004946*eta2 + 16087.433824017446*eta3 - 29331.524552164105*eta4)*chiPN2))

            rdcp3 = abs(delta*(2.920930733198033 - 3.038523690239521*eta1)*eta1 + chidiff1*delta*eta1*(6.3472251472354975*eta1 - 171.23657247338042*eta2 + 544.1978232314333*eta3) + chidiff1*delta*eta1*(1.9701247529688362*eta1 - 2.8616711550845575*eta2 - 0.7347258030219584*eta3)*chiPN1 + chidiff1*eta5*(-334.0969956136684 + 92.91301644484749*chiPN1 - 5.353399481074393*chiPN2) + delta*eta1*chiPN1*(-2.7294297839371824*(1.148166706456899 - 4.384077347340523*eta1 + 36.120093043420326*eta2 - 87.26454353763077*eta3) + 0.23949142867803436*(-0.6931516433988293 + 33.33372867559165*eta1 - 307.3404155231787*eta2 + 862.3123076782916*eta3)*chiPN1 + 0.1930861073906724*(3.7735099269174106 - 19.11543562444476*eta1 - 78.07256429516346*eta2 + 485.67801863289293*eta3)*chiPN2))

            lambda = abs(1.0092933052569775 - 0.2791855444800297*eta1 + 1.7110615047319937*eta2 + chidiff1*(-0.1054835719277311*eta1 + 7.506083919925026*eta2 - 30.1595680078279*eta3) + chidiff1*(2.078267611384239*eta1 - 10.166026002515457*eta2 - 1.2091616330625208*eta3)*chiPN1 + chiPN1*(0.17250873250247642*(1.0170226856985174 + 1.0395650952176598*eta1 - 35.73623734051525*eta2 + 403.68074286921444*eta3 - 1194.6152711219886*eta4) + 0.06850746964805364*(1.507796537056924 + 37.81075363806507*eta1 - 863.117144661059*eta2 + 6429.543634627373*eta3 - 15108.557419182316*eta4)*chiPN1))

            sigma = abs(1.374451177213076 - 0.1147381625630186*eta1 + chidiff1*(0.6646459256372743*eta1 - 5.020585319906719*eta2 + 9.817281653770431*eta3) + chidiff1*(3.8734254747587973*eta1 - 39.880716190740465*eta2 + 99.05511583518896*eta3)*S1 + S1*(0.013272603498067647*(1.809972721953344 - 12.560287006325837*eta1 - 134.597005438578*eta2 + 786.2235720637008*eta3) + 0.006850483944311038*(-6.478737679813189 - 200.29813775611166*eta1 + 2744.3629484255357*eta2 - 7612.096007280672*eta3)*S1))
        
        elseif ell_emm == 33

            alambda = 0.
            lambda = 0.
            sigma = 0.  #  the values put to 0. are not defined in the code for the version 122022

            rdcp1 = abs(delta*eta1*(12.439702602599235 - 4.436329538596615*eta1 + 22.780673360839497*eta2) + delta*eta1*(chidiff1*(-41.04442169938298*eta1 + 502.9246970179746*eta2 - 1524.2981907688634*eta3) + chidiff2*(32.23960072974939*eta1 - 365.1526474476759*eta2 + 1020.6734178547847*eta3)) + chidiff1*delta*eta1*(-52.85961155799673*eta1 + 577.6347407795782*eta2 - 1653.496174539196*eta3)*S1 + chidiff1*eta5*(257.33227387984863 - 34.5074027042393*chidiff2 - 21.836905132600755*S1 - 15.81624534976308*S2) + 13.499999999999998*delta*eta1*S1*(-0.13654149379906394*(2.719687834084113 + 29.023992126142304*eta1 - 742.1357702210267*eta2 + 4142.974510926698*eta3 - 6167.08766058184*eta4 - 3591.1757995710486*eta5) - 0.06248535354306988*(6.697567446351289 - 78.23231700361792*eta1 + 444.79350113344543*eta2 - 1907.008984765889*eta3 + 6601.918552659412*eta4 - 10056.98422430965*eta5)*S1)*(-3.9329308614837704 + S1)^(-1))

            rdcp2 = abs(delta*eta1*(8.425057692276933 + 4.543696144846763*eta1) + chidiff1*delta*eta1*(-32.18860840414171*eta1 + 412.07321398189293*eta2 - 1293.422289802462*eta3) + chidiff1*delta*eta1*(-17.18006888428382*eta1 + 190.73514518113845*eta2 - 636.4802385540647*eta3)*S1 + delta*eta1*S1*(0.1206817303851239*(8.667503604073314 - 144.08062755162752*eta1 + 3188.189172446398*eta2 - 35378.156133055556*eta3 + 163644.2192178668*eta4 - 265581.70142471837*eta5) + 0.08028332044013944*(12.632478544060636 - 322.95832000179297*eta1 + 4777.45310151897*eta2 - 35625.58409457366*eta3 + 121293.97832549023*eta4 - 148782.33687815256*eta5)*S1) + chidiff1*eta5*(159.72371180117415 - 29.10412708633528*chidiff2 - 1.873799747678187*S1 + 41.321480132899524*S2))

            rdcp3 = abs(delta*eta1*(2.485784720088995 + 2.321696430921996*eta1) + delta*eta1*(chidiff1*(-10.454376404653859*eta1 + 147.10344302665484*eta2 - 496.1564538739011*eta3) + chidiff2*(-5.9236399792925996*eta1 + 65.86115501723127*eta2 - 197.51205149250532*eta3)) + chidiff1*delta*eta1*(-10.27418232676514*eta1 + 136.5150165348149*eta2 - 473.30988537734174*eta3)*S1 + chidiff1*eta5*(32.07819766300362 - 3.071422453072518*chidiff2 + 35.09131921815571*S1 + 67.23189816732847*S2) + 13.499999999999998*delta*eta1*S1*(0.0011484326782460882*(4.1815722950796035 - 172.58816646768219*eta1 + 5709.239330076732*eta2 - 67368.27397765424*eta3 + 316864.0589150127*eta4 - 517034.11171277676*eta5) - 0.009496797093329243*(0.9233282181397624 - 118.35865186626413*eta1 + 2628.6024206791726*eta2 - 23464.64953722729*eta3 + 94309.57566199072*eta4 - 140089.40725211444*eta5)*S1)*(0.09549360183532198 - 0.41099904730526465*S1 + S2)^(-1))


        elseif ell_emm == 32

            alambda = chidiff2*(-3.4614418482110163*eta3 + 35.464117772624164*eta4 - 85.19723511005235*eta5) + chidiff1*delta*(2.0328561081997463*eta3 - 46.18751757691501*eta4 + 170.9266105597438*eta5) + chidiff2*(-0.4600401291210382*eta3 + 12.23450117663151*eta4 - 42.74689906831975*eta5)*S1 + chidiff1*delta*(5.786292428422767*eta3 - 53.60467819078566*eta4 + 117.66195692191727*eta5)*S1 + S1*(-0.0013330716557843666*(56.35538385647113*eta1 - 1218.1550992423377*eta2 + 16509.69605686402*eta3 - 102969.88022112886*eta4 + 252228.94931931415*eta5 - 150504.2927996263*eta6) + 0.0010126460331462495*(-33.87083889060834*eta1 + 502.6221651850776*eta2 - 1304.9210590188136*eta3 - 36980.079328277505*eta4 + 295469.28617550555*eta5 - 597155.7619486618*eta6)*S1 - 0.00043088431510840695*(-30.014415072587354*eta1 - 1900.5495690280086*eta2 + 76517.21042363928*eta3 - 870035.1394696251*eta4 + 3.9072674134789007e6*eta5 - 6.094089675611567e6*eta6)*S2) + (0.08408469319155859*eta1 - 1.223794846617597*eta2 + 6.5972460654253515*eta3 - 15.707327897569396*eta4 + 14.163264397061505*eta5)*(1 - 8.612447115134758*eta1 + 18.93655612952139*eta2)^(-1)

            lambda = 0.978510781593996 + 0.36457571743142897*eta1 - 12.259851752618998*eta2 + 49.19719473681921*eta3 + chidiff1*delta*(-188.37119473865533*eta3 + 2151.8731700399308*eta4 - 6328.182823770599*eta5) + chidiff2*(115.3689949926392*eta3 - 1159.8596972989067*eta4 + 2657.6998831179444*eta5) + S1*(0.22358643406992756*(0.48943645614341924 - 32.06682257944444*eta1 + 365.2485484044132*eta2 - 915.2489655397206*eta3) + 0.0792473022309144*(1.877251717679991 - 103.65639889587327*eta1 + 1202.174780792418*eta2 - 3206.340850767219*eta3)*S1)
        
            sigma = 1.3353917551819414 + 0.13401718687342024*eta1 + chidiff1*delta*(144.37065005786636*eta3 - 754.4085447486738*eta4 + 123.86194078913776*eta5) + chidiff2*(209.09202210427972*eta3 - 1769.4658099037918*eta4 + 3592.287297392387*eta5) + S1*(-0.012086025709597246*(-6.230497473791485 + 600.5968613752918*eta1 - 6606.1009717965735*eta2 + 17277.60594350428*eta3) - 0.06066548829900489*(-0.9208054306316676 + 142.0346574366267*eta1 - 1567.249168668069*eta2 + 4119.373703246675*eta3)*S1)

            rdaux1 = chidiff2*(-4.188795724777721*eta2 + 53.39200466700963*eta3 - 131.19660856923554*eta4) + chidiff1*delta*(14.284921364132623*eta2 - 321.26423637658746*eta3 + 1242.865584938088*eta4) + S1*(-0.022968727462555794*(83.66854837403105*eta1 - 3330.6261333413177*eta2 + 77424.12614733395*eta3 - 710313.3016672594*eta4 + 2.6934917075009225e6*eta5 - 3.572465179268999e6*eta6) + 0.0014795114305436387*(-1672.7273629876313*eta1 + 90877.38260964208*eta2 - 1.6690169155105734e6*eta3 + 1.3705532554135624e7*eta4 - 5.116110998398143e7*eta5 + 7.06066766311127e7*eta6)*S1) + (4.45156488896258*eta1 - 77.39303992494544*eta2 + 522.5070635563092*eta3 - 1642.3057499049708*eta4 + 2048.333892310575*eta5)*(1 - 9.611489164758915*eta1 + 24.249594730050312*eta2)^(-1)

            rdaux2 = chidiff2*(-18.550171209458394*eta2 + 188.99161055445936*eta3 - 440.26516625611*eta4) + chidiff1*delta*(13.132625215315063*eta2 - 340.5204040505528*eta3 + 1327.1224176812448*eta4) + chiPN1*(-0.16707403272774676*(6.678916447469937*eta1 + 1331.480396625797*eta2 - 41908.45179140144*eta3 + 520786.0225074669*eta4 - 3.1894624909922685e6*eta5 + 9.51553823212259e6*eta6 - 1.1006903622406831e7*eta7) + 0.015205286051218441*(108.10032279461095*eta1 - 16084.215590200103*eta2 + 462957.5593513407*eta3 - 5.635028227588545e6*eta4 + 3.379925277713386e7*eta5 - 9.865815275452062e7*eta6 + 1.1201307979786257e8*eta7)*S1) + (3.902154247490771*eta1 - 55.77521071924907*eta2 + 294.9496843041973*eta3 - 693.6803787318279*eta4 + 636.0141528226893*eta5)*(1 - 8.56699762573719*eta1 + 19.119341007236955*eta2)^(-1)

            rdcp1 = abs(chidiff2*(-261.63903838092017*eta3 + 2482.4929818200458*eta4 - 5662.765952006266*eta5) + chidiff1*delta*(200.3023530582654*eta3 - 3383.07742098347*eta4 + 11417.842708417566*eta5) + chidiff2*(-177.2481070662751*eta3 + 1820.8637746828358*eta4 - 4448.151940319403*eta5)*S1 + chidiff1*delta*(412.749304734278*eta3 - 4156.641392955615*eta4 + 10116.974216563232*eta5)*S1 + S1*(-0.07383539239633188*(40.59996146686051*eta1 - 527.5322650311067*eta2 + 4167.108061823492*eta3 - 13288.883172763119*eta4 - 23800.671572828596*eta5 + 146181.8016013141*eta6) + 0.03576631753501686*(-13.96758180764024*eta1 - 797.1235306450683*eta2 + 18007.56663810595*eta3 - 151803.40642097822*eta4 + 593811.4596071478*eta5 - 878123.747877138*eta6)*S1 + 0.01007493097350273*(-27.77590078264459*eta1 + 4011.1960424049857*eta2 - 152384.01804465035*eta3 + 1.7595145936445233e6*eta4 - 7.889230647117076e6*eta5 + 1.2172078072446395e7*eta6)*S2) + (4.146029818148087*eta1 - 61.060972560568054*eta2 + 336.3725848841942*eta3 - 832.785332776221*eta4 + 802.5027431944313*eta5)*(1 - 8.662174796705683*eta1 + 19.288918757536685*eta2)^(-1))

            rdcp2 = abs(chidiff2*(-220.42133216774002*eta3 + 2082.031407555522*eta4 - 4739.292554291661*eta5) + chidiff1*delta*(179.07548162694007*eta3 - 2878.2078963030094*eta4 + 9497.998559135678*eta5) + chidiff2*(-128.07917402087625*eta3 + 1392.4598433465628*eta4 - 3546.2644951338134*eta5)*S1 + chidiff1*delta*(384.31792882093424*eta3 - 3816.5687272960417*eta4 + 9235.479593415908*eta5)*S1 + S1*(-0.06144774696295017*(35.72693522898656*eta1 - 168.08433700852038*eta2 - 3010.678442066521*eta3 + 45110.034521934074*eta4 - 231569.4154711447*eta5 + 414234.84895584086*eta6) + 0.03663881822701642*(-22.057692852225696*eta1 + 223.9912685075838*eta2 - 1028.5261783449762*eta3 - 12761.957255385*eta4 + 141784.13567610556*eta5 - 328718.5349981628*eta6)*S1 + 0.004849853669413881*(-90.35491669965123*eta1 + 19286.158446325957*eta2 - 528138.5557827373*eta3 + 5.175061086459432e6*eta4 - 2.1142182400264673e7*eta5 + 3.0737963347449116e7*eta6)*S2) + (3.133378729082171*eta1 - 45.83572706555282*eta2 + 250.23275606463622*eta3 - 612.0498767005383*eta4 + 580.3574091493459*eta5)*(1 - 8.698032720488515*eta1 + 19.38621948411302*eta2)^(-1))

            rdcp3 = abs(chidiff2*(-79.14146757219045*eta3 + 748.8207876524461*eta4 - 1712.3401586150026*eta5) + chidiff1*delta*(65.1786095079065*eta3 - 996.4553252426255*eta4 + 3206.5675278160684*eta5) + chidiff2*(-36.474455088940225*eta3 + 421.8842792746865*eta4 - 1117.0227933265749*eta5)*S1 + chidiff1*delta*(169.07368933925878*eta3 - 1675.2562326502878*eta4 + 4040.0077967763787*eta5)*S1 + S1*(-0.01992370601225598*(36.307098892574196*eta1 - 846.997262853445*eta2 + 16033.60939445582*eta3 - 138800.53021166887*eta4 + 507922.88946543116*eta5 - 647376.1499824544*eta6) + 0.014207919520826501*(-33.80287899746716*eta1 + 1662.2913368534057*eta2 - 31688.885017467597*eta3 + 242813.43893659746*eta4 - 793178.4767168422*eta5 + 929016.897093022*eta6)*S1) + (0.9641853854287679*eta1 - 13.801372413989519*eta2 + 72.80610853168994*eta3 - 168.65551450831953*eta4 + 147.2372582604103*eta5)*(1 - 8.65963828355163*eta1 + 19.112920222001367*eta2)^(-1))

        
        elseif ell_emm == 44
            alambda = 0.
            lambda = 0.
            sigma = 0.

            rdcp1 = abs(eta1*(chidiff1*delta*(-8.51952446214978*eta1 + 117.76530248141987*eta2 - 297.2592736781142*eta3) + chidiff2*(-0.2750098647982238*eta1 + 4.456900599347149*eta2 - 8.017569928870929*eta3)) + eta1*(5.635069974807398 - 33.67252878543393*eta1 + 287.9418482197136*eta2 - 3514.3385364216438*eta3 + 25108.811524802128*eta4 - 98374.18361532023*eta5 + 158292.58792484726*eta6) + eta1*S1*(-0.4360849737360132*(-0.9543114627170375 - 58.70494649755802*eta1 + 1729.1839588870455*eta2 - 16718.425586396803*eta3 + 71236.86532610047*eta4 - 111910.71267453219*eta5) - 0.024861802943501172*(-52.25045490410733 + 1585.462602954658*eta1 - 15866.093368857853*eta2 + 35332.328181283*eta3 + 168937.32229060197*eta4 - 581776.5303770923*eta5)*S1 + 0.005856387555754387*(186.39698091707513 - 9560.410655118145*eta1 + 156431.3764198244*eta2 - 1.0461268207440731e6*eta3 + 3.054333578686424e6*eta4 - 3.2369858387064277e6*eta5)*S2))

            rdcp2 = abs(eta1*(chidiff1*delta*(-2.861653255976984*eta1 + 50.50227103211222*eta2 - 123.94152825700999*eta3) + chidiff2*(2.9415751419018865*eta1 - 28.79779545444817*eta2 + 72.40230240887851*eta3)) + eta1*(3.2461722686239307 + 25.15310593958783*eta1 - 792.0167314124681*eta2 + 7168.843978909433*eta3 - 30595.4993786313*eta4 + 49148.57065911245*eta5) + eta1*S1*(-0.23311779185707152*(-1.0795711755430002 - 20.12558747513885*eta1 + 1163.9107546486134*eta2 - 14672.23221502075*eta3 + 73397.72190288734*eta4 - 127148.27131388368*eta5) + 0.025805905356653*(11.929946153728276 + 350.93274421955806*eta1 - 14580.02701600596*eta2 + 174164.91607515427*eta3 - 819148.9390278616*eta4 + 1.3238624538095295e6*eta5)*S1 + 0.019740635678180102*(-7.046295936301379 + 1535.781942095697*eta1 - 27212.67022616794*eta2 + 201981.0743810629*eta3 - 696891.1349708183*eta4 + 910729.0219043035*eta5)*S2))

            rdcp3 = abs(eta1*(chidiff1*delta*(2.4286414692113816*eta1 - 23.213332913737403*eta2 + 66.58241012629095*eta3) + chidiff2*(3.085167288859442*eta1 - 31.60440418701438*eta2 + 78.49621016381445*eta3)) + eta1*(0.861883217178703 + 13.695204704208976*eta1 - 337.70598252897696*eta2 + 2932.3415281149432*eta3 - 12028.786386004691*eta4 + 18536.937955014455*eta5) + eta1*S1*(-0.048465588779596405*(-0.34041762314288154 - 81.33156665674845*eta1 + 1744.329802302927*eta2 - 16522.343895064576*eta3 + 76620.18243090731*eta4 - 133340.93723954144*eta5) + 0.024804027856323612*(-8.666095805675418 + 711.8727878341302*eta1 - 13644.988225595187*eta2 + 112832.04975245205*eta3 - 422282.0368440555*eta4 + 584744.0406581408*eta5)*S1))

        end


        ##### PHASE FITS ####

        if ell_emm == 21

            #Insp_Phase_21_lambda = 13.664473636545068 - 170.08866400251395*eta + 3535.657736681598*eta2 - 26847.690494515424*eta3 + 96463.68163125668*eta4 - 133820.89317471132*eta5+(S*(18.52571430563905 - 41.55066592130464*S + eta3*(83493.24265292779 + 16501.749243703132*S - 149700.4915210766*S2) + eta*(3642.5891077598003 + 1198.4163078715173*S - 6961.484805326852*S2) + 33.8697137964237*S2 + eta2*(-35031.361998480075 - 7233.191207000735*S + 62149.00902591944*S2)))/(6.880288191574696 + 1. .*S)+-134.27742343186577*dchi*Seta*eta2

            p1 = 4045.84 + 7.63226/eta - 1956.93*eta - 23428.1*eta2 + 369153. .*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6+- 347.273*S + 83.5428*S2 - 355.67*S3 + (4.44457*S + 16.5548*S2 + 13.6971*S3)/eta + eta*( - 79.761*S - 355.299*S2 + 1114.51*S3 - 1077.75*S4) + 92.6654*S4 + eta2*(- 619.837*S - 722.787*S2 + 2392.73*S3 + 2689.18*S4)+ ( 918.976*chi1*Seta - 918.976*chi2*Seta)*eta + ( 91.7679*chi1*Seta - 91.7679*chi2*Seta)*eta2
            p2 = 3509.09 + 0.91868/eta + 194.72*eta - 27556.2*eta2 + 369153. .*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6+((0.7083999999999999 - 60.1611*eta + 131.815*eta2 - 619.837*eta3)*S + (6.104720000000001 - 59.2068*eta + 278.588*eta2 - 722.787*eta3)*S2 + (5.7791 + 117.913*eta - 1180.4*eta2 + 2392.73*eta3)*S3 + eta*(92.6654 - 1077.75*eta + 2689.18*eta2)*S4)/(eta)+-91.7679*Seta*eta*(chi1*(-1.6012352903357276 - 1. .*eta) + chi2*(1.6012352903357276 + 1. .*eta))
            p3 = 3241.68 + 890.016*eta - 28651.9*eta2 + 369153. .*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6+(-2.2484 + 187.641*eta - 619.837*eta2)*S + (3.22603 + 166.323*eta - 722.787*eta2)*S2 + (117.913 - 1094.59*eta + 2392.73*eta2)*S3 + (92.6654 - 1077.75*eta + 2689.18*eta2)*S4+91.7679*(dchi)*delta*eta2
            p4 = 3160.88 + 974.355*eta - 28932.5*eta2 + 369780. .*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6+(26.3355 - 196.851*eta + 438.401*eta2)*S + (45.9957 - 256.248*eta + 117.563*eta2)*S2 + (-20.0261 + 467.057*eta - 1613. .*eta2)*S3 + (-61.7446 + 577.057*eta - 1096.81*eta2)*S4+65.3326*dchi*Seta*eta2
            p5 = 3102.36 + 315.911*eta - 1688.26*eta2 + 3635.76*eta3+(-23.0959 + 320.93*eta - 1029.76*eta2)*S + (-49.5435 + 826.816*eta - 3079.39*eta2)*S2 + (40.7054 - 365.842*eta + 1094.11*eta2)*S3 + (81.8379 - 1243.26*eta + 4689.22*eta2)*S4+119.014*(dchi)*Seta*eta2
            p6 = 3089.18 + 4.89194*eta + 190.008*eta2 - 255.245*eta3+(2.96997 + 57.1612*eta - 432.223*eta2)*S + (-18.8929 + 630.516*eta - 2804.66*eta2)*S2 + (-24.6193 + 549.085*eta2)*S3 + (-12.8798 - 722.674*eta + 3967.43*eta2)*S4+74.0984*(dchi)*Seta*eta2
            
        elseif ell_emm == 33 
            
            #Insp_Phase_33_lambda = 4.1138398568400705 + 9.772510519809892*eta - 103.92956504520747*eta2 + 242.3428625556764*eta3+((-0.13253553909611435 + 26.644159828590055*eta - 105.09339163109497*eta2)*S)/(1. + 0.11322426762297967*S)+-19.705359163581168*dchi*eta2*delta

            p1 = 4360.19 + 4.27128/eta - 8727.4*eta + 18485.9*eta2 + 371303.00000000006*eta3 - 3.22792e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6+((11.6635 - 251.579*eta - 3255.6400000000003*eta2 + 19614.6*eta3 - 34860.2*eta4)*S + (14.8017 + 204.025*eta - 5421.92*eta2 + 36587.3*eta3 - 74299.5*eta4)*S2)/(eta)+eta*(223.65100000000004*chi1*Seta*(3.9201300240106223 + 1. .*eta) - 223.65100000000004*chi2*Seta*(3.9201300240106223 + 1. .*eta))
            p2 = 3797.06 + 0.786684/eta - 2397.09*eta - 25514. .*eta2 + 518314.99999999994*eta3 - 3.41708e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6+((6.7812399999999995 + 39.4668*eta - 3520.37*eta2 + 19614.6*eta3 - 34860.2*eta4)*S + (4.80384 + 293.215*eta - 5914.61*eta2 + 36587.3*eta3 - 74299.5*eta4)*S2)/(eta)+-223.65100000000004*Seta*eta*(chi1*(-1.3095134830606614 - 1. .*eta) + chi2*(1.3095134830606614 + 1. .*eta))
            p3 = 3321.83 + 1796.03*eta - 52406.1*eta2 + 605028. .*eta3 - 3.52532e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6+(223.601 - 3714.77*eta + 19614.6*eta2 - 34860.2*eta3)*S + (314.317 - 5906.46*eta + 36587.3*eta2 - 74299.5*eta3)*S2+223.651*(dchi)*Seta*eta2
            p4 = 3239.44 - 661.15*eta + 5139.79*eta2 + 3456.2*eta3 - 248477. .*eta4 + 1.17255e6*eta5 - 1.70363e6*eta6+(225.859 - 4150.09*eta + 24364. .*eta2 - 46537.3*eta3)*S + (35.2439 - 994.971*eta + 8953.98*eta2 - 23603.5*eta3)*S2 + (-310.489 + 5946.15*eta - 35337.1*eta2 + 67102.4*eta3)*S3+30.484*dchi*Seta*eta2
            p5 = 3114.3 + 2143.06*eta - 49428.3*eta2 + 563997. .*eta3 - 3.35991e6*eta4 + 9.99745e6*eta5 - 1.17123e7*eta6+(190.051 - 3705.08*eta + 23046.2*eta2 - 46537.3*eta3)*S + (63.6615 - 1414.2*eta + 10166.1*eta2 - 23603.5*eta3)*S2 + (-257.524 + 5179.97*eta - 33001.4*eta2 + 67102.4*eta3)*S3+54.9833*dchi*Seta*eta2
            p6 = 3111.46 + 384.121*eta - 13003.6*eta2 + 179537. .*eta3 - 1.19313e6*eta4 + 3.79886e6*eta5 - 4.64858e6*eta6+(182.864 - 3834.22*eta + 24532.9*eta2 - 50165.9*eta3)*S + (21.0158 - 746.957*eta + 6701.33*eta2 - 17842.3*eta3)*S2 + (-292.855 + 5886.62*eta - 37382.4*eta2 + 75501.8*eta3)*S3+75.5162*dchi*Seta*eta2
            
        elseif ell_emm == 32  
            
            #Insp_Phase_32_lambda = (9.913819875501506 + 18.424900617803107*eta - 574.8672384388947*eta2 + 2671.7813055097877*eta3 - 6244.001932443913*eta4)/(1. - 0.9103118343073325*eta)+(-4.367632806613781 + 245.06757304950986*eta - 2233.9319708029775*eta2 + 5894.355429022858*eta3)*S + (-1.375112297530783 - 1876.760129419146*eta + 17608.172965575013*eta2 - 40928.07304790013*eta3)*S2 + (-1.28324755577382 - 138.36970336658558*eta + 708.1455154504333*eta2 - 273.23750933544176*eta3)*S3 + (1.8403161863444328 + 2009.7361967331492*eta - 18636.271414571278*eta2 + 42379.205045791656*eta3)*S4+dchi*Seta*eta2*(-105.34550407768225 - 1566.1242344157668*chi1*eta + 1566.1242344157668*chi2*eta + 2155.472229664981*eta*S)

            p1 = 4414.11 + 4.21564/eta - 10687.8*eta + 58234.6*eta2 - 64068.40000000001*eta3 - 704442. .*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6+((6.39833 - 610.267*eta + 2095.72*eta2 - 3970.89*eta3)*S + (22.956700000000005 - 99.1551*eta + 331.593*eta2 - 794.79*eta3)*S2 + (10.4333 + 43.8812*eta - 541.261*eta2 + 294.289*eta3)*S3 + eta*(106.047 - 1569.0299999999997*eta + 4810.61*eta2)*S4)/(eta)+132.244*Seta*eta*(chi1*(6.227738120444028 - 1. .*eta) + chi2*(-6.227738120444028 + 1. .*eta))
            p2 = 3980.7 + 0.956703/eta - 6202.38*eta + 29218.1*eta2 + 24484.2*eta3 - 807629. .*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6+((1.92692 - 226.825*eta + 75.246*eta2 + 1291.56*eta3)*S + (15.328700000000001 - 99.1551*eta + 608.328*eta2 - 2402.94*eta3)*S2 + (10.4333 + 43.8812*eta - 541.261*eta2 + 294.289*eta3)*S3 + eta*(106.047 - 1569.0299999999997*eta + 4810.61*eta2)*S4)/(eta)+132.244*Seta*eta*(chi1*(2.5769789177580837 - 1. .*eta) + chi2*(-2.5769789177580837 + 1. .*eta))
            p3 = 3416.57 + 2308.63*eta - 84042.9*eta2 + 1.01936e6*eta3 - 6.0644e6*eta4 + 1.76399e7*eta5 - 2.0065e7*eta6+(24.6295 - 282.354*eta - 2582.55*eta2 + 12750. .*eta3)*S + (433.675 - 8775.86*eta + 56407.8*eta2 - 114798. .*eta3)*S2 + (559.705 - 10627.4*eta + 61581. .*eta2 - 114029. .*eta3)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4+63.9466*dchi*Seta*eta2
            p4 = 3307.49 - 476.909*eta - 5980.37*eta2 + 127610. .*eta3 - 919108. .*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6+(-5.02553 - 282.354*eta + 1291.56*eta2)*S + (-43.8823 + 740.123*eta - 2402.94*eta2)*S2 + (43.8812 - 370.362*eta + 294.289*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4+-132.244*(dchi)*Seta*eta2
            p5 = 3259.03 - 3967.58*eta + 111203. .*eta2 - 1.81883e6*eta3 + 1.73811e7*eta4 - 9.56988e7*eta5 + 2.75056e8*eta6 - 3.15866e8*eta7+(19.7509 - 1104.53*eta + 3810.18*eta2)*S + (-230.07 + 2314.51*eta - 5944.49*eta2)*S2 + (-201.633 + 2183.43*eta - 6233.99*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4+112.714*dchi*Seta*eta2
            p6 = 3259.03 - 3967.58*eta + 111203. .*eta2 - 1.81883e6*eta3 + 1.73811e7*eta4 - 9.56988e7*eta5 + 2.75056e8*eta6 - 3.15866e8*eta7+(19.7509 - 1104.53*eta + 3810.18*eta2)*S + (-230.07 + 2314.51*eta - 5944.49*eta2)*S2 + (-201.633 + 2183.43*eta - 6233.99*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4+112.714*dchi*Seta*eta2
            
        elseif ell_emm == 44   
            
            #Insp_Phase_44_lambda = 5.254484747463392 - 21.277760168559862*eta + 160.43721442910618*eta2 - 1162.954360723399*eta3 + 1685.5912722190276*eta4 - 1538.6661348106031*eta5+(0.007067861615983771 - 10.945895160727437*eta + 246.8787141453734*eta2 - 810.7773268493444*eta3)*S + (0.17447830920234977 + 4.530539154777984*eta - 176.4987316167203*eta2 + 621.6920322846844*eta3)*S2+-8.384066369867833*dchi*Seta*eta2

            p1 = 4349.66 + 4.34125/eta - 8202.33*eta + 5534.1*eta2 + 536500. .*eta3 - 4.33197e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6+((12.0704 - 528.098*eta + 1822.9100000000003*eta2 - 9349.73*eta3 + 17900.9*eta4)*S + (10.4092 + 253.334*eta - 5452.04*eta2 + 35416.6*eta3 - 71523. .*eta4)*S2 + eta*(492.60300000000007 - 9508.5*eta + 57303.4*eta2 - 109418. .*eta3)*S3)/(eta)+-262.143*Seta*eta*(chi1*(-3.0782778864970646 - 1. .*eta) + chi2*(3.0782778864970646 + 1. .*eta))
            p2 = 3804.19 + 0.66144/eta - 2421.77*eta - 33475.8*eta2 + 665951. .*eta3 - 4.50145e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6+((5.83038 - 172.047*eta + 926.576*eta2 - 7676.87*eta3 + 17900.9*eta4)*S + (6.17601 + 253.334*eta - 5672.02*eta2 + 35722.1*eta3 - 71523. .*eta4)*S2 + eta*(492.60300000000007 - 9508.5*eta + 57303.4*eta2 - 109418. .*eta3)*S3)/(eta)+-262.143*Seta*eta*(chi1*(-1.0543062374352932 - 1. .*eta) + chi2*(1.0543062374352932 + 1. .*eta))
            p3 = 3308.97 + 2353.58*eta - 66340.1*eta2 + 777272. .*eta3 - 4.64438e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6+(-21.5697 + 926.576*eta - 7989.26*eta2 + 17900.9*eta3)*S + (353.539 - 6403.24*eta + 37599.5*eta2 - 71523. .*eta3)*S2 + (492.603 - 9508.5*eta + 57303.4*eta2 - 109418. .*eta3)*S3+262.143*(dchi)*Seta*eta2
            p4 = 3245.63 - 928.56*eta + 8463.89*eta2 - 17422.6*eta3 - 165169. .*eta4 + 908279. .*eta5 - 1.31138e6*eta6+(32.506 - 590.293*eta + 3536.61*eta2 - 6758.52*eta3)*S + (-25.7716 + 738.141*eta - 4867.87*eta2 + 9129.45*eta3)*S2 + (-15.7439 + 620.695*eta - 4679.24*eta2 + 9582.58*eta3)*S3+87.0832*dchi*Seta*eta2
            p5 = 3108.38 + 3722.46*eta - 119588. .*eta2 + 1.92148e6*eta3 - 1.69796e7*eta4 + 8.39194e7*eta5 - 2.17143e8*eta6 + 2.2829700000000003e8*eta7+(118.319 - 529.854*eta)*eta*S + (21.0314 - 240.648*eta + 516.333*eta2)*S2 + (20.3384 - 356.241*eta + 999.417*eta2)*S3+97.1364*(dchi)*Seta*eta2
            p6 = 3096.03 + 986.752*eta - 20371.1*eta2 + 220332. .*eta3 - 1.31523e6*eta4 + 4.29193e6*eta5 - 6.01179e6*eta6+(-9.96292 - 118.526*eta + 2255.76*eta2 - 6758.52*eta3)*S + (-14.4869 + 370.039*eta - 3605.8*eta2 + 9129.45*eta3)*S2 + (17.0209 + 70.1931*eta - 3070.08*eta2 + 9582.58*eta3)*S3+23.0759*dchi*Seta*eta2
        end    

        # end IMRPhenomXHM_FillPhaseFitsArray




        PNdominant = ampNorm * (2. /emm)^( -7/6.) # = Pi * Sqrt(2 eta/3) (2Pi /m)^(-7/6). Miss the f^(-7/6). The pi power included in ampNorm
        IMRPhenomXHMRingdownAmpVersion = 2



        if ell_emm == 21
            if q == 1
                PNdominantlmpower = 2;
                PNdominantlm = sqrt(2.) / 3. * 1.5 * dchi * 0.5 * (2 * pi / emm)^(2/3.);
            
            else
                PNdominantlmpower = 1;
                PNdominantlm = sqrt(2.) / 3. * delta * (2 * pi / emm)^(1/3.);
            end

        elseif ell_emm == 33
            if q == 1
                PNdominantlmpower = 4;
                PNdominantlm = 0.75 * sqrt(5. /7) * dchi * 0.5 * (65/24. - 28/3. * eta) * (2 * pi / emm)^(4/3.);
            
            else
                PNdominantlmpower = 1;
                PNdominantlm = 0.75 * sqrt(5. /7) * delta * (2 * pi / emm)^(1/3.);
            end
        
        elseif ell_emm == 32
            PNdominantlmpower = 2;
            PNdominantlm = 0.75 * sqrt(5. /7) * (2 * pi / emm)^(1/3.);
            IMRPhenomXHMRingdownAmpVersion = 1
        
        elseif ell_emm == 44
            PNdominantlmpower = 2;
            PNdominantlm = 4/9. * sqrt(10/7.) * (1 - 3 * eta) * (2 * pi / emm)^(1/3.);

        end
        PNdominantlm = abs(PNdominantlm) # We only care about the absolute value of it

        InspRescaleFactor = 0;
        RDRescaleFactor = 0;
        InterRescaleFactor = 0;

        if ell_emm ==32
            nCollocPtsInterPhase=6;
        
        else
            nCollocPtsInterPhase=5;
        end
        
        N_MAX_COEFFICIENTS_PHASE_INS = 13       #Maximun number of coefficients of the inspiral ansatz
        N_MAX_COEFFICIENTS_PHASE_INTER = 6      #Maximun number of coefficients of the intermediate ansatz
        N_MAX_COEFFICIENTS_PHASE_RING = 8       #Maximun number of coefficients of the ringdown ansatz. 5 + 3 optional for flattened phi'
        
        N_MAX_COEFFICIENTS_AMPLITUDE_INS = 3    #Maximun number of collocation points in the inspiral
        N_MAX_COEFFICIENTS_AMPLITUDE_INTER = 8  #Maximun number of collocation points in the intermediate. New release 4 coll points + 2x2 boundaries
        N_MAX_COEFFICIENTS_AMPLITUDE_RING = 6   #Maximun number of coefficients in the ringdown. We only use 3 degrees of freedom, but we have the double to store fits of 3 coefficients or 3 collocation points 
        N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX = 4 
            
        #= Transform IMRPhenomXHMIntermediateAmpVersion number to int array defining what to do for each collocation point =#
        #= 0: don't use coll point, 1: use point, 2: use point and derivative (this only applies for boundaries) =#
        # e.g. VersionCollocPtsInter = {1, 1, 1, 1, 0, 2} -> use the two boundaries, add derivative to the right one, skip third collocation point
        if ell_emm == 21
            VersionCollocPtsInter = [1, 1, 0, 1, 0, 2]
            CollocPtsInter = 4 # ANDREA even if some of them not considered 

        else 
            VersionCollocPtsInter = [2, 1, 1, 1, 1, 2]
            CollocPtsInter = 6
        end
        
        IMRPhenomXHMIntermediateAmpFreqsVersion = 0

        nCoefficientsInter = sum(VersionCollocPtsInter) # Number of coefficients in the intermediate ansatz # 5 or 10
        #= The number of coefficients in the intermediate ansatz cannot be larger than the number of available collocation points in IMRPhenomXHMIntermediateAmpVersion 
        If e.g. IMRPhenomXHMIntermediateAmpVersion has 6 slots, the maximum number of coefficients would be 6 + 2 because we count for the derivatives at the boundaries. =#
        


        #Take the cutting frequencies at the inspiral and ringdown
        fAmpMatchIN  = _fcutInsp(model, eta, chi1, emm, fMECOlm)
        fAmpMatchIM  = _fcutRD(model, fring, fdamp, ModeMixingOn, fring_22, fdamp_22)
        fAmpRDfalloff = fring + 2 * fdamp;


        ##### from IMRPhenomXHM_GetPNAmplitudeCoefficients


        #= 
        Take Frequency Domain Post-Newtonian Amplitude Coefficients =#
        #=** Post-Newtonian Inspiral Ansatz Coefficients **=#

        #= Fill pAmp with the coefficients of the power series in frequency of the Forier Domain Post-Newtonian Inspiral Ansatz.
        The 21 mode by default does not used the power series because it breaks down before the end of the inspiral, but it corresponding power series is available here. =#

        #= The ansatz in Fourier Domain is built as follows: we multiply the Time-domain Post-Newtonian series up to 3PN by the phasing factor given by the Stationary-Phase-Approximation,
        and the we reexpand in powers of f up to 3PN.
        The only difference with the 21 is that this does not perform the last reexpansion, so it is not a real power series but it is a quantity better behaved. =#
        #= The coefficients below correspond to those in eqs E10-E14 in arXiv:2001.10914 =#

        #int inspversion = pWFHM->IMRPhenomXHMInspiralAmpFitsVersion;

        #= The coefficients below correspond to those in eqs E10-E14 in arXiv:2001.10914 =#
        # TaylorF2 PN Amplitude Coefficients
        prefactors = [sqrt(2)/3., 0.75*sqrt(5/7.), sqrt(5/7.)/3., 4*sqrt(2)/9*sqrt(5/7.)]; #Global factors of each PN hlm
        PNglobalfactor = (2. /(emm))^(-7/6.)*prefactors[index]; #This is to compensate that we rescale data with the leading order of the 22
 
        
        useFAmpPN = 0
        if ell_emm == 21

            pnInitial = 0.;
            pnOneThird = delta*pi^(1. / 3.)*(2.)^(1. / 3.);
            pnTwoThirds = (-3*(chiA + chiS*delta))/2. *pi^(2. / 3.)*2. ^(2. / 3.);
            pnThreeThirds = (335*delta + 1404*delta*eta)/672. *pi*2.;
            pnFourThirds = (3427*chiA - (1im*672)*delta + 3427*chiS*delta - 8404*chiA*eta - 3860*chiS*delta*eta - 1344*delta*pi - (1im*672)*delta*log(16))/1344. *pi ^(4. / 3.)*2. ^(4. / 3.);
            pnFiveThirds = (-155965824*chiA*chiS - 964357*delta + 432843264*chiA*chiS*eta - 23670792*delta*eta + 24385536*chiA*pi + 24385536*chiS*delta*pi - 77982912*delta*chiA*chiA + 81285120*delta*eta*chiA*chiA - 77982912*delta*chiS*chiS + 39626496*delta*eta*chiS*chiS + 21535920*delta*eta*eta)/8.128512e6*pi ^(5. / 3.)*2. ^(5. / 3.);
            pnSixThirds = (143063173*chiA - (1im*1350720)*delta + 143063173*chiS*delta - 546199608*chiA*eta - (1im*72043776)*delta*eta - 169191096*chiS*delta*eta - 9898560*delta*pi + 20176128*delta*eta*pi - (1im*5402880)*delta*log(2) - (1im*17224704)*delta*eta*log(2) + 61725888*chiS*delta*chiA*chiA - 81285120*chiS*delta*eta*chiA*chiA + 20575296*(chiA)^(3) - 81285120*eta*(chiA)^(3) + 61725888*chiA*chiS*chiS - 165618432*chiA*eta*chiS*chiS + 20575296*delta*(chiS)^(3) - 1016064*delta*eta*chiS*chiS*chiS + 128873808*chiA*eta*eta - 3859632*chiS*delta*eta*eta)/5.419008e6*pi ^2 *2. ^2;


        elseif ell_emm == 33

            pnInitial = 0.;
            pnOneThird = delta*pi^(1. / 3.)*(2. / 3.)^(1. / 3.);
            pnTwoThirds = 0.;
            pnThreeThirds = (-1945*delta + 2268*delta*eta)/672. *pi*(2. / 3.);
            pnFourThirds = (325*chiA - (1im*504)*delta + 325*chiS*delta - 1120*chiA*eta - 80*chiS*delta*eta + 120*delta*pi + (1im*720)*delta*log(1.5))/120. *pi ^(4. / 3.)*(2. / 3.) ^(4. / 3.);
            pnFiveThirds = (-2263282560*chiA*chiS - 1077664867*delta + 9053130240*chiA*chiS*eta - 5926068792*delta*eta - 1131641280*delta*chiA*chiA + 4470681600*delta*eta*chiA*chiA - 1131641280*delta*chiS*chiS + 55883520*delta*eta*chiS*chiS + 2966264784*delta*eta*eta)/4.4706816e8*pi ^(5. / 3.)*(2. / 3.) ^(5. / 3.);
            pnSixThirds = (22007835*chiA + (1im*26467560)*delta + 22007835*chiS*delta - 80190540*chiA*eta - (1im*98774368)*delta*eta - 31722300*chiS*delta*eta - 9193500*delta*pi + 17826480*delta*eta*pi - (1im*37810800)*delta*log(1.5) + (1im*37558080)*delta*eta*log(1.5) - 12428640*chiA*eta*eta - 6078240*chiS*delta*eta*eta)/2.17728e6*pi ^2*(2. / 3.) ^2;
            
            
        elseif ell_emm == 32

            pnInitial = 0.;
            pnOneThird = 0.;
            pnTwoThirds = (-1 + 3*eta)*pi^(2. / 3.);
            pnThreeThirds = -4*chiS*eta*pi;
            pnFourThirds = (10471 - 61625*eta + 82460*eta*eta)/10080. *pi ^(4. / 3.);
            pnFiveThirds = ((1im*2520) - 3955*chiS - 3955*chiA*delta - (1im*11088)*eta + 10810*chiS*eta + 11865*chiA*delta*eta - 12600*chiS*eta*eta)/840. *pi ^(5. / 3.);
            pnSixThirds = (824173699 + 2263282560*chiA*chiS*delta - 26069649*eta - 15209631360*chiA*chiS*delta*eta + 3576545280*chiS*eta*pi + 1131641280*chiA*chiA - 7865605440*eta*chiA*chiA + 1131641280*chiS*chiS - 11870591040*eta*chiS*chiS - 13202119896*eta*eta + 13412044800*chiA*chiA*eta*eta + 5830513920*chiS*chiS*eta*eta + 5907445488*(eta)^(3))/4.4706816e8*pi ^2;
            
            
        elseif ell_emm == 44

            pnInitial = 0.;
            pnOneThird = 0.;
            pnTwoThirds = (1 - 3*eta)*pi^(2. / 3.)*(0.5)^(2. / 3.);
            pnThreeThirds = 0.;
            pnFourThirds = (-158383 + 641105*eta - 446460*eta*eta)/36960. *pi ^(4. / 3.)*(0.5) ^(4. / 3.);
            pnFiveThirds = ((1im*-1008) + 565*chiS + 565*chiA*delta + (1im*3579)*eta - 2075*chiS*eta - 1695*chiA*delta*eta + 240*pi - 720*eta*pi + (1im*960)*log(2) - (1im*2880)*eta*log(2) + 1140*chiS*eta*eta)/120. *pi ^(5. / 3.)*(0.5) ^(5. / 3.);
            pnSixThirds = (7888301437 - 147113366400*chiA*chiS*delta - 745140957231*eta + 441340099200*chiA*chiS*delta*eta - 73556683200*chiA*chiA + 511264353600*eta*chiA*chiA - 73556683200*chiS*chiS + 224302478400*eta*chiS*chiS + 2271682065240*eta*eta - 871782912000*chiA*chiA*eta*eta - 10897286400*chiS*chiS*eta*eta - 805075876080*(eta)^(3))/2.90594304e10*pi ^2*(0.5) ^2;
        end

        ##### end IMRPhenomXHM_GetPNAmplitudeCoefficients
        InspiralAmpVeto=0;
        IntermediateAmpVeto=0;
        RingdownAmpVeto=0;
        IMRPhenomXHMInspiralAmpVersion = 123
        nCollocPtsInspAmp = 3
        nCoefficientsRDAux = 0
        nCollocPtsRDPhase = 0


        timeshift =0.;
        phaseshift=0.;
        phiref22 = 0.; 
        fRDAux = 0;
        if ModeMixingOn == true
            nCollocPtsRDAux = 2;    
            nCoefficientsRDAux = 4;
            fRDAux = fring - fdamp; #v2
            nCollocPtsRDPhase = 4

        end
        
        if ModeMixingOn == true
            #// For the 32 we need the speroidal coefficients for the phase and the amplitude coefficients of the 22.
            #         Reconstruct the coefficients of the ringdown phase in a spheroidal-harmonics basis --only for the 32 mode    
            # Andrea the coeff for modemixing are alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S
            function GetSpheroidalCoefficients(mc, eta, chi1, chi2, nCollocPtsRDPhase, fring, fdamp, fring_22, fdamp_22, S, dchi, Seta, Phase22, ModeMixingCoeffs, fRef)

                nCollocationPts_RD_Phase = nCollocPtsRDPhase;
                CollocValuesPhaseRingdown= Vector{typeofFD}(undef, nCollocationPts_RD_Phase)
                CollocFreqsPhaseRingdown= Vector{typeofFD}(undef, nCollocationPts_RD_Phase)


                #********************* RINGDOWN PHASE COLLOCATION POINTS *****************

                #// this function initializes the frequencies of the collocation points for the spheroidal ringdown reconstruction
                # from IMRPhenomXHM_Ringdown_CollocPtsFreqs

                #//The frequencies are stored in the struct pPhase, which is initialized with new values as the code starts to work on each mode.
                
                fringlm = fring 
                fdamplm = fdamp
                fring22 = fring_22
                CollocationPointsFreqsPhaseRD = Vector{typeofFD}(undef, nCollocationPts_RD_Phase) #4

                # version 122019
                CollocationPointsFreqsPhaseRD[1] = fring22;
                CollocationPointsFreqsPhaseRD[3] = fringlm - 0.5*fdamplm;
                CollocationPointsFreqsPhaseRD[2] = fringlm - 1.5*fdamplm;
                CollocationPointsFreqsPhaseRD[4] = fringlm + 0.5*fdamplm;

                #version 122022 but not yet implemented
                # CollocationPointsFreqsPhaseRD[1] = fring22 - fdamp22;
                # CollocationPointsFreqsPhaseRD[2] = fring22;
                # CollocationPointsFreqsPhaseRD[3] = (fring22 + fringlm) * 0.5;
                # CollocationPointsFreqsPhaseRD[4] = fringlm;
                # CollocationPointsFreqsPhaseRD[5] = fringlm + fdamplm;

                
                #//32 Spheroidal

                RD_Phase_32_p1 = 3169.372056189274 + 426.8372805022653*eta - 12569.748101922158*eta2 + 149846.7281073725*eta3 - 817182.2896823225*eta4 + 1.5674053633767858e6*eta5+(19.23408352151287 - 1762.6573670619173*eta + 7855.316419853637*eta2 - 3785.49764771212*eta3)*S + (-42.88446003698396 + 336.8340966473415*eta - 5615.908682338113*eta2 + 20497.5021807654*eta3)*S2 + (13.918237996338371 + 10145.53174542332*eta - 91664.12621864353*eta2 + 201204.5096556517*eta3)*S3 + (-24.72321125342808 - 4901.068176970293*eta + 53893.9479532688*eta2 - 139322.02687945773*eta3)*S4 + (-61.01931672442576 - 16556.65370439302*eta + 162941.8009556697*eta2 - 384336.57477596396*eta3)*S5+dchi*Seta*eta2*(641.2473192044652 - 1600.240100295189*chi1*eta + 1600.240100295189*chi2*eta + 13275.623692212472*eta*S)
                RD_Phase_32_p2 = 3131.0260952676376 + 206.09687819102305*eta - 2636.4344627081873*eta2 + 7475.062269742079*eta3+(49.90874152040307 - 691.9815135740145*eta - 434.60154548208334*eta2 + 10514.68111669422*eta3)*S + (97.3078084654917 - 3458.2579971189534*eta + 26748.805404989867*eta2 - 56142.13736008524*eta3)*S2 + (-132.49105074500454 + 429.0787542102207*eta + 7269.262546204149*eta2 - 27654.067482558712*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4+477.51566939885424*dchi*Seta*eta2
                RD_Phase_32_p3 = 3082.803556599222 + 76.94679795837645*eta - 586.2469821978381*eta2 + 977.6115755788503*eta3+(45.08944710349874 - 807.7353772747749*eta + 1775.4343704616288*eta2 + 2472.6476419567534*eta3)*S + (95.57355060136699 - 2224.9613131172046*eta + 13821.251641893134*eta2 - 25583.314298758105*eta3)*S2 + (-144.96370424517866 + 2268.4693587493093*eta - 10971.864789147161*eta2 + 16259.911572457446*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4+378.2359918274837*dchi*Seta*eta2
                RD_Phase_32_p4 = 3077.0657367004565 + 64.99844502520415*eta - 357.38692756785395*eta2+(34.793450080444714 - 986.7751755509875*eta - 9490.641676924794*(eta)^(3) + 5700.682624203565*eta2)*S + (57.38106384558743 - 1644.6690499868596*eta - 19906.416384606226*(eta)^(3) + 11008.881935880598*eta2)*S2 + (-126.02362949830213 + 3169.3397351803583*eta + 62863.79877094988*(eta)^(3) - 26766.730897942085*eta2)*S3 + (-169.30909412804587 + 4900.706039920717*eta + 95314.99988114933*(eta)^(3) - 41414.05689348732*eta2)*S4+390.5443469721231*dchi*Seta*eta2
                RD_Phase_32_p5 = 0.

                b = [RD_Phase_32_p1, RD_Phase_32_p2, RD_Phase_32_p3, RD_Phase_32_p4] # 4 coefficients not 5
                A = Matrix{typeofFD}(undef, nCollocationPts_RD_Phase, nCollocationPts_RD_Phase)
                for i in 1:nCollocationPts_RD_Phase

                    CollocFreqsPhaseRingdown[i]  = CollocationPointsFreqsPhaseRD[i];
                    ff = CollocFreqsPhaseRingdown[i]
                    ffm1 = 1. /ff
                    ffm2 = ffm1^2
                    # lorentzian = fdamp / (fdamp * fdamp + (ff - fring) * (ff - fring));
                    # fpowers = [1., ffm1, ffm2, ffm2 * ffm2, lorentzian]
                    fpowers= [1., (fdamp)/((fdamp^2)+(ff-(fring))^(2)), ffm2, ffm2*ffm2]

                    for j in 1:nCollocationPts_RD_Phase
                        A[i,j] = fpowers[j] # Ansatz spheroidal
                    end
                end


                #// solve the linear system of Eq. (6.8)
                x = lu(A) \ b

                # ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)
                alpha0_S = x[1]
                alphaL_S = x[2]
                alpha2_S = x[3]
                alpha4_S = x[4]

                if debug == true
                    println("Spheroidal coefficients")
                    println("alpha0_S: ", alpha0_S)
                    println("alphaL_S: ", alphaL_S)
                    println("alpha2_S: ", alpha2_S)
                    println("alpha4_S: ", alpha4_S)
                end
                #//  we have reconstructed the "shape" of the spheroidal ringdown phase derivative, dphiS, but we need to adjust the relative time and phase shift of the final phiS wrt IMRPhenomX

                #************** time-shift of spheroidal ansatz ******************
                frefRD= fring_22 + fdamp_22
                #// here we call a fit for dphiS(fref)-dphi22(fref)
                #*************** 32 specific fits **************

                noSpin = 11.851438981981772 + 167.95086712701223*eta - 4565.033758777737*eta2 + 61559.132976189896*eta3 - 364129.24735853914*eta4 + 739270.8814129328*eta5;
                eqSpin = (9.506768471271634 + 434.31707030999445*eta - 8046.364492927503*eta2 + 26929.677144312944*eta3)*S + (-5.949655484033632 - 307.67253970367034*eta + 1334.1062451631644*eta2 + 3575.347142399199*eta3)*S2 + (3.4881615575084797 - 2244.4613237912527*eta + 24145.932943269272*eta2 - 60929.87465551446*eta3)*S3 + (15.585154698977842 - 2292.778112523392*eta + 24793.809334683185*eta2 - 65993.84497923202*eta3)*S4;
                uneqSpin = 465.7904934097202*dchi*sqrt(1. -4. *eta)*eta2;
                tshift_ = noSpin + eqSpin + uneqSpin;
                
                #// we compute dphi22(fref)
                #Phase22 = Phase_22_ConnectionCoefficients(mc, eta, chi1, chi2) # ANDREA maybe no need to redo this function after 22, use the same coefficients
                timeshift=TimeShift_22(model, eta, Seta, S, dchi, fring_22, fdamp_22, Phase22)
                
                phi0_S = 0;

                #// we impose that dphiS(fref)-dphi22(fref) has the value given by our fit
                dphi22ref = 1. /eta*_completePhaseDer(model, frefRD, Phase22, fdamp_22, fring_22)+timeshift; # ANDREA it is subtracted from timeshift, maybe do it properly linb - 2. *pi*(500. +psi4tostrain)
                ModeMixingCoeffs = [alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S]
                alpha0_S = alpha0_S +dphi22ref+tshift_-RD_Phase_Ansatz(model, frefRD, 0., 0., 0., fring, fdamp, ModeMixingCoeffs,  true)
                
                ModeMixingCoeffs = [alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S]

                #************** phase-shift of spheroidal ansatz ******************
                frefRD = fring_22 

                # we compute phi22(fref) 
                #phi0 = 1.
                
                #MfRef = fRef * M * GMsun_over_c3 

                #tmp = Phi(PhenomXAS(), [fRef, frefRD/(M * GMsun_over_c3 )] , mc, eta, chi1, chi2) # Andrea da rifare con coeffs in Phase22
                tmp =_completePhase(model, [MfRef, frefRD], Phase22, fdamp_22, fring_22)
                phiref22 = -1. ./ eta*tmp[1] - timeshift*MfRef - phaseshift + 2.0*phi0 + pi/4;
                phi22ref= 1. /eta*tmp[2] + timeshift*frefRD + phaseshift + phiref22;

                
                #// we call a fit for Mod[phiS(fref)-phi22(fref),2*Pi]

                ##### from IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift

                noSpin = -1.3328895897490733 - 22.209549522908667*eta + 1056.2426481245027*eta2 - 21256.376324666326*eta3 + 246313.12887984765*eta4 - 1.6312968467540336e6*eta5 + 5.614617173188322e6*eta6 - 7.612233821752137e6*eta7;
                eqSpin = (S*(-1.622727240110213 + 0.9960210841611344*S - 1.1239505323267036*S2 - 1.9586085340429995*S3 + eta2*(196.7055281997748 + 135.25216875394943*S + 1086.7504825459278*S2 + 546.6246807461155*S3 - 312.1010566468068*S4) + 0.7638287749489343*S4 + eta*(-47.475568056234245 - 35.074072557604445*S - 97.16014978329918*S2 - 34.498125910065156*S3 + 24.02858084544326*S4) + eta3*(62.632493533037625 - 22.59781899512552*S - 2683.947280170815*S2 - 1493.177074873678*S3 + 805.0266029288334*S4)))/(-2.950271397057221 + 1. .*S);
                uneqSpin = (Seta*(chi2*(eta)^(2.5)*(88.56162028006072 - 30.01812659282717*S) + chi2*eta2*(43.126266433486435 - 14.617728550838805*S) + chi1*eta2*(-43.126266433486435 + 14.617728550838805*S) + chi1*(eta)^(2.5)*(-88.56162028006072 + 30.01812659282717*S)))/(-2.950271397057221 + 1. .*S);
                total = noSpin + eqSpin + uneqSpin;


                phishift= total;

                #//we adjust the relative phase of our reconstruction
                phi0_S= phi22ref - RD_Phase_AnsatzInt(model, frefRD, 0., 0., 0., fring, fdamp, ModeMixingCoeffs,  true) + phishift;

                return alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S
            end

            alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S = GetSpheroidalCoefficients(mc, eta, chi1, chi2, nCollocPtsRDPhase, fring, fdamp, fring_22, fdamp_22, S, dchi, Seta, Phase22, ModeMixingCoeffs, fRef) # da fare 27/8, make ampl(xas) che dia in output tutti i coeff needed
            ModeMixingCoeffs = [alpha0_S, alphaL_S, alpha2_S, alpha4_S, phi0_S]

            if debug == true
                println()
                println("After GetSpheroidalCoefficients")
                println("alpha0_S: ", alpha0_S)
                println("alphaL_S: ", alphaL_S)
                println("alpha2_S: ", alpha2_S)
                println("alpha4_S: ", alpha4_S)
                println("phi0_S: ", phi0_S)
            end
            
        end

        ##### from IMRPhenomXHM_Get_Inspiral_Amp_Coefficients(pAmp, pWFHM, pWF22);

        CollocationPointsFreqsAmplitudeInsp = [0.5*fAmpMatchIN, 0.75*fAmpMatchIN, fAmpMatchIN]

        fcutInsp_seven_thirds = fAmpMatchIN^(7/3)
        fcutInsp_eight_thirds = fAmpMatchIN^(8/3)
        fcutInsp_three = fAmpMatchIN^3

        function Inspiral_PNAmp_Ansatz(f)   # input are CollocationPointsFreqsAmplitudeInsp
            #Return the Fourier Domain Post-Newtonian ansatz up to 3PN without the pseudoPN terms for a particular frequency
            #### from IMRPhenomXHM_Inspiral_PNAmp_Ansatz
            #// The 21 mode is special, is not a power series


            #This returns the amplitude strain rescaled with the prefactor of the 22 mode: divided by sqrt(2*eta/3.)/pi^(1/6)*Mf^(-7/6)
            
            pnAmp = abs(pnInitial + f^(1/3) * pnOneThird + f^(2/3)   * pnTwoThirds  + f * pnThreeThirds + f^(4/3)  * pnFourThirds + f^(5/3)  * pnFiveThirds + f^(2) * pnSixThirds);
            
            pnAmp *= PNglobalfactor * f^(-7/6) * ampNorm;
 

            return pnAmp;

        end

        PNAmplitudeInsp1 = Inspiral_PNAmp_Ansatz(CollocationPointsFreqsAmplitudeInsp[1])
        PNAmplitudeInsp2 = Inspiral_PNAmp_Ansatz(CollocationPointsFreqsAmplitudeInsp[2])
        PNAmplitudeInsp3 = Inspiral_PNAmp_Ansatz(CollocationPointsFreqsAmplitudeInsp[3])
        PNAmplitudeInsp = [PNAmplitudeInsp1, PNAmplitudeInsp2, PNAmplitudeInsp3]

        #### from IMRPhenomXHM_Get_Inspiral_Amp_Coefficients

        f1 = CollocationPointsFreqsAmplitudeInsp[1]
        f2 = CollocationPointsFreqsAmplitudeInsp[2]
        f3 = CollocationPointsFreqsAmplitudeInsp[3]
        finsp = CollocationPointsFreqsAmplitudeInsp[3]

        v1 = (Insp_v1 - PNAmplitudeInsp[1]) / PNdominant / f1^(-7/6)
        v2 = (Insp_v2 - PNAmplitudeInsp[2]) / PNdominant / f2^(-7/6)
        v3 = (Insp_v3 - PNAmplitudeInsp[3]) / PNdominant / f3^(-7/6)


        #InspiralCoefficient = zeros(3) #CollocationPointsValuesAmplitudeInsp
        InspiralCoefficient1 =  (finsp^(7/3)*(-(f1^3 * f3^(8/3) * v2) + f1^(8/3) * f3^3 * v2 + f2^3 * (f3^(8/3) * v1 - f1^(8/3) * v3) + f2^(8/3) * (-(f3^3 * v1) + f1^3 * v3))) / (f1^(7/3) * (f1^(1/3) - f2^(1/3)) * f2^(7/3) * (f1^(1/3) - f3^(1/3)) * (f2^(1/3) - f3^(1/3)) * f3^(7/3))
        InspiralCoefficient2 =  (finsp^(8/3)*(f1^3 * f3^(7/3) * v2 - f1^(7/3) * f3^3 * v2 + f2^3 * (-f3^(7/3) * v1 + f1^(7/3) * v3) + f2^(7/3) * (f3^3 * v1 - f1^3 * v3))) / (f1^(7/3) * (f1^(1/3) - f2^(1/3)) * f2^(7/3) * (f1^(1/3) - f3^(1/3)) * (f2^(1/3) - f3^(1/3)) * f3^(7/3))
        InspiralCoefficient3 = (finsp^3*(f1^(7/3)*(-f1^(1/3) + f3^(1/3))*f3^(7/3)*v2 + f2^(7/3)*(-(f3^(8/3)*v1) + f1^(8/3)*v3) + f2^(8/3)*(f3^(7/3)*v1 - f1^(7/3)*v3)))/(f1^(7/3)*(f1^(1/3) - f2^(1/3))*f2^(7/3)*(f1^(1/3) - f3^(1/3))*(f2^(1/3) - f3^(1/3))*f3^(7/3));
        InspiralCoefficient = [InspiralCoefficient1, InspiralCoefficient2, InspiralCoefficient3]


        ##### from IMRPhenomXHM_RD_Amp_Coefficients

        CollocationPointsFreqsAmplitudeRD_1 = fring_array[index] - fdamp_array[index]
        CollocationPointsFreqsAmplitudeRD_2 = fring_array[index] 
        CollocationPointsFreqsAmplitudeRD_3 = fring_array[index] + fdamp_array[index]

        if  rdcp3 >= rdcp2 * rdcp2 / rdcp1 
            rdcp3 = 0.5 * rdcp2 * rdcp2 / rdcp1;
        end

        if  rdcp3 > rdcp2 
            rdcp3 = 0.5 * rdcp2;
        end

        if rdcp1 < rdcp2 && rdcp3 > rdcp1
            rdcp3 = rdcp1;
        end


        RDCoefficient = Vector{typeofFD}(undef, 5)
        if ell_emm == 32    # IMRPhenomXHMRingdownAmpVersion = 1
            RDCoefficient[1] = abs(alambda)
            RDCoefficient[2] = lambda
            RDCoefficient[3] = sigma
        else    # IMRPhenomXHMRingdownAmpVersion = 2
            deno = (sqrt(rdcp1 / rdcp3) - (rdcp1 / rdcp2));
            if deno <= 0
                deno = 1e-16;
            end
            RDCoefficient[1] = rdcp1 * fdamp / deno;
            RDCoefficient[3] = sqrt(RDCoefficient[1] / (rdcp2 * fdamp));
            RDCoefficient[2] = 0.5 * RDCoefficient[3] * log(rdcp1 / rdcp3);
        end

        function RD_Amp_Ansatz(Mf)

            #************* Amplitude Ringdown Ansatz ************

            # For the modes with mixing this is the ansatz of the spheroidal part.
            ##### from IMRPhenomXHM_RD_Amp_Ansatz
            #fAmpRDfalloff = 0.0;    
            ff = Mf
            frd = fring
            fda = fdamp
            dfr = ff - frd;
            ampRD = 0. #zeros(typeof(ff))

            if nCoefficientsRDAux > 0 && ff .< fRDAux   #  fRDAux exists only for mode 32 (mode mixing)
                    # Polynomial 
                fpower = 1.
                for i in 1:nCoefficientsRDAux
                    ampRD += fpower * RDAuxCoefficient[i];
                    fpower  *= ff;
                end

            elseif  fAmpRDfalloff > 0 && ff .> fAmpRDfalloff # high freq part
                ampRD =  RDCoefficient[4] * exp(- RDCoefficient[5] * (ff - fAmpRDfalloff));

            else # Lorentzian with exponential falloff 
                dfd = fda * RDCoefficient[3];
                ampRD =  RDCoefficient[1] * fda / ( exp(RDCoefficient[2] / dfd * dfr) * (dfr * dfr + dfd * dfd)); # * pampNorm * factor;

            end

            return ampRD

        end

        function RD_Amp_DAnsatz(Mf)

        
            #** Derivative of the RD Ansatz for modes without mixing **

            #static IMRPhenomXHM_RD_Amp_DAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

            ff = Mf
            frd = fring
            fda = fdamp
            dfr = ff - frd;
            numerator =  RDCoefficient[1] * (dfr * dfr * RDCoefficient[2] + 2 * fda * dfr * RDCoefficient[3] + fda * fda * RDCoefficient[2] * RDCoefficient[3]^2);
            denom =  (dfr * dfr + fda * fda * RDCoefficient[3]^2);
            denom =  RDCoefficient[3] * denom * denom * exp(dfr * RDCoefficient[2] / (fda * RDCoefficient[3]));
            DampRD =  - numerator / denom;

            return DampRD

        end
    

        tmp = fAmpRDfalloff;
        fAmpRDfalloff = 0;
        RDCoefficient[4] = RD_Amp_Ansatz(tmp)   # this takes as input the value of fAmpRDfalloff, but we use tmp bc this function will be reused for a different purpose (calculate the amplitude and not the coeff of RD)
        # need to make this an array since is an array of freq
        RDCoefficient[5] =  -1.0 * RD_Amp_DAnsatz(tmp) / RDCoefficient[4];
        fAmpRDfalloff = tmp;

        if nCoefficientsRDAux > 0 && ModeMixingOn == true  # only mode 32 has mixing, so it enters here

            # from IMRPhenomXHM_RDAux_Amp_Coefficients
            CollocationPointsFreqsAmplitudeRDAux = Vector{typeofFD}(undef, 4)
            CollocationPointsValuesAmplitudeRDAux = Vector{typeofFD}(undef, 4)

            CollocationPointsValuesAmplitudeRDAux[1] = abs(rdaux1)
            CollocationPointsValuesAmplitudeRDAux[2] = abs(rdaux2)
            

            ampRD = RD_Amp_Ansatz(fRDAux)
            DampRD = RD_Amp_DAnsatz(fRDAux)
            

            CollocationPointsValuesAmplitudeRDAux[3] = ampRD 
            CollocationPointsValuesAmplitudeRDAux[4] = DampRD 

            
            CollocationPointsFreqsAmplitudeRDAux[1] = fAmpMatchIM;
            CollocationPointsFreqsAmplitudeRDAux[2] = 0.5 * (fAmpMatchIM + fRDAux); #// First Chebyshev node
            CollocationPointsFreqsAmplitudeRDAux[3] = fRDAux;
            CollocationPointsFreqsAmplitudeRDAux[4] = fRDAux;
        
        

            # Define linear system of equations 
        
            #// b is the vector with the values of collocation points
            b = CollocationPointsValuesAmplitudeRDAux
            A = Matrix{typeofFD}(undef, nCoefficientsRDAux, nCoefficientsRDAux)
            for i in 1:nCoefficientsRDAux

                #    A = (1, f1, f1^2, f1^3, f1^4)
                #          (1, f2, f2^2, f2^3, f2^4)
                #          (1, f3, f3^2, f3^3, f3^4)
                #          (0,  1,   2f3, 3f3^2, 4f3^3) #andrea: fixed the last row
                #          Until number of collocation points
                #   
                fcollpoint = CollocationPointsFreqsAmplitudeRDAux[i];
                fpower = 1.; #// 1, f, f^2, f^3, f^4, ...
                if i < nCoefficientsRDAux 
                    for j in 1:nCoefficientsRDAux
                        A[i,j] = fpower
                        fpower *= fcollpoint;
                    end
                else # Last row of the matrix for the derivative
                    fpower = 1.;
                    A[i,1] = 0.;
                    for j in 2: nCoefficientsRDAux
                        A[i,j] = (j-1) * fpower;
                        fpower *= fcollpoint;
                    end
                end
            end
    
        
            # We now solve the system A x = b via an LU decomposition. x is the solution vector 

            RDAuxCoefficient = lu(A) \ b

        end


        ##### end IMRPhenomXHM_RD_Amp_Coefficients

                  
        ##### from IMRPhenomXHM_Intermediate_Amp_Coefficients

        # Previously we checked that nCollocPtsInterAmp read from the IMRPhenomXHMIntermediateAmpVersion is equal to the number of free coefficients in the ansatz. 
        nCollocPtsInterAmp = nCoefficientsInter; # ANDREA, nCollocPtsInterAmp is different from pWFHM->nCollocPtsInterAmp, no idea why
        nCollocPtsInterAmp_pWHM = length(VersionCollocPtsInter)

        # Define set of collocation points 
        # from IMRPhenomXHM_Intermediate_Amp_CollocationPoints

        CollocationPointsFreqsAmplitudeInter = Vector{typeofFD}(undef, nCollocPtsInterAmp)
        # Define collocation points frequencies 
        if IMRPhenomXHMIntermediateAmpFreqsVersion == 0
            #// Equispaced. Get boundaries too
            deltaf = (fAmpMatchIM - fAmpMatchIN) / (nCollocPtsInterAmp_pWHM - 1);
            idx = 1;
            for i in 1:nCollocPtsInterAmp_pWHM
                if VersionCollocPtsInter[i] == 1
                    #// Add point
                    CollocationPointsFreqsAmplitudeInter[idx] = fAmpMatchIN + deltaf * (i-1);
                
                elseif VersionCollocPtsInter[i] == 2
                    #// Add point + derivative
                    CollocationPointsFreqsAmplitudeInter[idx] = fAmpMatchIN + deltaf * (i-1);
                    CollocationPointsFreqsAmplitudeInter[idx + 1] = CollocationPointsFreqsAmplitudeInter[idx];
                end
                idx += VersionCollocPtsInter[i];
            end
        elseif IMRPhenomXHMIntermediateAmpFreqsVersion == 1
            # // Chebyshev. Get boundaries too
            semisum = 0.5 * (fAmpMatchIN + fAmpMatchIM);
            semidif = 0.5 * (fAmpMatchIM - fAmpMatchIN);
            for i in (nCollocPtsInterAmp_pWHM + 2):-1:1
                CollocationPointsFreqsAmplitudeInter[i] = semisum + semidif * cos((i-1) * π / nCollocPtsInterAmp_pWHM);
            end
        end

        # Define values 

        #// Be careful with TGR thing. Inspiral affecting intermediate region.
        tmp_factor = InspRescaleFactor; # 0
        InspRescaleFactor = InterRescaleFactor; # 0
        tmpnCollocPts = 1; # julia index


        ##### from IMRPhenomXHM_Inspiral_Amp_Ansatz
        IMRPhenomXHMInspiralAmpFreqsVersion     = 122022

        function Inspiral_Amp_Ansatz(Mf_)

            #InspAmp is the amplitude strain rescaled with the prefactor of the 22 mode: divided by [sqrt(2*eta/3.)/pi^(1/6) * f^(-7/6)]
    
            # New release. 
            InspAmp = Inspiral_PNAmp_Ansatz(Mf_)
            # println("\nMf_ = ", Mf_)
            # println("InspAmp: ", InspAmp)
            pseudoterms = (Mf_^(7/3) / fcutInsp_seven_thirds * InspiralCoefficient[1]) +
                          + (Mf_^(8/3) / fcutInsp_eight_thirds * InspiralCoefficient[2]) +
                          + (Mf_^3     / fcutInsp_three * InspiralCoefficient[3]);
            pseudoterms *= Mf_^(-7/6) * PNdominant;

            InspAmp += pseudoterms;
                    
        
            return InspAmp
        end

        # Automatic differentiation
        function Inspiral_Amp_NDAnsatz_FD(Mff)
            return ForwardDiff.derivative(Inspiral_Amp_Ansatz, Mff)
        end

        function Inspiral_Amp_NDAnsatz(Mff)
            # Numerical derivative of the inspiral (4th order finite differences)
            #It is used for reconstructing the intermediate region 

            df = 10e-10;
            centralfreq = Mff;

            fun2R = Inspiral_Amp_Ansatz(centralfreq + 2*df);
            funR  = Inspiral_Amp_Ansatz(centralfreq + df);
            funL  = Inspiral_Amp_Ansatz(centralfreq - df);
            fun2L = Inspiral_Amp_Ansatz(centralfreq - 2*df);
            Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

            return Nder;
        

        end

        CollocationPointsValuesAmplitudeInter = Vector{typeofFD}(undef, nCollocPtsInterAmp)
        if VersionCollocPtsInter[1] == 1
            CollocationPointsValuesAmplitudeInter[1] = Inspiral_Amp_Ansatz(fAmpMatchIN)
            tmpnCollocPts += 1;
        
        elseif VersionCollocPtsInter[1] == 2
            CollocationPointsValuesAmplitudeInter[1] = Inspiral_Amp_Ansatz(fAmpMatchIN)
            if debug == true
                println(" No FD: ", Inspiral_Amp_NDAnsatz(fAmpMatchIN))
                println(" FD: ", Inspiral_Amp_NDAnsatz_FD(fAmpMatchIN))
            end

            CollocationPointsValuesAmplitudeInter[2] = Inspiral_Amp_NDAnsatz_FD(fAmpMatchIN);#//CollocationPointsValuesAmplitudeInsp[2];
            tmpnCollocPts += 2;
        end
        InspRescaleFactor = tmp_factor;


        #  Call parameter space fits 
        #  pWFHM->nCollocPtsInterAmp also includes the boundaries, for the parameter space fits we just need the in between points 
        IntermediateAmpFits = [int1, int2, int3, int4]
        for  i in 2:(nCollocPtsInterAmp_pWHM -1)  # 

            if VersionCollocPtsInter[i] == 1
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = abs(IntermediateAmpFits[i-1]); 
                tmpnCollocPts+=1;
            end
        end
        
        tmp_factor = RDRescaleFactor; # 0
        RDRescaleFactor = InterRescaleFactor; # 0 
        

        function RD_Amp_NDAnsatz_FD(Mf)
            return ForwardDiff.derivative(x -> abs(SpheroidalToSpherical(x)), Mf)
        end

        function RD_Amp_NDAnsatz(Mf)

            # ** Derivative of the RD Ansatz for modes with mixing **
            # // It can not be obtained analytically, so we use finite difference of 4th order

                ff = Mf
                df = 10e-10;
                fun2R = ff + 2*df;
                funR  = ff + df;
                funL  = ff - df;
                fun2L = ff - 2*df;

                fun2R = abs(SpheroidalToSpherical(fun2R));

                funR  = abs(SpheroidalToSpherical(funR));  

                funL  = abs(SpheroidalToSpherical(funL));

                fun2L = abs(SpheroidalToSpherical(fun2L));

                Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

            return Nder;
        end


        function SpheroidalToSpherical(Mf)



        # ************************************
        #                                     
        #   Spheroidal -> Spherical rotation  
        #                                     
        # ************************************
        #  The rotation consists of a linear transformation using the mixing coefficients given by Berti [10.1103/PhysRevD.90.064012].

        # h32_spherical = a1 * h22_spheroidal + a2 * h32_spheroidal,  where a1, a2 are the mixing coefficients.

        # Since the 22 is the most dominant mode, it will not show a significant mixing with any other mode,
        # so we can assume that h22_spheroidal = h22_spherical
        # 

        #// In principle this could be generalized to the 43 mode: for the time being, assume the mode solved for is only the 32.
        
        #// Compute the 22 mode using PhenomX functions. This gives the 22 mode rescaled with the leading order.  This is because the 32 is also rescaled.

            ##### from XLALSimIMRPhenomXRingdownAmplitude22AnsatzAnalytical
            # *
            #     "Analytical" phenomenological ringdown ansatz for amplitude. This is used by the higher mode functions but
            #     can also be used to prototype or test model. Convenient wrapper exposed via XLAL.
            # 
            ff = Mf
            gammaD13 = fdamp_22 * Amp22.gamma1 * Amp22.gamma3;
            gammaR   = Amp22.gamma2 / (Amp22.gamma3 * fdamp_22);
            gammaD2  = (fdamp_22 * Amp22.gamma3) * (fdamp_22 * Amp22.gamma3);
            dfr      = ff - fring_22;

            
            amp22 =  exp(- dfr * gammaR ) * (gammaD13) / (dfr*dfr + gammaD2);
            #phiXAS = Phi(PhenomXAS(), ff/(M * GMsun_over_c3 ), mc, eta, chi1, chi2) # the phase of the 22 mode directly from PhenomXAS # from IMRPhenomX_Phase_22
            phiXAS = _completePhase(model, ff, Phase22, fdamp_22, fring_22)
            phi22=1. /eta * phiXAS + timeshift*Mf + phaseshift + phiref22;
            wf22R = amp22 * exp(1im * phi22);
            if IMRPhenomXHMRingdownAmpVersion != 0
                wf22R *= ampNorm * Mf^(-7/6);
            end
            #// Compute 32 mode in spheroidal.
            amplm=RD_Amp_Ansatz(Mf)#, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn);
            philm=RD_Phase_AnsatzInt(model, Mf, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn); # from Phase, if mode mixing is on, we need to compute the phase of the 32 mode
            #// Do the rotation.
            sphericalWF_32=conj(mixingCoeffs[3]) * wf22R + conj(mixingCoeffs[4])*amplm*exp(1im*philm);
            return sphericalWF_32

        
        end


        
        fRD = fAmpMatchIM
        if VersionCollocPtsInter[end] == 1

            if ModeMixingOn == false
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = RD_Amp_Ansatz(fRD)#, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn);#//CollocationPointsValuesAmplitudeRD[0];
            
            else
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = abs(SpheroidalToSpherical(fRD));
            end
            tmpnCollocPts +=1

        else # VersionCollocPtsInter[end] == 2 ANDREA we use only this one
            # Add point + derivative
            if ModeMixingOn == false
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = RD_Amp_Ansatz(fRD);#//CollocationPointsValuesAmplitudeRD[0];
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts + 1] = RD_Amp_DAnsatz(fRD);#//CollocationPointsValuesAmplitudeRD[0];
            
            else

                CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = abs(SpheroidalToSpherical(fRD));
                CollocationPointsValuesAmplitudeInter[tmpnCollocPts + 1] = RD_Amp_NDAnsatz_FD(fRD);
            end
            tmpnCollocPts += 2;
        end

        if debug == true
            println("CollocationPointsValuesAmplitudeInter: ", CollocationPointsValuesAmplitudeInter)
        end

        RDRescaleFactor = tmp_factor;

        # tmpnCollocPts must be the same tahn pWFHM->nCollocPtsInterAmp + 2 (=number of free coefficients in intermediate ansatz) 
        if tmpnCollocPts != (nCoefficientsInter+1)
            error("IMRPhenomXHM_Intermediate_Amp_Coefficients failed. Inconsistent number of collocation points and free parameters.");
        end


        b = Vector{typeofFD}(undef, nCollocPtsInterAmp)
        x = Vector{typeofFD}(undef, nCollocPtsInterAmp)
        A = Matrix{typeofFD}(undef, nCollocPtsInterAmp, nCollocPtsInterAmp)


        #  Define linear system of equations: A x = b 
        #  x is the solution vector: the coefficients of the intermediate ansatz 
        #  b is the vector of collocation points for a set of frequencies 
        #  A is the matrix of multiplicative factors to each coefficient of the ansatz.
        # Each row gives the ansatz evaluated at a collocation point frequency. 
        tmpnCollocPts = 1;
        for i in 1:nCollocPtsInterAmp_pWHM
            # Skip the 0 cases, means that collocation point is not used 
            if VersionCollocPtsInter[i] > 0
                #// Set b vector
                b[tmpnCollocPts] = CollocationPointsValuesAmplitudeInter[tmpnCollocPts];
                #gsl_vector_set(b, tmpnCollocPts, CollocationPointsValuesAmplitudeInter[tmpnCollocPts]);

                #// Set system matrix: Polynomial/f^7/6 at the collocation points frequencies.
                #  A = (1, f1, f1^2, f1^3, f1^4, ...) * f1^(-7/6)
                #     (1, f2, f2^2, f2^3, f2^4, ...) * f2^(-7/6)
                #     ....
                #     Until number of collocation points
                # If one of the conditions is a derivative we substitute by
                #     (0, 1, 2*fi, 3*fi^2, 4*fi^3, ...)
                # 

                fcollpoint = CollocationPointsFreqsAmplitudeInter[tmpnCollocPts];
                fcollpoint_m_seven_sixths = fcollpoint^(-7/6)
                fpower = 1.  #// 1, f, f^2, f^3, f^4, ...

                #// Add equation for Point. Assuming I will always use point, not the derivative alone
                for j in 1:nCollocPtsInterAmp

                    A[tmpnCollocPts,j] = fpower * fcollpoint_m_seven_sixths;
                    fpower *= fcollpoint;
                end
                tmpnCollocPts+=1
                #// Add equation for Derivative
                if VersionCollocPtsInter[i] == 2 
                #gsl_vector_set(b, tmpnCollocPts, CollocationPointsValuesAmplitudeInter[tmpnCollocPts]);
                    b[tmpnCollocPts] = CollocationPointsValuesAmplitudeInter[tmpnCollocPts];
                    fpower = 1. /fcollpoint;
                    for j in 1:nCollocPtsInterAmp
                        derivative = ((j-1.) - 7/6) * fpower * fcollpoint_m_seven_sixths;
                        A[tmpnCollocPts,j] = derivative;
                        fpower *= fcollpoint;
                    end
                    tmpnCollocPts+=1
                end
                
            end #  End non-zero if statement 
        end # End of loop over number of free coefficients. System of equations setup. 

        # tmpnCollocPts must be the same than the number of free coefficients and to nCoefficientsInter 

        if (tmpnCollocPts-1) != nCoefficientsInter 
            error("IMRPhenomXHM_Intermediate_Amp_Coefficients failed. Inconsistent number of collocation points (%i) and free parameters (%i).");
        end

        # We now solve the system A x = b via an LU decomposition. x is the solution vector 
        x = lu(A) \ b

        InterCoefficient = Vector{typeofFD}(undef, nCoefficientsInter)
        # The solution corresponds to the coefficients of the ansatz 
        for i in 1:nCoefficientsInter
            InterCoefficient[i] = x[i];
        end
        if debug == true
            println("b: ", b)
            println("CollocationPointsFreqsAmplitudeInter: ", CollocationPointsFreqsAmplitudeInter)
            println("InterCoefficient: ", InterCoefficient)
        end

        ##### end of IMRPhenomXHM_Intermediate_Amp_Coefficients
        ##### end of IMRPhenomXHM_GetAmplitudeCoefficients

        ##### from IMRPhenomXHM_GetPhaseCoefficients
    
        #  compute phase coefficients 
        #function GetPhaseCoefficients()
        
        
    
        # //intermediate
        # // will be determined by solving linear system at collocation points (see Eqs. (5.5)-(5.6) in the paper)
        c0 = 0.0;
        c1 = 0.0;
        c2 = 0.0;
        c3 = 0.0;
        c4 = 0.0;
        cL = 0.0;
        # // ringdown spherical
        # // these will be obtained by rescaling the 22 ringdown parameters (see Eq. (6.5) in the paper)
        alpha0 = 0.0;
        alpha2 = 0.0;
        alphaL = 0.0;
    
        
        
        # // set number of collocation points used (depends on the mode)
        nCollocationPts_inter=nCollocPtsInterPhase;
        # // define mass-ratio that discriminates between EMR and rest of cases
        eta_m1=1. /eta;
        
    
    
        #//initialize frequencies of colloc points in intermediate region

        #********************* INTERMEDIATE PHASE COLLOCATION POINTS *****************
        #function Intermediate_CollocPtsFreqs()

        #//The frequencies are stored in the struct pPhase, which is initialized with new values as the code starts to work on each mode.


        eta_factor=1.0+0.001*(0.25/eta-1.);

        fcut = eta_factor*emm*0.5*fMECO

        CollocationPointsFreqsPhaseInter = Vector{typeofFD}(undef, 6)
        
        CollocationPointsFreqsPhaseInter[1]=fcut;

        if ell_emm==32

            fRD22 = fring_22
            fdamp22 = fdamp_22
            fEnd= fRD22- 0.5*fdamp22;

            CollocationPointsFreqsPhaseInter[2]=(sqrt(3)*(fcut - fEnd) + 2*(fcut + fEnd))/4.;
            CollocationPointsFreqsPhaseInter[3]=(3*fcut + fEnd)/4.;
            CollocationPointsFreqsPhaseInter[4]=(fcut + fEnd)/2.;
            #// we use first and second derivative at fEnd, so this frequency is duplicated here
            CollocationPointsFreqsPhaseInter[5]= fEnd;
            CollocationPointsFreqsPhaseInter[6]= fEnd;
            fPhaseMatchIM = fEnd;
            #// correct cutting frequency for EMR with negative spins
            if eta<0.01 &&chi1<0 && IMRPhenomXHMIntermediatePhaseFreqsVersion==122019
                fPhaseMatchIM=fPhaseMatchIM*(1.2-0.25*chi1)
            end
        
        else

            CollocationPointsFreqsPhaseInter[2]=(sqrt(3)*(fcut - fring) + 2*(fcut + fring))/4.;
            CollocationPointsFreqsPhaseInter[3]=(3*fcut + fring)/4.;
            CollocationPointsFreqsPhaseInter[4]=(fcut + fring)/2.;
            CollocationPointsFreqsPhaseInter[5]=(fcut + 3*fring)/4.;
            CollocationPointsFreqsPhaseInter[6]=(fcut + 7*fring)/8.;
            fPhaseMatchIM = fring-fdamp;
    
        end

        fPhaseMatchIN = fMECOlm;

        # following part not used here
    
        # #********************* RINGDOWN PHASE COLLOCATION POINTS *****************

        # #// this function initializes the frequencies of the collocation points for the spheroidal ringdown reconstruction
        # # function Ringdown_CollocPtsFreqs()

        # #//The frequencies are stored in the struct pPhase, which is initialized with new values as the code starts to work on each mode.

        # fringlm = fring
        # fdamplm = fdamp
        # fring22 = fring_22

        # if IMRPhenomXHMRingdownPhaseFreqsVersion == 122019
            
        #     CollocationPointsFreqsPhaseRD = zeros(4)

        #     CollocationPointsFreqsPhaseRD[1] = fring22;
        #     CollocationPointsFreqsPhaseRD[3] = fringlm - 0.5*fdamplm;
        #     CollocationPointsFreqsPhaseRD[2] = fringlm - 1.5*fdamplm;
        #     CollocationPointsFreqsPhaseRD[4] = fringlm + 0.5*fdamplm;

        # else
        #     fdamp22 = fdamp_22;
        #     CollocationPointsFreqsPhaseRD = zeros(5)

        #     CollocationPointsFreqsPhaseRD[1] = fring22 - fdamp22;
        #     CollocationPointsFreqsPhaseRD[2] = fring22;
        #     CollocationPointsFreqsPhaseRD[3] = (fring22 + fringlm) * 0.5;
        #     CollocationPointsFreqsPhaseRD[4] = fringlm;
        #     CollocationPointsFreqsPhaseRD[5] = fringlm + fdamplm
        
        # end

        # ends here   

        #//for each collocation point, call fit giving the value of the phase derivative at that point
        CollocationPointsValuesPhaseInter = Vector{typeofFD}(undef, 6)
        CollocationPointsValuesPhaseInter = [p1, p2, p3, p4, p5, p6]
        CollocationPointsValuesPhaseInter .+= DeltaT
    
        fcutRD=fPhaseMatchIM;
        fcutInsp=fPhaseMatchIN;
    

        #************* Inspiral: rescale PhenX and apply PN corrections ********
    
        #//collect all the phenX inspiral coefficients
        # phi0, phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, 0., 0, 0, 0 from 22 # TaylorF2 PN Coefficients
        phenXnonLog = [Phase22.phi0,Phase22.phi1,Phase22.phi2,Phase22.phi3,Phase22.phi4,Phase22.phi5,Phase22.phi6,Phase22.phi7,Phase22.phi8,Phase22.phi9,0.,0.,0.,0.]
        phenXLog = [0.,0.,0.,0.,0.,Phase22.phi5L,Phase22.phi6L,0.,Phase22.phi8L,Phase22.phi9L,0,0,0,0]
        pseudoPN = [0.,0.,0.,0.,0.,0.,0.,0.,Phase22.sigma1,Phase22.sigma2,Phase22.sigma3,Phase22.sigma4,Phase22.sigma5]
        
        
        
        # //rescale the coefficients of phenX by applying phi_lm(f)~m/2 *phi_22(2/m f)
        # // for more details, see Appendix D of the paper
        m_over_2 = emm*0.5 
        two_over_m=1. /m_over_2;
    
        fact=phiNorm/eta;

        phi = Vector{typeofFD}(undef, 13)
        phiL = Vector{typeofFD}(undef, 13)


    
        for i in 1:N_MAX_COEFFICIENTS_PHASE_INS
    
            phenXnonLog[i]=phenXnonLog[i]*fact;
            phenXLog[i]=phenXLog[i]*fact;
            pseudoPN[i]=pseudoPN[i]*fact;
    
            #// scaling the logarithmic terms introduces some extra contributions in the non-log terms
            phi[i]=(phenXnonLog[i]+pseudoPN[i]-phenXLog[i]*log(m_over_2))* m_over_2^((8-(i-1))/3.);
            phiL[i]=phenXLog[i]* m_over_2^((8-(i-1))/3.);
        end
        
        
        # // if the mass-ratio is not extreme, use the complex PN amplitudes to correct the orbital phase at linear order in f
        # // see Eq. (4.12) in the paper
        if eta>0.01

            #// Returns linear-in-f term coming from the complex part of each mode's amplitude that will be added to the orbital phase, this is an expansion of Eq. (4.9) truncated at linear order
            function IMRPhenomXHM_Insp_Phase_LambdaPN(eta, ell_emm)

                if ell_emm==21
                    
                    #//2 f \[Pi] (-(1/2) - Log[16]/2)
                    output=(2. .*pi*(-0.5-2. .*log(2)));

                elseif ell_emm==33
                    #//2/3 f \[Pi] (-(21/5) + 6 Log[3/2])
                    output=2. ./3. .*pi*(-21. ./5. +6. .*log(1.5));

                elseif ell_emm==32
                    #{ //-((2376 f \[Pi] (-5 + 22 \[Eta]))/(-3960 + 11880 \[Eta]))
                    output=-((2376. .*pi*(-5. +22. .*eta))/(-3960. +11880*eta));

                elseif ell_emm==44
                    #{ //(45045 f \[powers_of_lalpiHM.itself] (336 - 1193 \[Eta] +320 (-1 + 3 \[Eta]) Log[2]))/(2 (-1801800 + 5405400 \[Eta]))
                    output=45045. .*pi*(336. .-1193. .*eta+320. .*(-1. +3. .*eta)*log(2))/(2. .*(-1801800. + 5405400. .*eta));
                        
                end

                return -1. .*output;
            end

            LambdaPN=IMRPhenomXHM_Insp_Phase_LambdaPN(eta, ell_emm);
        

        #// else use some phenomenological fits: the fits give the coefficient of a linear-in-f term to be added to the orbital phase
        else # Andrea to do
    
            LambdaPN= 0. #InspiralPhaseFits[pWFHM->modeInt](pWF22,pWFHM->IMRPhenomXHMInspiralPhaseVersion);
        end
    
        phi[9]= phi[9]+ LambdaPN;
    
    
        # ***************************************
        # *********** Intermediate-ringdown region ********
        # ***************************************
    
        # // below we set up the linear system to be solved in the intermediate region
    
    
        # GSL objects for solving system of equations via LU decomposition 

        b = Vector{typeofFD}(undef, nCollocationPts_inter)
        x = Vector{typeofFD}(undef, nCollocationPts_inter)
        A = Matrix{typeofFD}(undef, nCollocationPts_inter, nCollocationPts_inter)
    
    
        #// for high-spin cases: avoid sharp transitions for 21 mode by extending the inspiral region into the intermediate one
    
        if ell_emm==21 && STotR>=0.8
    
            insp_vals = Vector{typeofFD}(undef, 3)
            for i in 1:3
                FF=two_over_m*CollocationPointsFreqsPhaseInter[i];
                # _completePhaseDer is IMRPhenomX_dPhase_22
                insp_vals[i]=1. /eta* _completePhaseDer(model,FF, Phase22, fdamp_22, fring_22);
            end
    
            diff12=insp_vals[1]-insp_vals[2];
            diff23=insp_vals[2]-insp_vals[3];
    
            CollocationPointsValuesPhaseInter[2]=CollocationPointsValuesPhaseInter[3]+diff23;
            CollocationPointsValuesPhaseInter[1]=CollocationPointsValuesPhaseInter[2]+diff12;
    
        end
    
    
        # // choose collocation points according to spin/mass ratio
        # // current catalogue of simulations include some cases that create unphysical effects in the fits -> we need to use different subset of collocation points according to the parameters (we have to pick 5 out of 6 available fits)
        #  cpoints_indices is an array of integers labelling the collocation points chosen in each case, e.g.
        #  cpoints_indices={0,1,3,4,5} would mean that we are discarding the 3rd collocation points in the reconstructio 
    
        cpoints_indices = zeros(Int64,nCollocationPts_inter)
        
        cpoints_indices[1]=0;
        cpoints_indices[2]=1;
        cpoints_indices[5]=5;
    
    
        if eta < etaEMR || (emm==ell && STotR>=0.8) || (ell_emm==33 && STotR<0)
        
            cpoints_indices[3]=3;
            cpoints_indices[4]=4;
        
        elseif STotR>=0.8 && ell_emm==21
    
            cpoints_indices[3]=2;
            cpoints_indices[4]=4;
        
    
        else
            cpoints_indices[3]=2;
            cpoints_indices[4]=3;
        end
    
    
    
        if ModeMixingOn == false
            # // mode-mixing off: mode does not show significant mode-mixing
            # // the MixingOn flag is automatically switched off for the modes 21,33,44

    
            #   *********** Ringdown *************
    
            #    ansatz:
            #    alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)
    
            # //compute alpha2 by rescaling the coefficient of the 22-reconstruction
            if ell == emm
                wlm=2;
            else 
                wlm=emm/3.;
            end

            
            noSpin = 0.2088669311744758 - 0.37138987533788487*eta + 6.510807976353186*eta2 - 31.330215053905395*eta3 + 55.45508989446867*eta4;
            eqSpin = ((0.2393965714370633 + 1.6966740823756759*eta - 16.874355161681766*eta2 + 38.61300158832203*eta3)*S)/(1. - 0.633218538432246*S);
            uneqSpin = dchi*(0.9088578269496244*(eta)^(2.5) + 15.619592332008951*dchi*(eta)^(3.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            alpha2=1. /fring^2*wlm* total;
    
    
            #// compute alphaL
            
            noSpin = eta*(-1.1926122248825484 + 2.5400257699690143*eta - 16.504334734464244*eta2 + 27.623649807617376*eta3);
            eqSpin = eta3*S*(35.803988443700824 + 9.700178927988006*S - 77.2346297158916*S2) + eta*S*(0.1034526554654983 - 0.21477847929548569*S - 0.06417449517826644*S2) + eta2*S*(-4.7282481007397825 + 0.8743576195364632*S + 8.170616575493503*S2) + eta4*S*(-72.50310678862684 - 39.83460092417137*S + 180.8345521274853*S2);
            uneqSpin = (-0.7428134042821221*chi1*(eta)^(3.5) + 0.7428134042821221*chi2*(eta)^(3.5) + 17.588573345324154*(chi1)^(2)*(eta)^(4.5) - 35.17714669064831*chi1*chi2*(eta)^(4.5) + 17.588573345324154*(chi2)^(2)*(eta)^(4.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            alphaL=eta_m1*total;
    
            
    
            #// compute spherical-harmonic phase and its first derivative at matching frequency --used to glue RD and intermediate regions--
            phi0RD=RD_Phase_AnsatzInt(model, fcutRD, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)
            dphi0RD=RD_Phase_Ansatz(model, fcutRD, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)
    
    
            #*********** Intermediate *************
    
    
            #// set up the linear system by calling the fits of the intermediate-region collocation points
            for i in 1:nCollocationPts_inter
    
                ind=cpoints_indices[i] + 1 # Andrea: cpoints_indices is 0-based, ind is 1-based
                b[i] = CollocationPointsValuesPhaseInter[ind]
                ff = CollocationPointsFreqsPhaseInter[ind]
                ffm1=1. ./ff
                ffm2=ffm1*ffm1;
                fpowers = [1., (fdamp)/(fdamp^2+(ff-fring)^2),ffm1,ffm2,ffm2*ffm2, ffm1*ffm2]
                for j in 1:nCollocationPts_inter
                A[i,j] = fpowers[j]
                end
            end

        else
    
        
            #// mode-mixing on (the MixingOn flag is automatically switched on when reconstructing the 32 mode)
            
            # // for the 32 mode, the ringdown waveform is computed outside of the dedicated amplitude and phase routines, by calling SpheroidalToSpherical
            # // compute dphi/df and d2phi/df2 at fcutRD starting from rotation of spheroidal ansatz: these values will then enter the linear system that determines the coefficients in the intermediate region
    
            # // we compute derivatives by applying finite-difference schemes. In total we will need the value of the phase at three points around fcutRD
    
            # // we first compute the full spherical-harmonic-basis waveforms at three points
            function SpheroidalToSpherical_Der(fcutRD)
                
                SphericalWF = Vector{Union{ComplexF64, ForwardDiff.Dual{T} where T}}(undef, 3)
                fstep=0.0000001;
                
                for i in 1:3
        
                    FF=fcutRD+(i-2)*fstep;
                    SphericalWF[i]=SpheroidalToSpherical(FF)
                end
        
                phase_args = [
                    mod(angle(SphericalWF[1]), 2 * pi),
                    mod(angle(SphericalWF[2]), 2 * pi),
                    mod(angle(SphericalWF[3]), 2 * pi)
                ]
                # // make sure that all the three points belong to the same branch of mod
                for i in 1:3
        
                    if phase_args[i]>0
                        phase_args[i] -=2. .*pi
                    end
                end
            
                if debug == true
                    println("SphericalWF: ", SphericalWF)
                    println("phase_args: ", phase_args)
                    println("fcutRD: ", fcutRD)

                end

    
                # // store spherical-harmonic phase at fcutRD
                phi0RD=phase_args[1];
                fstep_m1= 1. ./fstep;
                # // we apply the FD schemes to get first and second phase derivatives at fcutRD
                dphi0RD=0.5*fstep_m1*(phase_args[3]-phase_args[1]);
                d2phi0RD=fstep_m1*fstep_m1*(phase_args[3]-2. .*phase_args[2]+phase_args[1]);
                return phi0RD, dphi0RD, d2phi0RD
            end


            function SpheroidalToSpherical_Der_FD(fcutRD)
                return ForwardDiff.derivative(x -> mod(angle(SpheroidalToSpherical(x)), 2 * pi), fcutRD)
            end

            function SpheroidalToSpherical_Der_FD2(fcutRD)
                return ForwardDiff.derivative(SpheroidalToSpherical_Der_FD, fcutRD)
            end


            #phi0RD, dphi0RD, d2phi0RD = SpheroidalToSpherical_Der(fcutRD)
            phi0RD = mod(angle(SpheroidalToSpherical(fcutRD)), 2 * pi)
            # in LAL they use a finite difference scheme to compute the derivative and need to move the phase to the correct branch
            # not needed here but we do it anyway to keep the code consistent
            if phi0RD > 0
                phi0RD -= 2 * pi
            end
            dphi0RD = SpheroidalToSpherical_Der_FD(fcutRD)
            d2phi0RD = SpheroidalToSpherical_Der_FD2(fcutRD)
            if debug == true
                println("phi0RD: ", phi0RD)
                println("dphi0RD: ", dphi0RD)
                println("d2phi0RD: ", d2phi0RD)
            end
            # // To achieve a smoother transition with the intermediate region (IR), we feed into the intermediate-region reconstruction the first and second derivative of the reconstructed ringdown phase, see Eq. (5.6)
    
            # // first derivative
            #   // if the mass-ratio is not extreme, we use the derivative computed using the FD scheme above, else we keep the original value of the fit
            #   // in the second case, we will therefore need to glue intermediate and ringdown region later
            if eta>etaEMR
                CollocationPointsFreqsPhaseInter[nCollocationPts_inter-2+1]=fcutRD;  # ANdrea +1 due to julia
                CollocationPointsValuesPhaseInter[nCollocationPts_inter-2+1]=dphi0RD;
            end
            # // second derivative
            CollocationPointsFreqsPhaseInter[nCollocationPts_inter-1+1]=fcutRD;
            CollocationPointsValuesPhaseInter[nCollocationPts_inter-1+1]=d2phi0RD;
    
    
            #// set up the linear system to determine the intermediate region coefficients
            for i in 1:nCollocationPts_inter
    
                b[i] = CollocationPointsValuesPhaseInter[i]
                ff=CollocationPointsFreqsPhaseInter[i] 
                ffm1=1. ./ff
                ffm2=ffm1*ffm1;
                fpowers = [1., fdamp / (fdamp^2 + (ff - fring)^2),ffm1,ffm2,ffm2*ffm2, ffm1*ffm2]
                for j in 1:nCollocationPts_inter
                    A[i,j] = fpowers[j]
                end
            end
    
            #// set the last point to equal the second derivative of the rotated spheroidal ansatz
            cpoint_ind=nCollocationPts_inter-1 +1 # Andrea +1 due to julia
            b[cpoint_ind] = CollocationPointsValuesPhaseInter[cpoint_ind]
            ff=CollocationPointsFreqsPhaseInter[cpoint_ind];
            ffm1=1. ./ff
            ffm2=ffm1*ffm1
            ffm3=ffm2*ffm1
            ffm4=ffm2*ffm2
            ffm5=ffm3*ffm2;
    
            fpowers = [0.,-2*(fdamp)*(ff-(fring))/(((fdamp)^(2)+(ff-(fring))^(2)))^(2),-ffm2, -2. .* ffm3, -4. .*ffm5,-3* ffm4]
            for j in 1:nCollocationPts_inter
                A[cpoint_ind,j] = fpowers[j]    
            end     
            
            if debug == true
                println("\n Derivative:")
                println("ff: ", CollocationPointsFreqsPhaseInter)
                println("b: ", b)
            end
    
        end


    
        
        # // at this point we have initizialized all the collocation points in the intermediate region, so we can solve to get the coefficients of the ansatz
        
        #  In the intermediate region , the most general ansatz is
        # (c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +cL fdamp/((fdamp)^2 + (f - fRD )^2))
    
        #  We now solve the system A x = b via an LU decomposition 


        x = lu(A) \ b
    
        c0 = x[1]
        cL = x[2]
        c1 = x[3]
        c2 = x[4]
        c4 = x[5]

    
        #// currently the 32 mode is calibrated using one extra point
        if ell_emm == 32
            c3 = x[6]
        end
    
        #// if the mass-ratio is extreme, we need to glue intermediate and ringdown region, as explained above
        if eta<etaEMR
        
            c0=c0+dphi0RD - Inter_Phase_Ansatz(fcutRD);
    

        end
    



        #************ PHASE ANSATZ ***************

        function Inter_Phase_AnsatzInt(ff)

            invf  = 1. /ff;
            invf2 = invf*invf;
            invf3 = invf2*invf;
            logfv= log(ff);

            fda=fdamp;
            frd=fring;

            if ell_emm!=32
            #  5 coefficients 

                phaseIR = c0 *ff + (c1)*logfv - (c2)*invf -1. ./3. * (c4)*invf3 + (cL)* atan((ff - frd)/fda);
            else           #  6 coefficients 

                invf3 = invf*invf2;
                # (c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +
                # cL fdamp/((fdamp)^2 + (f - fRD )^2))
                phaseIR = c0 *ff + (c1)*logfv - (c2)*invf -1. ./3. * (c4)*invf3 -0.5*c3*invf2+ (cL)* atan((ff - frd)/fda);
            end

            return phaseIR;


        end


        function Inter_Phase_Ansatz(ff)

            invf  = 1. ./ff;
            invf2 = invf*invf;
            invf4 = invf2*invf2;
            fda=fdamp;
            frd=fring;

            if ell_emm!=32
            #  5 coefficients 

            #(c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 +
            #cL fdamp/((fdamp)^2 + (f - fRD )^2))
            dphaseIR = ( c0 + (c1)*invf + (c2)*invf2 + (c4)*invf4 + ( (cL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            else #            6 coefficients 

                invf3 = invf*invf2;
                # (c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +
                # cL fdamp/((fdamp)^2 + (f - fRD )^2))
                dphaseIR = ( c0 + (c1)*invf + (c2)*invf2 + (c4)*invf4 + (c3)*invf3 + ( (cL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            end

            return dphaseIR

            
        end

        #************ INSPIRAL PHASE ANSATZ ************

        function Inspiral_Phase_AnsatzInt(Mf)
            

            # 
            # The function loads the rescaled phenX coefficients for the phase and correct the result with a linear-in-f term coming from the complex part of the PN amplitude

            # The rescaling of the 22-coefficients is explained in App. D
            # 
            
            #//compute the orbital phase by laoding the rescaled phenX coefficients
            philm=0.;
            freqs = [ Mf^(-5/3), Mf^(-4/3), Mf^(-1), Mf^(-2/3), Mf^(-1/3), 1., Mf^(1/3), Mf^(2/3), Mf, Mf^(4/3), Mf^(5/3), Mf^2, Mf^(7/3)]
            logMf= log(Mf);
            for i in 1:N_MAX_COEFFICIENTS_PHASE_INS
                philm+=(phi[i]+phiL[i]*(logMf))*freqs[i];
            end
            return philm;
            
        end     
        
        
        function Inspiral_Phase_Ansatz(Mf)


            #//compute the orbital phase by laoding the rescaled phenX coefficients
            dphilm=0.;
            coeffs = [-5. ./3,-4. ./3,-1.,-2. ./3,-1. ./3,0.,1. ./3, 2. ./3, 1., 4. ./3, 5. ./3, 2., 7. ./3]
            freqs = [Mf^(-8 /3), Mf^(-7 /3), Mf^(-2), Mf^(-5 /3), Mf^(-4 /3), Mf^(-1), Mf^(-2 /3), Mf^(-1 /3), 1., Mf^(1 /3), Mf^(2 /3), Mf, Mf^(4 /3)]
                
            logMf= log(Mf);

            for i in 1:N_MAX_COEFFICIENTS_PHASE_INS
                dphilm+=((phi[i]+phiL[i]*(logMf))*coeffs[i]+phiL[i])*freqs[i];
            end
            return dphilm;
        
        end
    
        #***************** end of ringdown and inspiral reconstruction *************
    
        #//glue inspiral and intermediate

        C1INSP = Inter_Phase_Ansatz(fcutInsp) - Inspiral_Phase_Ansatz(fcutInsp)
        CINSP = -(C1INSP*fcutInsp)+Inter_Phase_AnsatzInt(fcutInsp) - Inspiral_Phase_AnsatzInt(fcutInsp)
        # //glue ringdown and intermediate regions
    
        C1RD=Inter_Phase_Ansatz(fcutRD)-dphi0RD;
        CRD=-C1RD*fcutRD+Inter_Phase_AnsatzInt(fcutRD)-phi0RD;
        
        
        # // we now have a C1 reconstruction of the phase
        # // below we align each mode so that at low-f its relative phase wrt the 22 agrees with PN
        # //somehow arbitary cutoff to pick frequency for alignment: must be in the inspiral region
        if eta>etaEMR
            falign=0.6*m_over_2*fMECO;
        else
            falign=m_over_2*fMECO;
        end

    
    
        # compute explicitly the phase normalization to be applied to IMRPhenomX, when mode-mixing is on this will have been already computed in GetSpheroidalCoefficients 
        tmp = _completePhase(model,[MfRef,two_over_m*falign], Phase22, fdamp_22, fring_22)
        phi_MfRef = tmp[1]
        phi_falign = tmp[2]

        #phi_MfRef = _completePhase(fRef)
        if ModeMixingOn== false

            # Phase22 = Phase_22_ConnectionCoefficients(mc, eta, chi1, chi2)
            timeshift= TimeShift_22(model, eta, Seta, S, dchi, fring_22, fdamp_22, Phase22)
            phiref22 = -1 ./eta*phi_MfRef - timeshift*MfRef - phaseshift + 2.0*phi0 + pi/4.0 

        end
    
        #// we determine the phase normalization of each mode by imposing Eq. (4.13), i.e. phi_lm(fref)-m/2 phi_22(2/m fref)~3. .*pi_4*(1-m_over_2)
        deltaphiLM = m_over_2*(1 ./eta*phi_falign + phaseshift + phiref22)+timeshift*falign-3. .*pi/4. *(1-m_over_2)-(Inspiral_Phase_AnsatzInt(falign)+C1INSP*falign+CINSP);
        deltaphiLM=mod(deltaphiLM, 2. .*pi);
    
    
        #// for the 21, we need to make sure the sign of the PN amplitude is positive, else we'll need to flip its phase by Pi
        if ell_emm == 21
            
            ff = 0.008
            output=(-16*delta*eta*ff*(pi)^(1.5))/(3. *sqrt(5)) + (4*(2)^(0.3333333333333333)*(chi1 - chi2 + delta *(chi1+ chi2))*eta*(ff)^(1.3333333333333333)*(pi)^(1.8333333333333333))/sqrt(5) + (2*(2)^(0.6666666666666666)*eta*(306*delta - 360*delta*eta)*(ff)^(1.6666666666666667)*(pi)^(2.1666666666666665))/(189. *sqrt(5));

            if output>=0  
                ampsign=1;
            else 
                ampsign=-1;
            end

            
            #// the sign of the 21 amplitude changes across the parameter space, we add Pi to the phase to account for this - remember that the amplitude of the model is positive by construction
            if ampsign>0
                deltaphiLM+=pi;
            end
    
        end


        ##### end of IMRPhenomXHM_GetPhaseCoefficients

        function IMRPhenomXHMFDAddMode(htildelm)
            # Function to sum one mode (htildelm) to hp/c tilde 
            if ell % 2 !=0
        
                minus1l = -1; # (-1)^l 
        
            else
        
                minus1l = 1;
        
            end
        
            Ym = _spinWeightedSphericalHarmonic(model, iota, ell_emm, -emm)    # since _spinWeightedSphericalHarmonic() is first called with -m
                
            # Equatorial symmetry: add in -m and m mode 
            Ystar = conj(_spinWeightedSphericalHarmonic(model, iota, ell_emm, emm))
            factorp = 0.5 * (Ym + minus1l * Ystar);
            factorc = 1im * 0.5 * ( Ym - minus1l * Ystar);
        
        
            for j in eachindex(htildelm)    
              hlm = htildelm[j];
              hp[j] += (factorp * hlm);
              hc[j] += (factorc * hlm);
            end

            if debug == true
                println()
                println("ell_emm = $ell_emm")
                println("Ym = $Ym")
                println("Ystar = $Ystar")
                println("factorp = $factorp")
                println("factorc = $factorc")
                println("hlm[1] =", htildelm[1])
            end
        
        end
        
        
        
        # ************************************
        #                                     
        #   RECONSTRUCTION FUNCTIONS THROUGH  
        #      INSPIRAL, MERGER & RINGDOWN    
        #                                     
        # ************************************
        
        #  These functions return the amplitude/phase for a given input frequency.
        #    They are piece-wise functions that distiguish between inspiral, intermediate and ringdown
        #    and call the corresponding reconstruction function.
        # 
        
        # *****************************************************************************
        #   Compute IMRPhenomXHM AMPLITUDE given an amplitude coefficients struct pAmp 
        # *****************************************************************************
        
        # // WITHOUT mode mixing. It returns the whole amplitude (in NR units) without the normalization factor of the 22: sqrt[2 * eta / (3 * pi^(1/3))]
        function Amplitude_noModeMixing(Mf)
          #// If it is an odd mode and equal black holes case this mode is zero.
            if Ampzero==1
                return 0.;
            end
            
            #// Use step function to only calculate IMR regions in approrpiate frequency regime
            #// Inspiral range
            if Mf < fAmpMatchIN
                Amp = Inspiral_Amp_Ansatz(Mf);
            
            #// MRD range
            elseif Mf > fAmpMatchIM        
                Amp = RD_Amp_Ansatz(Mf)
                if  IMRPhenomXHMRingdownAmpVersion == 0
                    Amp *= Mf^(-7/6) * ampNorm;
                end

            else #  //Second intermediate region
                Amp = Intermediate_Amp_Ansatz(Mf)
            end
            return Amp
        end
           
        
        #// If the 22 mode has been previously computed, we use it here for the rotation.
        function SpheroidalToSphericalRecycle(Mf,idx)
        
            # // The input 22 in the whole 22, and for the rotation we have to rescaled with the leading order. This is because the 32 is also rescaled.
            # //complex wf22R = wf22/(powers_of_f->m_seven_sixths * pWF22->amp0);
            wf22R = wf22[idx] / Amp0;

            if IMRPhenomXHMRingdownAmpVersion == 0
                wf22R /= ampNorm * Mf^(-7/6);
            end
            #// Compute 32 mode in spheroidal.
            amplm=RD_Amp_Ansatz(Mf)#, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)
            philm=RD_Phase_AnsatzInt(model, Mf, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)
            #// Do the rotation
            sphericalWF_32 =  conj(mixingCoeffs[3])*wf22R  +  conj(mixingCoeffs[4])*amplm*exp(1im*philm);
            return sphericalWF_32;
        
        end
        
        
        function Amplitude_ModeMixingRecycle(Mf,idx)
        
            #// WITH mode mixing and recycling the previously computed 22 mode. It returns the whole amplitude (in NR units) without the normalization factor of the 22: sqrt[2 * eta / (3 * pi^(1/3))].
              #// Use step function to only calculate IMR regions in approrpiate frequency regime
              #// Inspiral range
            if Mf < fAmpMatchIN
                Amp =  Inspiral_Amp_Ansatz(Mf)
            
                #// MRD range
            elseif Mf > fAmpMatchIM
                Amp = abs(SpheroidalToSphericalRecycle(Mf,idx));
                if IMRPhenomXHMRingdownAmpVersion == 0 
                    Amp *= Mf^(-7/6) * ampNorm;
                end

            else
                Amp = Intermediate_Amp_Ansatz(Mf);
            end

            if Amp < 0 && IMRPhenomXHMReleaseVersion != 122019
                Amp = 0.;
            end
        
            return Amp;
        
        end
        
        #    / WITH mode mixing and recycling the previously computed 22 mode.
        function Phase_ModeMixingRecycle(Mf,idx)
          #// Inspiral range, f < fPhaseInsMax
            if Mf < fPhaseMatchIN
          
                PhiIns = Inspiral_Phase_AnsatzInt(Mf);
                return PhiIns + C1INSP*Mf + CINSP + deltaphiLM;
            end
            #// MRD range, f > fPhaseIntMax
            if Mf > fPhaseMatchIM
          
                PhiMRD = angle(SpheroidalToSphericalRecycle(Mf,idx));
                return PhiMRD + C1RD*Mf + CRD + deltaphiLM;
            end
        
            #//Intermediate range, fPhaseInsMax < f < fPhaseIntMax
            PhiInt = Inter_Phase_AnsatzInt(Mf);
            return PhiInt + deltaphiLM;
        
        end
        
        
        function Phase_noModeMixing(Mf)
            #// WITHOUT mode mixing.
            
            #// Inspiral range, f < fPhaseInsMax
            if Mf < fPhaseMatchIN
                PhiIns = Inspiral_Phase_AnsatzInt(Mf)
                return PhiIns + C1INSP*Mf + CINSP + deltaphiLM;
            end
            #// MRD range, f > fPhaseIntMax
            if Mf > fPhaseMatchIM
                PhiMRD = RD_Phase_AnsatzInt(model, Mf, alpha0, alpha2, alphaL, fring, fdamp, ModeMixingCoeffs,  ModeMixingOn)
                return PhiMRD + C1RD*Mf + CRD + deltaphiLM;
            end
            #//Intermediate range, fPhaseInsMax < f < fPhaseIntMax
            PhiInt = Inter_Phase_AnsatzInt(Mf);
            return PhiInt + deltaphiLM;
            
        end
        
        
        
        # *******************************************
        #              INTERMEDIATE AMPLITUDE ANSATZ     
        #     *******************************************
            
        #     // Build the polynomial with the coefficients given and return the inverse of the polynomial (this is the ansatz)
        function Intermediate_Amp_Ansatz(Mf)
        
            result = 0.
            fpower = 1.;
            # Ansatz = f^(-7/6) * polynomial 
            for i in 1:nCoefficientsInter
                result += InterCoefficient[i] * fpower;
                fpower *= Mf;
            end
            result *= Mf^(-7/6);
            return result;
                
        end



        #  Loop over frequencies to generate waveform 
        #  With mode mixng 
        if debug == true
            println("MODE $ell $emm")
            println("Mf amp phi")
        end
        if ModeMixingOn==true
          #// If the 22 mode has been already computed we use it for the mixing of the 32.
            for idx in 1:len
            
                #wf22 = data[idx + offset]; #//This will be rescaled inside SpheroidalToSphericalRecycle for the rotation

                amp = Amplitude_ModeMixingRecycle(Mf[idx],idx);
                Phi = Phase_ModeMixingRecycle(Mf[idx],idx);
                
                
                if debug == true
                    fff = Mf[idx]
                    ccc = Amp0 * amp
                    println("data$idx = [$fff, $ccc, $Phi]")
                end
                # Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] 
                htildelm__[idx+offset] = Amp0 * amp * exp(1im * Phi);
            end
          
          #// If the 22 has not been computed, its ringdown part is computed internally using pAmp22 and pPhase22.
           #   No mode mixing 
        else
            for idx in 1:len
          
                amp = Amplitude_noModeMixing(Mf[idx])
                Phi = Phase_noModeMixing(Mf[idx])
                
                if debug == true
                    fff = Mf[idx]
                    ccc = Amp0 * amp
                    println("data$idx = [$fff, $ccc, $Phi]")
                end
                /#* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] 
                htildelm__[idx+offset] = Amp0 * amp * exp(1im * Phi);
            end
        end # End of loop over frequencies 

        

        if debug == true 
            println("fring: ", fring)
            println("fdamp: ", fdamp)
            println("fAmpMatchIN: ", fAmpMatchIN)
            println("fAmpMatchIM: ", fAmpMatchIM)
            println("fPhaseMatchIN: ", fPhaseMatchIN)
            println("fPhaseMatchIM: ", fPhaseMatchIM)
            println("alpha2: " , alpha2)
            println("alphaL: " , alphaL)
            println("c0: ", c0) 
            println("c1: ", c1)
            println("c2: ", c2)
            println("c4: ", c4)
            println("cL: ", cL)

            println()
            println("LambdaPN = ", LambdaPN)

            inter_fit = [p1, p2, p3, p4, p5, p6]
            for i in 1:6
                println("CollocationPointsFreqsPhaseInter[$i] = ", CollocationPointsFreqsPhaseInter[i])
                println("CollocationPointsValuesPhaseInter[$i] = ", inter_fit[i])
            end

            println()
            println("alambda: ", RDCoefficient[1])
            println("lambda: ", RDCoefficient[2])
            println("sigma: ", RDCoefficient[3])
            println("lc: ", RDCoefficient[4])
            println("RDCoeff5: ", RDCoefficient[5])

            println()

            println("fInsp1: ", CollocationPointsFreqsAmplitudeInsp[1])
            println("fInsp2: ", CollocationPointsFreqsAmplitudeInsp[2])
            println("fInsp3: ", CollocationPointsFreqsAmplitudeInsp[3])

            println("fInt1: ", CollocationPointsFreqsAmplitudeInter[1])
            println("fInt2: ", CollocationPointsFreqsAmplitudeInter[2])
    

            println("A_fInsp1: ", Insp_v1)
            println("A_fInsp2: ", Insp_v2)
            println("A_fInsp3: ", Insp_v3)
            println("A_fInt1: ", CollocationPointsValuesAmplitudeInter[1])
            println("A_fInt2: ", CollocationPointsValuesAmplitudeInter[2])


            println("CollocationPointsValuesAmplitudeInter: ", CollocationPointsValuesAmplitudeInter)


            println()

            println("PN_f1: ", PNAmplitudeInsp[1])
            println("PN_f2: ", PNAmplitudeInsp[2])
            println("PN_f3: ", PNAmplitudeInsp[3])
        

            println()

            println("C1INSP: ", C1INSP)
            println("CINSP: ", CINSP)
            println("C1RD: ", C1RD)
            println("CRD: ", CRD)
            println()
            println("timeshift: ", timeshift)
            println("deltaphiLM: ", deltaphiLM)
            println("phaseshift: ", phaseshift)
            println("phiref22: ", phiref22)

            println()
        end

        if return_modes
            coprecessing_modes[(ell, -emm)] = copy(htildelm__)
        end

        IMRPhenomXHMFDAddMode(htildelm__); #// add both positive and negative modes

        if debug ==true && ell_emm == 21

            println(" Parameters already checked: ")
            println("a0coloc: ", Phase22.a0coloc)
            println("a1coloc: ", Phase22.a1coloc)
            println("a2coloc: ", Phase22.a2coloc)
            println("a3coloc: ", Phase22.a3coloc)
            println("a4coloc: ", Phase22.a4coloc)

            println("b0coloc: ", Phase22.b0coloc)
            println("b1coloc: ", Phase22.b1coloc)
            println("b2coloc: ", Phase22.b2coloc)
            println("b3coloc: ", Phase22.b3coloc)
            println("b4coloc: ", Phase22.b4coloc)

            println("c0coloc: ", Phase22.c0coloc)
            println("c1coloc: ", Phase22.c1coloc)
            println("c2coloc: ", Phase22.c2coloc)
            println("cLcoloc: ", Phase22.cLcoloc)
            println("c4coloc: ", Phase22.c4coloc)

            println("C2Int: ", Phase22.C2Int)
            println("C2MRD: ", Phase22.C2MRD)

            println()
            println("TaylorF2 PN Coefficients: ")
            println("phi0: ", Phase22.phi0)
            println("phi1: ", Phase22.phi1)
            println("phi2: ", Phase22.phi2)
            println("phi3: ", Phase22.phi3)
            println("phi4: ", Phase22.phi4)
            println("phi5: ", Phase22.phi5)
            println("phi6: ", Phase22.phi6)
            println("phi7: ", Phase22.phi7)
            println("phi8: ", Phase22.phi8)
            println("phi5L: ", Phase22.phi5L)
            println("phi6L: ", Phase22.phi6L)
            println("phi8L: ", Phase22.phi8L)

            println("phi8P: ", Phase22.sigma1)
            println("phi9P: ", Phase22.sigma2)
            println("phi10P: ", Phase22.sigma3)
            println("phi11P: ", Phase22.sigma4)
            println("phi12P: ", Phase22.sigma5)

            println("TaylorF2 PN Derivative Coefficients: ")
            println("dphi0: ", Phase22.dphi0)
            println("dphi1: ", Phase22.dphi1)
            println("dphi2: ", Phase22.dphi2)
            println("dphi3: ", Phase22.dphi3)
            println("dphi4: ", Phase22.dphi4)
            println("dphi5: ", Phase22.dphi5)
            println("dphi6: ", Phase22.dphi6)
            println("dphi7: ", Phase22.dphi7)
            println("dphi8: ", Phase22.dphi8)
            println("dphi9: ", Phase22.dphi9)
            println("dphi6L: ", Phase22.dphi6L)
            println("dphi8L: ", Phase22.dphi8L)
            println("dphi9L: ", Phase22.dphi9L)
        end
        
    end # End of loop over modes 




    return return_modes ? coprecessing_modes : (hp, hc)
end


function _completePhaseDer(model, infreqs, Phase22, fdamp, fring)
    return @. ifelse(infreqs <= Phase22.fPhaseMatchIN, (infreqs^(-8. /3.))*Phase22.dphase0*(Phase22.dphi0 + Phase22.dphi1*(infreqs^(1. /3.)) + Phase22.dphi2*(infreqs^(2. /3.)) + Phase22.dphi3*infreqs + Phase22.dphi4*(infreqs^(4. /3.)) + Phase22.dphi5*(infreqs^(5. /3.)) + (Phase22.dphi6 + Phase22.dphi6L*log(infreqs))*infreqs*infreqs + Phase22.dphi7*(infreqs^(7. /3.)) + (Phase22.dphi8 + Phase22.dphi8L*log(infreqs))*(infreqs^(8. /3.)) + (Phase22.dphi9  + Phase22.dphi9L*log(infreqs))*infreqs*infreqs*infreqs + Phase22.a0coloc*(infreqs^(8. /3.)) + Phase22.a1coloc*infreqs*infreqs*infreqs + Phase22.a2coloc*(infreqs^(10. /3.)) + Phase22.a3coloc*(infreqs^(11. /3.)) + Phase22.a4coloc*(infreqs^4)), ifelse(infreqs <= Phase22.fPhaseMatchIM, Phase22.b0coloc + Phase22.b1coloc/infreqs + Phase22.b2coloc/(infreqs*infreqs) + Phase22.b3coloc/(infreqs*infreqs*infreqs) + Phase22.b4coloc/(infreqs*infreqs*infreqs*infreqs) + (4. *Phase22.cLcoloc) / ((4. *fdamp*fdamp) + (infreqs - fring)*(infreqs - fring)) + Phase22.C2Int, (Phase22.c0coloc + Phase22.c1coloc*(infreqs^(-1. /3.)) + Phase22.c2coloc/(infreqs*infreqs) + Phase22.c4coloc/(infreqs*infreqs*infreqs*infreqs) + (Phase22.cLcoloc / (fdamp*fdamp + (infreqs - fring)*(infreqs - fring)))) + Phase22.C2MRD))
end

function _completePhase(model, infreqs, Phase22, fdamp, fring, fcutPar = 0.3)
    phiNorm = - (3. * (pi^(-5. /3.)))/ 128.
    c4ov3   = Phase22.c4coloc / 3.
    cLovfda = Phase22.cLcoloc / fdamp
    return @. ifelse(infreqs <= Phase22.fPhaseMatchIN, phiNorm*(infreqs^(-5. /3.))*(Phase22.phi0 + Phase22.phi1*(infreqs^(1. /3.)) + Phase22.phi2*(infreqs^(2. /3.)) + Phase22.phi3*infreqs + Phase22.phi4*(infreqs^(4. /3.)) + (Phase22.phi5 + Phase22.phi5L*log(infreqs))*(infreqs^(5. /3.)) + (Phase22.phi6 + Phase22.phi6L*log(infreqs))*infreqs*infreqs + Phase22.phi7*(infreqs^(7. /3.)) + (Phase22.phi8 + Phase22.phi8L*log(infreqs))*(infreqs^(8. /3.)) + (Phase22.phi9  + Phase22.phi9L*log(infreqs))*infreqs*infreqs*infreqs + Phase22.sigma1*(infreqs^(8. /3.)) + Phase22.sigma2*(infreqs*infreqs*infreqs) + Phase22.sigma3*(infreqs^(10. /3.)) + Phase22.sigma4*(infreqs^(11. /3.)) + Phase22.sigma5*(infreqs^4)), ifelse(infreqs <= Phase22.fPhaseMatchIM, Phase22.b0coloc*infreqs + Phase22.b1coloc*log(infreqs) - Phase22.b2coloc/infreqs - Phase22.b3coloc/(infreqs*infreqs)/2. - (Phase22.b4coloc/(infreqs*infreqs*infreqs)/3.) + (2. * Phase22.cLcoloc * atan((infreqs - fring) / (2. * fdamp)))/fdamp + Phase22.C1Int + Phase22.C2Int*infreqs, ifelse(infreqs < fcutPar, (Phase22.c0coloc*infreqs + 1.5*Phase22.c1coloc*(infreqs^(2. /3.)) - Phase22.c2coloc/infreqs - c4ov3/(infreqs*infreqs*infreqs) + (cLovfda * atan((infreqs - fring)/fdamp))) + Phase22.C1MRD + Phase22.C2MRD*infreqs, 0.)))
end

function TimeShift_22(model, eta, Seta, totchi, dchi, fring_22, fdamp_22, Phase22)

    eta2 = eta^2
    totchi2 = totchi^2
    etaInv = 1. / eta
    delta = Seta
    lina = 0.
    linb         = ((3155.1635543201924 + 1257.9949740608242*eta - 32243.28428870599*eta2 + 347213.65466875216*eta2*eta - 1.9223851649491738e6*eta2*eta2 + 5.3035911346921865e6*eta2*eta2*eta - 5.789128656876938e6*eta2*eta2*eta2) + ((-24.181508118588667 + 115.49264174560281*eta - 380.19778216022763*eta2)*totchi + (24.72585609641552 - 328.3762360751952*eta + 725.6024119989094*eta2)*totchi2 + (23.404604124552 - 646.3410199799737*eta + 1941.8836639529036*eta2)*totchi2*totchi + (-12.814828278938885 - 325.92980012408367*eta + 1320.102640190539*eta2)*totchi2*totchi2) + (-148.17317525117338*dchi*delta*eta2))
    dphi22Ref    = etaInv * _completePhaseDer(model, fring_22-fdamp_22, Phase22, fdamp_22, fring_22)
    psi4tostrain = ((13.39320482758057 - 175.42481512989315*eta + 2097.425116152503*eta2 - 9862.84178637907*eta2*eta + 16026.897939722587*eta2*eta2) + ((4.7895602776763 - 163.04871764530466*eta + 609.5575850476959*eta2)*totchi + (1.3934428041390161 - 97.51812681228478*eta + 376.9200932531847*eta2)*totchi2 + (15.649521097877374 + 137.33317057388916*eta - 755.9566456906406*eta2)*totchi2*totchi + (13.097315867845788 + 149.30405703643288*eta - 764.5242164872267*eta2)*totchi2*totchi2) + (105.37711654943146*dchi*Seta*eta2))
    
    tshift         = linb - dphi22Ref - 2. *pi*(500. +psi4tostrain)
    # // here we align the model to the hybrids, for which psi4 peaks 500M before the end of the waveform
    # // linb is a parameter-space fit of dphi22(fring22-fdamp22), evaluated on the calibration dataset

    return tshift

end
