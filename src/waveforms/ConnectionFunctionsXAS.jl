# Andrea needs to fill the structure Phase22

function Phase_22_ConnectionCoefficients(mc,
    eta,
    chi1,
    chi2;
    fcutPar=0.3,  # ANDREA in phenomD is 0.2 
    GMsun_over_c3= uc.GMsun_over_c3,
    fInsJoin_PHI = 0.018,
    InsPhaseVersion=104,
    IntPhaseVersion=105,
    final_spin_override=nothing
)  
    

    
    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    pi2 = pi * pi


    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)

    # We work in dimensionless frequency M*f, not f

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    q = 0.5 * (1.0 + Seta - 2.0 * eta) / eta
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    m1ByMSq = m1ByM * m1ByM
    m2ByMSq = m2ByM * m2ByM
    # PN symmetry coefficient
    delta = Seta
    totchi = ((m1ByMSq * chi1 + m2ByMSq * chi2) / (m1ByMSq + m2ByMSq))
    totchi2 = totchi * totchi
    dchi = chi1 - chi2
    # Normalised PN reduced spin parameter
    chi_eff = (m1ByM * chi1 + m2ByM * chi2)
    chiPN = (chi_eff - (38.0 / 113.0) * eta * (chi1 + chi2)) / (1.0 - (76.0 * eta / 113.0))
    chiPN2 = chiPN * chiPN

    dphase0 = 5.0 / (128.0 * (pi^(5.0 / 3.0)))

    gpoints4 = [0.0, 1.0 / 4.0, 3.0 / 4.0, 1.0]
    gpoints5 = [
        0.0,
        1.0 / 2.0 - 1.0 / (2.0 * sqrt(2.0)),
        1.0 / 2.0,
        1.0 / 2 + 1.0 / (2.0 * sqrt(2.0)),
        1.0,
    ]

    # Compute final spin and radiated energy
    aeff   = (((3.4641016151377544*eta + 20.0830030082033*eta2 - 12.333573402277912*eta2*eta)/(1 + 7.2388440419467335*eta)) + ((m1ByMSq + m2ByMSq)*totchi + ((-0.8561951310209386*eta - 0.09939065676370885*eta2 + 1.668810429851045*eta2*eta)*totchi + (0.5881660363307388*eta - 2.149269067519131*eta2 + 3.4768263932898678*eta2*eta)*totchi2 + (0.142443244743048*eta - 0.9598353840147513*eta2 + 1.9595643107593743*eta2*eta)*totchi2*totchi) / (1 + (-0.9142232693081653 + 2.3191363426522633*eta - 9.710576749140989*eta2*eta)*totchi)) + (0.3223660562764661*dchi*Seta*(1 + 9.332575956437443*eta)*eta2 - 0.059808322561702126*dchi*dchi*eta2*eta + 2.3170397514509933*dchi*Seta*(1 - 3.2624649875884852*eta)*eta2*eta*totchi))
    aeff_nonprec = aeff
    aeff = isnothing(final_spin_override) ? aeff : final_spin_override
    Erad   = ((((0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta2*eta + 3.145145224278187*eta2*eta2)*(1. + (-0.13084389181783257 - 1.1387311580238488*eta + 5.49074464410971*eta2)*totchi + (-0.17762802148331427 + 2.176667900182948*eta2)*totchi2 + (-0.6320191645391563 + 4.952698546796005*eta - 10.023747993978121*eta2)*totchi*totchi2)) / (1. + (-0.9919475346968611 + 0.367620218664352*eta + 4.274567337924067*eta2)*totchi)) + (- 0.09803730445895877*dchi*Seta*(1. - 3.2283713377939134*eta)*eta2 + 0.01118530335431078*dchi*dchi*eta2*eta - 0.01978238971523653*dchi*Seta*(1. - 4.91667749015812*eta)*eta*totchi))
    # Compute ringdown and damping frequencies from fits
    fring = ((0.05947169566573468 - 0.14989771215394762*aeff + 0.09535606290986028*aeff*aeff + 0.02260924869042963*aeff*aeff*aeff - 0.02501704155363241*aeff*aeff*aeff*aeff - 0.005852438240997211*(aeff^5) + 0.0027489038393367993*(aeff^6) + 0.0005821983163192694*(aeff^7))/(1 - 2.8570126619966296*aeff + 2.373335413978394*aeff*aeff - 0.6036964688511505*aeff*aeff*aeff*aeff + 0.0873798215084077*(aeff^6)))/(1. - Erad)
    fdamp = ((0.014158792290965177 - 0.036989395871554566*aeff + 0.026822526296575368*aeff*aeff + 0.0008490933750566702*aeff*aeff*aeff - 0.004843996907020524*aeff*aeff*aeff*aeff - 0.00014745235759327472*(aeff^5) + 0.0001504546201236794*(aeff^6))/(1 - 2.5900842798681376*aeff + 1.8952576220623967*aeff*aeff - 0.31416610693042507*aeff*aeff*aeff*aeff + 0.009002719412204133*(aeff^6)))/(1. - Erad)
    
    aeff = aeff_nonprec
    # Fitting function for hybrid minimum energy circular orbit (MECO) function and computation of ISCO frequency
    Z1tmp = 1. + cbrt((1. - aeff*aeff) ) * (cbrt(1. + aeff) + cbrt(1. - aeff))
    Z1tmp = ifelse(Z1tmp>3., 3., Z1tmp)
    Z2tmp = sqrt(3. *aeff*aeff + Z1tmp*Z1tmp)
    fISCO  = (1. / ((3. + Z2tmp - sign(aeff)*sqrt((3. - Z1tmp) * (3. + Z1tmp + 2. *Z2tmp)))^(3. /2.) + aeff))/pi
    
    fMECO    = (((0.018744340279608845 + 0.0077903147004616865*eta + 0.003940354686136861*eta2 - 0.00006693930988501673*eta2*eta)/(1. - 0.10423384680638834*eta)) + ((chiPN*(0.00027180386951683135 - 0.00002585252361022052*chiPN + eta2*eta2*(-0.0006807631931297156 + 0.022386313074011715*chiPN - 0.0230825153005985*chiPN2) + eta2*(0.00036556167661117023 - 0.000010021140796150737*chiPN - 0.00038216081981505285*chiPN2) + eta*(0.00024422562796266645 - 0.00001049013062611254*chiPN - 0.00035182990586857726*chiPN2) + eta2*eta*(-0.0005418851224505745 + 0.000030679548774047616*chiPN + 4.038390455349854e-6*chiPN2) - 0.00007547517256664526*chiPN2))/(0.026666543809890402 + (-0.014590539285641243 - 0.012429476486138982*eta + 1.4861197211952053*eta2*eta2 + 0.025066696514373803*eta2 + 0.005146809717492324*eta2*eta)*chiPN + (-0.0058684526275074025 - 0.02876774751921441*eta - 2.551566872093786*eta2*eta2 - 0.019641378027236502*eta2 - 0.001956646166089053*eta2*eta)*chiPN2 + (0.003507640638496499 + 0.014176504653145768*eta + 1. *eta2*eta2 + 0.012622225233586283*eta2 - 0.00767768214056772*eta2*eta)*chiPN2*chiPN)) + (dchi*dchi*(0.00034375176678815234 + 0.000016343732281057392*eta)*eta2 + dchi*Seta*eta*(0.08064665214195679*eta2 + eta*(-0.028476219509487793 - 0.005746537021035632*chiPN) - 0.0011713735642446144*chiPN)))
    
    
    fIMmatch = 0.6 * (0.5 * fring + fISCO)
    fINmatch = fMECO
    deltaf   = (fIMmatch - fINmatch) * 0.03
    fPhaseMatchIN  = fINmatch - 1.0*deltaf
    fPhaseMatchIM  = fIMmatch + 0.5*deltaf
    fPhaseInsMin  = 0.0026
    fPhaseInsMax  = 1.020 * fMECO
    fPhaseRDMin   = fIMmatch
    fPhaseRDMax   = fring + 1.25*fdamp
    
    #fRef = ifelse(amin(fgrid, axis=0) > fPhaseInsMin, amin(fgrid, axis=0), fPhaseInsMin)
    
    phiNorm = - (3. * (pi^(-5. /3.)))/ 128.
    
    # Ringdown phase collocation points:
    # functionault is to use 5 pseudo-PN coefficients and hence 5 collocation points.
    deltax = fPhaseRDMax - fPhaseRDMin
    xmin   = fPhaseRDMin
    
    CollocationPointsPhaseRD0 = gpoints5[1]* deltax + xmin
    CollocationPointsPhaseRD1 = gpoints5[2]* deltax + xmin
    CollocationPointsPhaseRD2 = gpoints5[3]* deltax + xmin
    # Collocation point 4 is set to the ringdown frequency ~ dip in Lorentzian
    CollocationPointsPhaseRD3 = fring
    CollocationPointsPhaseRD4 = gpoints5[5]* deltax + xmin
    
    CollocationValuesPhaseRD0 = (((eta*(0.7207992174994245 - 1.237332073800276*eta + 6.086871214811216*eta2))/(0.006851189888541745 + 0.06099184229137391*eta - 0.15500218299268662*eta2 + 1. *eta2*eta)) + (((0.06519048552628343 - 25.25397971063995*eta - 308.62513664956975*eta2*eta2 + 58.59408241189781*eta2 + 160.14971486043524*eta2*eta)*totchi + eta*(-5.215945111216946 + 153.95945758807616*eta - 693.0504179144295*eta2 + 835.1725103648205*eta2*eta)*totchi2 + (0.20035146870472367 - 0.28745205203100666*eta - 47.56042058800358*eta2*eta2)*totchi2*totchi + eta*(5.7756520242745735 - 43.97332874253772*eta + 338.7263666984089*eta2*eta)*totchi2*totchi2 + (-0.2697933899920511 + 4.917070939324979*eta - 22.384949087140086*eta2*eta2 - 11.61488280763592*eta2)*totchi2*totchi2*totchi)/(1. - 0.6628745847248266*totchi)) + (-23.504907495268824*dchi*delta*eta2))
    CollocationValuesPhaseRD1 = (((eta*(-9.460253118496386 + 9.429314399633007*eta + 64.69109972468395*eta2))/(-0.0670554310666559 - 0.09987544893382533*eta + 1. *eta2)) + ((17.36495157980372*eta*totchi + eta2*eta*totchi*(930.3458437154668 + 808.457330742532*totchi) + eta2*eta2*totchi*(-774.3633787391745 - 2177.554979351284*totchi - 1031.846477275069*totchi2) + eta2*totchi*(-191.00932194869588 - 62.997389062600035*totchi + 64.42947340363101*totchi2) + 0.04497628581617564*totchi2*totchi)/(1. - 0.7267610313751913*totchi)) + (dchi*delta*(-36.66374091965371 + 91.60477826830407*eta)*eta2))
    CollocationValuesPhaseRD2 = (((eta*(-8.506898502692536 + 13.936621412517798*eta))/(-0.40919671232073945 + 1. *eta)) + ((eta*(1.7280582989361533*totchi + 18.41570325463385*totchi2*totchi - 13.743271480938104*totchi2*totchi2) + eta2*(73.8367329022058*totchi - 95.57802408341716*totchi2*totchi + 215.78111099820157*totchi2*totchi2) + 0.046849371468156265*totchi2 + eta2*eta*totchi*(-27.976989112929353 + 6.404060932334562*totchi - 633.1966645925428*totchi2*totchi + 109.04824706217418*totchi2))/(1. - 0.6862449113932192*totchi)) + (641.8965762829259*dchi*delta*eta2*eta2*eta))
    CollocationValuesPhaseRD3 = (((-85.86062966719405 - 4616.740713893726*eta - 4925.756920247186*eta2 + 7732.064464348168*eta2*eta + 12828.269960300782*eta2*eta2 - 39783.51698102803*eta2*eta2*eta)/(1. + 50.206318806624004*eta)) + ((totchi*(33.335857451144356 - 36.49019206094966*totchi + eta2*eta*(1497.3545918387515 - 101.72731770500685*totchi)*totchi - 3.835967351280833*totchi2 + 2.302712009652155*totchi2*totchi + eta2*(93.64156367505917 - 18.184492163348665*totchi + 423.48863373726243*totchi2 - 104.36120236420928*totchi2*totchi - 719.8775484010988*totchi2*totchi2) + 1.6533417657003922*totchi2*totchi2 + eta*(-69.19412903018717 + 26.580344399838758*totchi - 15.399770764623746*totchi2 + 31.231253209893488*totchi2*totchi + 97.69027029734173*totchi2*totchi2) + eta2*eta2*(1075.8686153198323 - 3443.0233614187396*totchi - 4253.974688619423*totchi2 - 608.2901586790335*totchi2*totchi + 5064.173605639933*totchi2*totchi2)))/(-1.3705601055555852 + 1. *totchi)) + (dchi*delta*eta*(22.363215261437862 + 156.08206945239374*eta)))
    CollocationValuesPhaseRD4 = (((eta*(7.05731400277692 + 22.455288821807095*eta + 119.43820622871043*eta2))/(0.26026709603623255 + 1. *eta)) + ((eta2*(134.88158268621922 - 56.05992404859163*totchi)*totchi + eta*totchi*(-7.9407123129681425 + 9.486783128047414*totchi) + eta2*eta*totchi*(-316.26970506215554 + 90.31815139272628*totchi))/(1. - 0.7162058321905909*totchi)) + (43.82713604567481*dchi*delta*eta2*eta))

    CollocationValuesPhaseRD4 = CollocationValuesPhaseRD4 + CollocationValuesPhaseRD3
    CollocationValuesPhaseRD2 = CollocationValuesPhaseRD2 + CollocationValuesPhaseRD3
    CollocationValuesPhaseRD1 = CollocationValuesPhaseRD1 + CollocationValuesPhaseRD3
    CollocationValuesPhaseRD0 = CollocationValuesPhaseRD0 + CollocationValuesPhaseRD1
    
    A0i = [1., CollocationPointsPhaseRD0^(-1. /3.), CollocationPointsPhaseRD0^(-2), CollocationPointsPhaseRD0^(-4), -(dphase0) / (fdamp*fdamp + (CollocationPointsPhaseRD0 - fring)*(CollocationPointsPhaseRD0 - fring))]
    A1i = [1., CollocationPointsPhaseRD1^(-1. /3.), CollocationPointsPhaseRD1^(-2), CollocationPointsPhaseRD1^(-4), -(dphase0) / (fdamp*fdamp + (CollocationPointsPhaseRD1 - fring)*(CollocationPointsPhaseRD1 - fring))]
    A2i = [1., CollocationPointsPhaseRD2^(-1. /3.), CollocationPointsPhaseRD2^(-2), CollocationPointsPhaseRD2^(-4), -(dphase0) / (fdamp*fdamp + (CollocationPointsPhaseRD2 - fring)*(CollocationPointsPhaseRD2 - fring))]
    A3i = [1., CollocationPointsPhaseRD3^(-1. /3.), CollocationPointsPhaseRD3^(-2), CollocationPointsPhaseRD3^(-4), -(dphase0) / (fdamp*fdamp + (CollocationPointsPhaseRD3 - fring)*(CollocationPointsPhaseRD3 - fring))]
    A4i = [1., CollocationPointsPhaseRD4^(-1. /3.), CollocationPointsPhaseRD4^(-2), CollocationPointsPhaseRD4^(-4), -(dphase0) / (fdamp*fdamp + (CollocationPointsPhaseRD4 - fring)*(CollocationPointsPhaseRD4 - fring))]
    
    Acoloc = permutedims(cat(A0i, A1i, A2i, A3i, A4i, dims=2), (2,1))
    bcoloc = [CollocationValuesPhaseRD0, CollocationValuesPhaseRD1, CollocationValuesPhaseRD2, CollocationValuesPhaseRD3, CollocationValuesPhaseRD4]
    coeffscoloc = Acoloc \ bcoloc # solve linear siystems in julia
    c0coloc  = coeffscoloc[1]
    c1coloc  = coeffscoloc[2]
    c2coloc  = coeffscoloc[3]
    c4coloc  = coeffscoloc[4]
    cRDcoloc = coeffscoloc[5]
    cLcoloc  = -(dphase0 * cRDcoloc)
    phaseRD  = CollocationValuesPhaseRD0
    
    # Inspiral phase collocation points
    # functionault is to use 4 pseudo-PN coefficients and hence 4 collocation points.
    
    deltax      = fPhaseInsMax - fPhaseInsMin
    xmin        = fPhaseInsMin
    
    if InsPhaseVersion == 104
        CollocationPointsPhaseIns0 = gpoints4[1]* deltax + xmin
        CollocationPointsPhaseIns1 = gpoints4[2]* deltax + xmin
        CollocationPointsPhaseIns2 = gpoints4[3]* deltax + xmin
        CollocationPointsPhaseIns3 = gpoints4[4]* deltax + xmin

        CollocationValuesPhaseIns0 = (((-17294.000000000007 - 19943.076428555978*eta + 483033.0998073767*eta2)/(1. + 4.460294035404433*eta)) + ((chiPN*(68384.62786426462 + 67663.42759836042*chiPN - 2179.3505885609297*chiPN2 + eta*(-58475.33302037833 + 62190.404951852535*chiPN + 18298.307770807573*chiPN2 - 303141.1945565486*chiPN2*chiPN) + 19703.894135534803*chiPN2*chiPN + eta2*(-148368.4954044637 - 758386.5685734496*chiPN - 137991.37032619823*chiPN2 + 1.0765877367729193e6*chiPN2*chiPN) + 32614.091002011017*chiPN2*chiPN2))/(2.0412979553629143 + 1. *chiPN)) + (12017.062595934838*dchi*delta*eta))
        CollocationValuesPhaseIns1 = (((-7579.300000000004 - 120297.86185566607*eta + 1.1694356931282217e6*eta2 - 557253.0066989232*eta2*eta)/(1. + 18.53018618227582*eta)) + ((chiPN*(-27089.36915061857 - 66228.9369155027*chiPN + eta2*(150022.21343386435 - 50166.382087278434*chiPN - 399712.22891153296*chiPN2) - 44331.41741405198*chiPN2 + eta*(50644.13475990821 + 157036.45676788126*chiPN + 126736.43159783827*chiPN2) + eta2*eta*(-593633.5370110178 - 325423.99477314285*chiPN + 847483.2999508682*chiPN2)))/(-1.5232497464826662 - 3.062957826830017*chiPN - 1.130185486082531*chiPN2 + 1. *chiPN2*chiPN)) + (3843.083992827935*dchi*delta*eta))
        CollocationValuesPhaseIns2 = (((15415.000000000007 + 873401.6255736464*eta + 376665.64637025696*eta2 - 3.9719980569125614e6*eta2*eta + 8.913612508054944e6*eta2*eta2)/(1. + 46.83697749859996*eta)) + ((chiPN*(397951.95299014193 - 207180.42746987*chiPN + eta2*eta*(4.662143741417853e6 - 584728.050612325*chiPN - 1.6894189124921719e6*chiPN2) + eta*(-1.0053073129700898e6 + 1.235279439281927e6*chiPN - 174952.69161683554*chiPN2) - 130668.37221912303*chiPN2 + eta2*(-1.9826323844247842e6 + 208349.45742548333*chiPN + 895372.155565861*chiPN2)))/(-9.675704197652225 + 3.5804521763363075*chiPN + 2.5298346636273306*chiPN2 + 1. *chiPN2*chiPN)) + (-1296.9289110696955*dchi*dchi*eta + dchi*delta*eta*(-24708.109411857182 + 24703.28267342699*eta + 47752.17032707405*chiPN)))
        CollocationValuesPhaseIns3 = ((2439.000000000001 - 31133.52170083207*eta + 28867.73328134167*eta2)/(1. + 0.41143032589262585*eta)) + ((chiPN*(16116.057657391262 + eta2*eta*(-375818.0132734753 - 386247.80765802023*chiPN) + eta*(-82355.86732027541 - 25843.06175439942*chiPN) + 9861.635308837876*chiPN + eta2*(229284.04542668918 + 117410.37432997991*chiPN)))/(-3.7385208695213668 + 0.25294420589064653*chiPN + 1. *chiPN2)) + (194.5554531509207*dchi*delta*eta)

        CollocationValuesPhaseIns0 = CollocationValuesPhaseIns0 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns1 = CollocationValuesPhaseIns1 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns3 = CollocationValuesPhaseIns3 + CollocationValuesPhaseIns2

        A0i = [1., CollocationPointsPhaseIns0^(1. /3.), CollocationPointsPhaseIns0^(2. /3.), CollocationPointsPhaseIns0]
        A1i = [1., CollocationPointsPhaseIns1^(1. /3.), CollocationPointsPhaseIns1^(2. /3.), CollocationPointsPhaseIns1]
        A2i = [1., CollocationPointsPhaseIns2^(1. /3.), CollocationPointsPhaseIns2^(2. /3.), CollocationPointsPhaseIns2]
        A3i = [1., CollocationPointsPhaseIns3^(1. /3.), CollocationPointsPhaseIns3^(2. /3.), CollocationPointsPhaseIns3]

        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseIns0, CollocationValuesPhaseIns1, CollocationValuesPhaseIns2, CollocationValuesPhaseIns3]

        coeffscoloc = Acoloc \ bcoloc
        a0coloc  = coeffscoloc[1]
        a1coloc  = coeffscoloc[2]
        a2coloc  = coeffscoloc[3]
        a3coloc  = coeffscoloc[4]
        a4coloc  = 0.

    
    elseif InsPhaseVersion == 114
        
        CollocationPointsPhaseIns0 = gpoints4[1]* deltax + xmin
        CollocationPointsPhaseIns1 = gpoints4[2]* deltax + xmin
        CollocationPointsPhaseIns2 = gpoints4[3]* deltax + xmin
        CollocationPointsPhaseIns3 = gpoints4[4]* deltax + xmin
        
        CollocationValuesPhaseIns0 = (((-36664.000000000015 + 277640.10051158903*eta - 581120.4916255298*eta2 + 1.415628418251648e6*eta2*eta - 7.640937162029471e6*eta2*eta2 + 1.1572710625359124e7*eta2*eta2*eta)/(1. - 4.011038704323779*eta)) + ((chiPN*(-38790.01253014577 - 50295.77273512981*chiPN + 15182.324439704937*chiPN2 + eta2*(57814.07222969789 + 344650.11918139807*chiPN + 17020.46497164955*chiPN2 - 574047.1384792664*chiPN2*chiPN) + 24626.598127509922*chiPN2*chiPN + eta*(23058.264859112394 - 16563.935447608965*chiPN - 36698.430436426395*chiPN2 + 105713.91549712936*chiPN2*chiPN)))/(-1.5445637219268247 - 0.24997068896075847*chiPN + 1. *chiPN2)) + (74115.77361380383*dchi*delta*eta2))
        CollocationValuesPhaseIns1 = (((-17762.000000000007 - 1.6929191194109183e6*eta + 8.420903644926643e6*eta2)/(1. + 98.061533474615*eta)) + ((chiPN*(-46901.6486082098 - 83648.57463631754*chiPN + eta2*(1.2502334322912344e6 + 1.4500798116821344e6*chiPN - 1.4822181506831646e6*chiPN2) - 41234.966418619966*chiPN2 + eta*(-24017.33452114588 - 15241.079745314566*chiPN + 136554.48806839858*chiPN2) + eta2*eta*(-3.584298922116994e6 - 3.9566921791790277e6*chiPN + 4.357599992831832e6*chiPN2)))/(-3.190220646817508 - 3.4308485421201387*chiPN - 0.6347932583034377*chiPN2 + 1. *chiPN2*chiPN)) + (24906.33337911219*dchi*delta*eta2))
        CollocationValuesPhaseIns2 = (((68014.00000000003 + 1.1513072539654972e6*eta - 2.725589921577228e6*eta2 + 312571.92531733884*eta2*eta)/(1. + 17.48539665509149*eta)) + ((chiPN*(-34467.00643820664 + 99693.81839115614*eta + 144345.24343461913*eta2*eta2 + (23618.044919850676 - 89494.69555164348*eta + 725554.5749749158*eta2*eta2 - 103449.15865381068*eta2)*chiPN + (10350.863429774612 - 73238.45609787296*eta + 3.559251543095961e6*eta2*eta2 + 888228.5439003729*eta2 - 3.4602940487291473e6*eta2*eta)*chiPN2))/(1. - 0.056846656084188936*chiPN - 0.32681474740130184*chiPN2 - 0.30562055811022015*chiPN2*chiPN)) + (-1182.4036752941936*dchi*dchi*eta + dchi*delta*eta*(-0.39185419821851025 - 99764.21095663306*eta + 41826.177356107364*chiPN)))
        CollocationValuesPhaseIns3 = (((5749.000000000003 - 37877.95816426952*eta)/(1. + 1.1883386102990128*eta)) + (((-4285.982163759047 + 24558.689969419473*eta - 49270.2296311733*eta2)*chiPN + eta*(-24205.71407420114 + 70777.38402634041*eta)*chiPN2 + (2250.661418551257 + 187.95136178643946*eta - 11976.624134935797*eta2)*chiPN2*chiPN)/(1. - 0.7220334077284601*chiPN)) + (dchi*delta*eta*(339.69292150803585 - 3459.894150148715*chiPN)))
        
        CollocationValuesPhaseIns0 = CollocationValuesPhaseIns0 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns1 = CollocationValuesPhaseIns1 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns3 = CollocationValuesPhaseIns3 + CollocationValuesPhaseIns2

        A0i = [1., CollocationPointsPhaseIns0^(1. /3.), CollocationPointsPhaseIns0^(2. /3.), CollocationPointsPhaseIns0]
        A1i = [1., CollocationPointsPhaseIns1^(1. /3.), CollocationPointsPhaseIns1^(2. /3.), CollocationPointsPhaseIns1]
        A2i = [1., CollocationPointsPhaseIns2^(1. /3.), CollocationPointsPhaseIns2^(2. /3.), CollocationPointsPhaseIns2]
        A3i = [1., CollocationPointsPhaseIns3^(1. /3.), CollocationPointsPhaseIns3^(2. /3.), CollocationPointsPhaseIns3]

        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseIns0, CollocationValuesPhaseIns1, CollocationValuesPhaseIns2, CollocationValuesPhaseIns3]

        coeffscoloc = Acoloc \ bcoloc
        a0coloc  = coeffscoloc[1]
        a1coloc  = coeffscoloc[2]
        a2coloc  = coeffscoloc[3]
        a3coloc  = coeffscoloc[4]
        a4coloc  = 0.

    elseif InsPhaseVersion == 105
        
        CollocationPointsPhaseIns0 = gpoints5[1]* deltax + xmin
        CollocationPointsPhaseIns1 = gpoints5[2]* deltax + xmin
        CollocationPointsPhaseIns2 = gpoints5[3]* deltax + xmin
        CollocationPointsPhaseIns3 = gpoints5[4]* deltax + xmin
        CollocationPointsPhaseIns4 = gpoints5[5]* deltax + xmin
        
        CollocationValuesPhaseIns0 = (((-14234.000000000007 + 16956.107542097994*eta + 176345.7518697656*eta2)/(1. + 1.294432443903631*eta)) + ((chiPN*(814.4249470391651 + 539.3944162216742*chiPN + 1985.8695471257474*chiPN2 + eta*(-918.541450687484 + 2531.457116826593*chiPN - 14325.55532310136*chiPN2 - 19213.48002675173*chiPN2*chiPN) + 1517.4561843692861*chiPN2*chiPN + eta2*(-517.7142591907573 - 14328.567448748548*chiPN + 21305.033147575057*chiPN2 + 50143.99945676916*chiPN2*chiPN)))/(0.03332712934306297 + 0.0025905919215826172*chiPN + (0.07388087063636305 - 0.6397891808905215*eta + 1. *eta2)*chiPN2)) + (dchi*delta*eta*(0.09704682517844336 + 69335.84692284222*eta)))
        CollocationValuesPhaseIns1 = (((-7520.900000000003 - 49463.18828584058*eta + 437634.8057596484*eta2)/(1. + 9.10538019868398*eta)) + ((chiPN*(25380.485895523005 + 30617.553968012628*chiPN + 5296.659585425608*chiPN2 + eta*(-49447.74841021405 - 94312.78229903466*chiPN - 5731.614612941746*chiPN2*chiPN) + 2609.50444822972*chiPN2*chiPN + 5206.717656940992*chiPN2*chiPN2 + eta2*(54627.758819129864 + 157460.98527210607*chiPN - 69726.85196686552*chiPN2 + 4674.992397927943*chiPN2*chiPN + 20704.368650323784*chiPN2*chiPN2)))/(1.5668927528319367 + 1. *chiPN)) + (-95.38600275845481*dchi*dchi*eta + dchi*delta*eta*(3271.8128884730654 + 12399.826237672185*eta + 9343.380589951552*chiPN)))
        CollocationValuesPhaseIns2 = (((11717.332402222377 + 4.972361134612872e6*eta + 2.137585030930089e7*eta2 - 8.882868155876668e7*eta2*eta + 2.4104945956043008e8*eta2*eta2 - 2.3143426271719798e8*eta2*eta2*eta)/(1. + 363.5524719849582*eta)) + ((chiPN*(52.42001436116159 - 50.547943589389966*chiPN + eta2*eta*chiPN*(-15355.56020802297 + 20159.588079899433*chiPN) + eta2*(-286.9576245212502 + 2795.982637986682*chiPN - 2633.1870842242447*chiPN2) - 1.0824224105690476*chiPN2 + eta*(-123.78531181532225 + 136.1961976556154*chiPN - 7.534890781638552*chiPN2*chiPN) + 5.973206330664007*chiPN2*chiPN + eta2*eta2*(1777.2176433016125 + 24069.288079063674*chiPN - 44794.9522164669*chiPN2 + 1584.1515277998406*chiPN2*chiPN)))/(-0.0015307616935628491 + (0.0010676159178395538 - 0.25*eta2*eta + 1. *eta2*eta2)*chiPN)) + (-1357.9794908614106*dchi*dchi*eta + dchi*delta*eta*(-23093.829989687543 + 21908.057881789653*eta + 49493.91485992256*chiPN)))
        CollocationValuesPhaseIns3 = (((4085.300000000002 + 62935.7755506329*eta - 1.3712743918777364e6*eta2 + 5.024685134555112e6*eta2*eta - 3.242443755025284e6*eta2*eta2)/(1. + 20.889132970603523*eta - 99.39498823723363*eta2)) + ((chiPN*(-299.6987332025542 - 106.2596940493108*chiPN + eta2*eta*(2383.4907865977148 - 13637.11364447208*chiPN - 14808.138346145908*chiPN2) + eta*(1205.2993091547498 - 508.05838536573464*chiPN - 1453.1997617403304*chiPN2) + 132.22338129554674*chiPN2 + eta2*(-2438.4917103042208 + 5032.341879949591*chiPN + 7026.9206794027405*chiPN2)))/(0.03089183275944264 + 1. *eta2*eta*chiPN - 0.010670764224621489*chiPN2)) + (-1392.6847421907178*dchi*delta*eta))
        CollocationValuesPhaseIns4 = (((5474.400000000003 + 131008.0112992443*eta - 1.9692364337640922e6*eta2 + 1.8732325307375633e6*eta2*eta)/(1. + 32.90929274981482*eta)) + ((chiPN*(18609.016486281424 - 1337.4947536109685*chiPN + eta2*(219014.98908698096 - 307162.33823247004*chiPN - 124111.02067626518*chiPN2) - 7394.9595046977365*chiPN2 + eta*(-87820.37490863055 + 53564.4178831741*chiPN + 34070.909093771494*chiPN2) + eta2*eta*(-248096.84257893753 + 536024.5354098587*chiPN + 243877.65824670633*chiPN2)))/(-1.5282904337787517 + 1. *chiPN)) + (-429.1148607925461*dchi*delta*eta))
        
        CollocationValuesPhaseIns0 = CollocationValuesPhaseIns0 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns1 = CollocationValuesPhaseIns1 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns3 = CollocationValuesPhaseIns3 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns4 = CollocationValuesPhaseIns4 + CollocationValuesPhaseIns2
        
        A0i = [1., CollocationPointsPhaseIns0^(1. /3.), CollocationPointsPhaseIns0^(2. /3.), CollocationPointsPhaseIns0, CollocationPointsPhaseIns0^(4. /3.)]
        A1i = [1., CollocationPointsPhaseIns1^(1. /3.), CollocationPointsPhaseIns1^(2. /3.), CollocationPointsPhaseIns1, CollocationPointsPhaseIns1^(4. /3.)]
        A2i = [1., CollocationPointsPhaseIns2^(1. /3.), CollocationPointsPhaseIns2^(2. /3.), CollocationPointsPhaseIns2, CollocationPointsPhaseIns2^(4. /3.)]
        A3i = [1., CollocationPointsPhaseIns3^(1. /3.), CollocationPointsPhaseIns3^(2. /3.), CollocationPointsPhaseIns3, CollocationPointsPhaseIns3^(4. /3.)]
        A4i = [1., CollocationPointsPhaseIns4^(1. /3.), CollocationPointsPhaseIns4^(2. /3.), CollocationPointsPhaseIns4, CollocationPointsPhaseIns4^(4. /3.)]
        
        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i, A4i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseIns0, CollocationValuesPhaseIns1, CollocationValuesPhaseIns2, CollocationValuesPhaseIns3, CollocationValuesPhaseIns4]
    
        coeffscoloc = Acoloc \ bcoloc
        a0coloc  = coeffscoloc[1]
        a1coloc  = coeffscoloc[2]
        a2coloc  = coeffscoloc[3]
        a3coloc  = coeffscoloc[4]
        a4coloc  = coeffscoloc[5]
        
    elseif InsPhaseVersion == 115
        
        CollocationPointsPhaseIns0 = gpoints5[1]* deltax + xmin
        CollocationPointsPhaseIns1 = gpoints5[2]* deltax + xmin
        CollocationPointsPhaseIns2 = gpoints5[3]* deltax + xmin
        CollocationPointsPhaseIns3 = gpoints5[4]* deltax + xmin
        CollocationPointsPhaseIns4 = gpoints5[5]* deltax + xmin
        
        CollocationValuesPhaseIns0 = (((-29240.00000000001 - 12488.41035199958*eta + 1.3911845288427814e6*eta2 - 3.492477584609041e6*eta2*eta)/(1. + 2.6711462529779824*eta - 26.80876660227278*eta2)) + ((chiPN*(-29536.155624432842 - 40024.5988680615*chiPN + 11596.401177843705*chiPN2 + eta2*(122185.06346551726 + 351091.59147835104*chiPN - 37366.6143666202*chiPN2 - 505834.54206320125*chiPN2*chiPN) + 20386.815769841945*chiPN2*chiPN + eta*(-9638.828456576934 - 30683.453790630676*chiPN - 15694.962417099561*chiPN2 + 91690.51338194775*chiPN2*chiPN)))/(-1.5343852108869265 - 0.2215651087041266*chiPN + 1. *chiPN2)) + (68951.75344813892*dchi*delta*eta2))
        CollocationValuesPhaseIns1 = (((-18482.000000000007 - 1.2846128476247871e6*eta + 4.894853535651343e6*eta2 + 3.1555931338015324e6*eta2*eta)/(1. + 82.79386070797756*eta)) + ((chiPN*(-19977.10130179636 - 24729.114403562427*chiPN + 10325.295899053815*chiPN2 + eta*(30934.123894659646 + 58636.235226102894*chiPN - 32465.70372990005*chiPN2 - 38032.16219587224*chiPN2*chiPN) + 15485.725898689267*chiPN2*chiPN + eta2*(-38410.1127729419 - 87088.84559983511*chiPN + 61286.73536122325*chiPN2 + 42503.913487705235*chiPN2*chiPN)))/(-1.5148031011828222 - 0.24267195338525768*chiPN + 1. *chiPN2)) + (5661.027157084334*dchi*delta*eta))
        CollocationValuesPhaseIns2 = (((60484.00000000003 + 4.370611564781374e6*eta - 5.785128542827255e6*eta2 - 8.82141241633613e6*eta2*eta + 1.3406996319926713e7*eta2*eta2)/(1. + 70.89393713617065*eta)) + ((chiPN*(21.91241092620993 - 32.57779678272056*chiPN + eta2*(-102.4982890239095 + 2570.7566494633033*chiPN - 2915.1250015652076*chiPN2) + 8.130585173873232*chiPN2 + eta*(-28.158436727309365 + 47.42665852824481*chiPN2) + eta2*eta*(-1635.6270690726785 - 13745.336370568011*chiPN + 19539.310912464192*chiPN2) + 1.2792283911312285*chiPN2*chiPN + eta2*eta2*(5558.382039622131 + 21398.7730201213*chiPN - 37773.40511355719*chiPN2 + 768.6183344184254*chiPN2*chiPN)))/(-0.0007758753818017038 + (0.0005304023864415552 - 0.25000000000000006*eta2*eta + 1. *eta2*eta2)*chiPN)) + (-1223.769262298142*dchi*dchi*eta + dchi*delta*eta*(-16.705471562129436 - 93771.93750060834*eta + 43675.70151058481*chiPN)))
        CollocationValuesPhaseIns3 = (((9760.400000000005 + 9839.852773121198*eta - 398521.0434645335*eta2 + 267575.4709475981*eta2*eta)/(1. + 6.1355249449135005*eta)) + ((chiPN*(-1271.406488219572 + eta2*(-9641.611385554736 - 9620.333878140807*chiPN) - 1144.516395171019*chiPN + eta*(5155.337817255137 + 4450.755534615418*chiPN)))/(0.1491519640750958 + (-0.0008208549820159909 - 0.15468508831447628*eta + 0.7266887643762937*eta2)*chiPN + (0.02282542856845755 - 0.445924460572114*eta + 1. *eta2)*chiPN2)) + (-1366.7949288045616*dchi*delta*eta))
        CollocationValuesPhaseIns4 = (((12971.000000000005 - 93606.05144508784*eta + 102472.4473167639*eta2)/(1. - 0.8909484992212859*eta)) + ((chiPN*(16182.268123259992 + 3513.8535400032874*chiPN + eta2*(343388.99445324624 - 240407.0282222587*chiPN - 312202.59917289804*chiPN2) - 10814.056847109632*chiPN2 + eta*(-94090.9232151429 + 35305.66247590705*chiPN + 65450.36389642103*chiPN2) + eta2*eta*(-484443.15601144277 + 449511.3965208116*chiPN + 552355.592066788*chiPN2)))/(-1.4720837917195788 + 1. *chiPN)) + (-494.2754225110706*dchi*delta*eta))
        
        CollocationValuesPhaseIns0 = CollocationValuesPhaseIns0 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns1 = CollocationValuesPhaseIns1 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns3 = CollocationValuesPhaseIns3 + CollocationValuesPhaseIns2
        CollocationValuesPhaseIns4 = CollocationValuesPhaseIns4 + CollocationValuesPhaseIns2
        
        A0i = [1., CollocationPointsPhaseIns0^(1. /3.), CollocationPointsPhaseIns0^(2. /3.), CollocationPointsPhaseIns0, CollocationPointsPhaseIns0^(4. /3.)]
        A1i = [1., CollocationPointsPhaseIns1^(1. /3.), CollocationPointsPhaseIns1^(2. /3.), CollocationPointsPhaseIns1, CollocationPointsPhaseIns1^(4. /3.)]
        A2i = [1., CollocationPointsPhaseIns2^(1. /3.), CollocationPointsPhaseIns2^(2. /3.), CollocationPointsPhaseIns2, CollocationPointsPhaseIns2^(4. /3.)]
        A3i = [1., CollocationPointsPhaseIns3^(1. /3.), CollocationPointsPhaseIns3^(2. /3.), CollocationPointsPhaseIns3, CollocationPointsPhaseIns3^(4. /3.)]
        A4i = [1., CollocationPointsPhaseIns4^(1. /3.), CollocationPointsPhaseIns4^(2. /3.), CollocationPointsPhaseIns4, CollocationPointsPhaseIns4^(4. /3.)]
        
        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i, A4i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseIns0, CollocationValuesPhaseIns1, CollocationValuesPhaseIns2, CollocationValuesPhaseIns3, CollocationValuesPhaseIns4]
    
        coeffscoloc = Acoloc \ bcoloc
        a0coloc  = coeffscoloc[1]
        a1coloc  = coeffscoloc[2]
        a2coloc  = coeffscoloc[3]
        a3coloc  = coeffscoloc[4]
        a4coloc  = coeffscoloc[5]
    
    else
        error("Inspiral phase version not implemented.")
    end
        
        
    sigma1 = (-5. /3.) * a0coloc
    sigma2 = (-5. /4.) * a1coloc
    sigma3 = (-5. /5.) * a2coloc
    sigma4 = (-5. /6.) * a3coloc
    sigma5 = (-5. /7.) * a4coloc
    
    # TaylorF2 PN Coefficients
    # Newtonian
    phi0   = 1.
    # .5 PN
    phi1   = 0.
    # 1 PN
    phi2   = (3715. /756. + (55. *eta)/9.) * (pi^(2. /3.))
    # 1.5 PN, Non-Spinning and Spin-Orbit
    phi3   = -16.0*pi*pi + ((113. *(chi1 + chi2 + chi1*delta - chi2*delta) - 76. *(chi1 + chi2)*eta)/6.) * pi
    # 2 PN, Non-Spinning and Spin-Spin
    phi4   = (15293365. /508032. + (27145. *eta)/504. + (3085*eta2)/72.)*(pi^(4. /3.)) + ((-5. *(81. *chi12*(1. + delta - 2. *eta) + 316. *chi1dotchi2*eta - 81. *chi22*(-1. + delta + 2. *eta)))/16.)*(pi^(4. /3.))
    # 2.5 PN, Non-Spinning and Spin-Orbit
    phi5L  = (((5. *(46374. - 6552. *eta)*pi)/4536.) + ((-732985. *(chi1 + chi2 + chi1*delta - chi2*delta) - 560. *(-1213. *(chi1 + chi2) + 63. *(chi1 - chi2)*delta)*eta + 85680. *(chi1 + chi2)*eta2)/4536.)) * (pi^(5. /3.))
    phi5   = 0.
    # 3 PN, Non-Spinning, Spin-Orbit, Spin-Spin
    phi6   = (11583231236531. /4.69421568e9 - (5. *eta*(3147553127. + 588. *eta*(-45633. + 102260. *eta)))/3.048192e6 - (6848. *MathConstants.eulergamma)/21. -(640. *pi*pi)/3. + (2255. *eta*pi*pi)/12. - (13696*log(2.))/21. - (6848. *log(pi))/63.)*pi*pi + ((5. *(227. *(chi1 + chi2 + chi1*delta - chi2*delta) - 156. *(chi1 + chi2)*eta)*pi)/3.)*pi*pi + ((5. *(20. *chi1dotchi2*eta*(11763. + 12488. *eta) + 7. *chi22*(-15103. *(-1. + delta) + 2. *(-21683. + 6580. *delta)*eta - 9808. *eta2) - 7. *chi12*(-15103. *(1. + delta) + 2. *(21683. + 6580. *delta)*eta + 9808. *eta2)))/4032.)*pi*pi
    # 3 PN log, Non-Spinning, Spin-Orbit and Spin-Spin
    phi6L  = (-6848/63.) * pi*pi
    # 3.5 PN, Non-Spinning, Spin-Orbit, Spin-Spin and Cubic-in-Spin
    phi7   = ((5. *(15419335. + 168. *(75703. - 29618. *eta)*eta)*pi)/254016.) * (pi^(7. /3.)) + ((5. *(-5030016755. *(chi1 + chi2 + chi1*delta - chi2*delta) + 4. *(2113331119. *(chi1 + chi2) + 675484362. *(chi1 - chi2)*delta)*eta - 1008. *(208433. *(chi1 + chi2) + 25011. *(chi1 - chi2)*delta)*eta2 + 90514368. *(chi1 + chi2)*eta2*eta))/6.096384e6) * (pi^(7. /3.)) + (-5. *(57. *chi12*(1. + delta - 2. *eta) + 220. *chi1dotchi2*eta - 57. *chi22*(-1. + delta + 2. *eta))*pi)*(pi^(7. /3.)) + ((14585. *(-(chi22*chi2*(-1. + delta)) + chi12*chi1*(1. + delta)) - 5. *(chi22*chi2*(8819. - 2985. *delta) + 8439. *chi1*chi22*(-1. + delta) - 8439. *chi12*chi2*(1. + delta) + chi12*chi1*(8819. + 2985. *delta))*eta + 40. *(chi1 + chi2)*(17. *chi12 - 14. *chi1dotchi2 + 17. *chi22)*eta2)/48.) * (pi^(7. /3.))
    # 4 PN, Non-Spinning and Spin-Orbit
    phi8   = ((-5. *(1263141. *(chi1 + chi2 + chi1*delta - chi2*delta) - 2. *(794075. *(chi1 + chi2) + 178533. *(chi1 - chi2)*delta)*eta + 94344. *(chi1 + chi2)*eta2)*pi*(-1. + log(pi)))/9072. ) * (pi^(8. /3.))
    # 4 PN log, Non-Spinning and Spin-Orbit
    phi8L  = ((-5. *(1263141. *(chi1 + chi2 + chi1*delta - chi2*delta) - 2. *(794075. *(chi1 + chi2) + 178533. *(chi1 - chi2)*delta)*eta + 94344. *(chi1 + chi2)*eta2)*pi)/9072.) * (pi^(8. /3.))
    
    if InsPhaseVersion in [114, 115]
        # 3.5PN, Leading Order Spin-Spin Tail Term
        phi7   = phi7 + ((5. *(65. *chi12*(1. + delta - 2. *eta) + 252. *chi1dotchi2*eta - 65. *chi22*(-1. + delta + 2. *eta))*pi)/4.) * (pi^(7. /3.))
        # 4.5PN, Tail Term
        phi9   = ((5. *(-256. + 451. *eta)*pi*pi*pi)/6. + (pi*(105344279473163. + 700. *eta*(-298583452147. + 96. *eta*(99645337. + 14453257. *eta)) - 12246091038720. *MathConstants.eulergamma - 24492182077440*log(2.)))/1.877686272e10 - (13696*pi*log(pi))/63. ) * pi*pi*pi
        # 4.5PN, Log Term
        phi9L  = ((-13696*pi)/63.) * pi*pi*pi
    else
        phi9   = 0.
        phi9L  = 0.
    end

    # Normalized Phase Derivatives
    dphi0  = 1.0
    dphi1  = 0.0
    dphi2  = (743. /252. + (11. *eta)/3. )*(pi^(2. /3.))
    dphi3  = ((chi2*(113. - 113. *delta - 76. *eta) + chi1*(113. *(1. + delta) - 76. *eta) - 96. *pi)/15.) * pi
    dphi4  = (3058673/508032. - (81*chi12*(1 + delta))/16. - (79*chi1dotchi2*eta)/4. + (81*chi22*(-1 + delta + 2*eta))/16. + (eta*(5429 + 5103*chi12 + 4319*eta))/504. ) * (pi^(4. /3.))
    dphi5  = ( (-146597*chi2*delta + 146597*(chi1 + chi2 + chi1*delta) + 112*(chi1*(-1213 + 63*delta) - chi2*(1213 + 63*delta))*eta - 17136*(chi1 + chi2)*eta2 + 6*(-7729 + 1092*eta)*pi)/1512. ) * (pi^(5. /3.))
    dphi6  = ( (-10052469856691 + 24236159077900*eta)/2.34710784e10 + (6848*MathConstants.eulergamma)/105. + (-951489*chi12*(1 + delta) - 180*chi1dotchi2*eta*(11763 + 12488*eta) + 63*chi22*(15103*(-1 + delta) + 2*(21683 - 6580*delta)*eta + 9808*eta2) + 7*eta*(18*chi12*(21683 + 6580*delta + 4904*eta) + eta*(-45633 + 102260*eta)) - 12096*(chi2*(227 - 227*delta - 156*eta) + chi1*(227*(1 + delta) - 156*eta))*pi - 3024*(-512 + 451*eta)*pi*pi )/36288. + (13696*log(2.))/105. + (6848*log(pi))/315.0   ) * pi*pi
    dphi6L = (  6848 / 315. ) * pi*pi
    dphi7  = (-(chi12*chi2*eta*(2813*(1 + delta) + 8*eta))/8. + (chi12*chi1*(-2917*(1 + delta) + (8819 + 2985*delta)*eta - 136*eta2))/24. + (chi2*(-5030016755*(-1 + delta) + 4*(-2113331119 + 675484362*delta)*eta + 1008*(208433 - 25011*delta)*eta2 - 90514368*eta2*eta))/3.048192e6 - (chi22*chi2*(2917 + eta*(-8819 + 136*eta) + delta*(-2917 + 2985*eta)))/24. + (chi1*(5030016755 - 8453324476*eta + 1008*eta*((208433 - 89796*eta)*eta - 378*chi22*(2813 + 8*eta)) + delta*(5030016755 + 504*eta*(-5360987 + 2126628*chi22 + 50022*eta))))/3.048192e6 + (-15419335/127008. + 114*chi12*(1 + delta - 2*eta) - (75703*eta)/756. + 440*chi1dotchi2*eta + (14809*eta2)/378. - 114*chi22*(-1 + delta + 2*eta))*pi ) * (pi^(7. /3.))
    dphi8  = ( ((chi1*(1263141 - 1588150*eta + 94344*eta2 - 9*delta*(-140349 + 39674*eta)) + chi2*(1263141 - 1588150*eta + 94344*eta2 + 9*delta*(-140349 + 39674*eta)))*pi)/3024. ) * (pi^(8. /3.))* log(pi)
    dphi8L  = ( ((chi1*(1263141 - 1588150*eta + 94344*eta2 - 9*delta*(-140349 + 39674*eta)) + chi2*(1263141 - 1588150*eta + 94344*eta2 + 9*delta*(-140349 + 39674*eta)))*pi)/3024. ) * (pi^(8. /3.))


    
    if InsPhaseVersion in [114, 115]
        # 3.5PN, Leading Order Spin-Spin Tail Term
        dphi7   = dphi7 + (((-65*chi12*(1 + delta - 2*eta) - 252*chi1dotchi2*eta + 65*chi22*(-1 + delta + 2*eta))* pi)/2. ) * (pi^(7. /3.))
        # 4.5PN, Tail Term
        dphi9 = ( (512/3. - (902*eta)/3.)*pi*pi*pi + pi * (-102282756713483/2.34710784e10 + (298583452147*eta)/3.3530112e7 - (9058667*eta2)/31752. - (2064751*eta2*eta)/49896. + (54784*MathConstants.eulergamma)/105. + (109568*log(2.))/105. + (54784*log(pi))/315.) ) * pi*pi*pi
        # 4.5PN, Log Term
        dphi9L = ( (54784*pi)/315. ) * pi*pi*pi
    else
        dphi9   = 0.
        dphi9L  = 0.
    end

    # Calculate phase at fmatchIN
    # First the standard TaylorF2 part
    phaseIN = dphi0 + dphi1*(fPhaseMatchIN^(1. /3.)) + dphi2*(fPhaseMatchIN^(2. /3.)) + dphi3*fPhaseMatchIN + dphi4*(fPhaseMatchIN^(4. /3.)) + dphi5*(fPhaseMatchIN^(5. /3.)) + (dphi6 + dphi6L*log(fPhaseMatchIN))*fPhaseMatchIN*fPhaseMatchIN + dphi7*(fPhaseMatchIN^(7. /3.)) + (dphi8 + dphi8L*log(fPhaseMatchIN))*(fPhaseMatchIN^(8. /3.)) + (dphi9 + dphi9L*log(fPhaseMatchIN))*fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN
    # Then the pseudo-PN
    phaseIN = phaseIN + a0coloc*(fPhaseMatchIN^(8. /3.)) + a1coloc*fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN + a2coloc*(fPhaseMatchIN^(10. /3.)) + a3coloc*(fPhaseMatchIN^(11. /3.)) + a4coloc*(fPhaseMatchIN^4)
    # Finally the overall phase
    phaseIN = phaseIN*(fPhaseMatchIN^(-8. /3.))*dphase0
    
    # Intermediate phase collocation points:
    # functionault is to use 5 collocation points
    
    deltax      = fPhaseMatchIM - fPhaseMatchIN
    xmin        = fPhaseMatchIN
    
    if IntPhaseVersion == 104
        CollocationPointsPhaseInt0 = gpoints4[1]* deltax + xmin
        CollocationPointsPhaseInt1 = gpoints4[2]* deltax + xmin
        CollocationPointsPhaseInt2 = gpoints4[3]* deltax + xmin
        CollocationPointsPhaseInt3 = gpoints4[4]* deltax + xmin
        
        CollocationValuesPhaseInt0 = phaseIN
        CollocationValuesPhaseInt3 = phaseRD
        
        v2IMmRDv4 = (((eta*(-8.244230124407343 - 182.80239160435949*eta + 638.2046409916306*eta2 - 578.878727101827*eta2*eta))/(-0.004863669418916522 - 0.5095088831841608*eta + 1. *eta2)) + ((totchi*(0.1344136125169328 + 0.0751872427147183*totchi + eta2*(7.158823192173721 + 25.489598292111104*totchi - 7.982299108059922*totchi2) + eta*(-5.792368563647471 + 1.0190164430971749*totchi + 0.29150399620268874*totchi2) + 0.033627267594199865*totchi2 + eta2*eta*(17.426503081351534 - 90.69790378048448*totchi + 20.080325133405847*totchi2)))/(0.03449545664201546 - 0.027941977370442107*totchi + (0.005274757661661763 + 0.0455523144123269*eta - 0.3880379911692037*eta2 + 1. *eta2*eta)*totchi2)) + (160.2975913661124*dchi*delta*eta2))
        v3IMmRDv4 = (((0.3145740304678042 + 299.79825045000655*eta - 605.8886581267144*eta2 - 1302.0954503758007*eta2*eta)/(1. + 2.3903703875353255*eta - 19.03836730923657*eta2)) + ((totchi*(1.150925084789774 - 0.3072381261655531*totchi + eta2*eta2*(12160.22339193134 - 1459.725263347619*totchi - 9581.694749116636*totchi2) + eta2*(1240.3900459406875 - 289.48392062629966*totchi - 1218.1320231846412*totchi2) - 1.6191217310064605*totchi2 + eta*(-41.38762957457647 + 60.47477582076457*totchi2) + eta2*eta*(-7540.259952257055 + 1379.3429194626635*totchi + 6372.99271204178*totchi2)))/(-1.4325421225106187 + 1. *totchi)) + (dchi*delta*eta2*eta*(-444.797248674011 + 1448.47758082697*eta + 152.49290092435044*totchi)))
        v2IM      = (((-84.09400000000004 - 1782.8025405571802*eta + 5384.38721936653*eta2)/(1. + 28.515617312596103*eta + 12.404097877099353*eta2)) + ((totchi*(22.5665046165141 - 39.94907120140026*totchi + 4.668251961072*totchi2 + 12.648584361431245*totchi2*totchi + eta2*(-298.7528127869681 + 14.228745354543983*totchi + 398.1736953382448*totchi2 + 506.94924905801673*totchi2*totchi - 626.3693674479357*totchi2*totchi2) - 5.360554789680035*totchi2*totchi2 + eta*(152.07900889608595 - 121.70617732909962*totchi2 - 169.36311036967322*totchi2*totchi + 182.40064586992762*totchi2*totchi2)))/(-1.1571201220629341 + 1. *totchi)) + (dchi*delta*eta2*eta*(5357.695021063607 - 15642.019339339662*eta + 674.8698102366333*totchi)))
        
        CollocationValuesPhaseInt1 = 0.75*(v2IMmRDv4 + CollocationValuesPhaseRD3) + 0.25*v2IM
        CollocationValuesPhaseInt2 = v3IMmRDv4 + CollocationValuesPhaseRD3
        
        A0i = [1., fring/CollocationPointsPhaseInt0, (fring/CollocationPointsPhaseInt0)*(fring/CollocationPointsPhaseInt0), (fring/CollocationPointsPhaseInt0)^3]
        A1i = [1., fring/CollocationPointsPhaseInt1, (fring/CollocationPointsPhaseInt1)*(fring/CollocationPointsPhaseInt1), (fring/CollocationPointsPhaseInt1)^3]
        A2i = [1., fring/CollocationPointsPhaseInt2, (fring/CollocationPointsPhaseInt2)*(fring/CollocationPointsPhaseInt2), (fring/CollocationPointsPhaseInt2)^3]
        A3i = [1., fring/CollocationPointsPhaseInt3, (fring/CollocationPointsPhaseInt3)*(fring/CollocationPointsPhaseInt3), (fring/CollocationPointsPhaseInt3)^3]
        
        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseInt0 - ((4. * cLcoloc) / ((2. *fdamp)*(2. *fdamp) + (CollocationPointsPhaseInt0 - fring)*(CollocationPointsPhaseInt0 - fring))), CollocationValuesPhaseInt1 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt1 - fring)*(CollocationPointsPhaseInt1 - fring))), CollocationValuesPhaseInt2 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt2 - fring)*(CollocationPointsPhaseInt2 - fring))), CollocationValuesPhaseInt3 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt3 - fring)*(CollocationPointsPhaseInt3 - fring)))]
        
        coeffscoloc = Acoloc \ bcoloc

        b0coloc = coeffscoloc[1]
        b1coloc = coeffscoloc[2] * fring
        b2coloc = coeffscoloc[3] * fring * fring
        b3coloc = 0.
        b4coloc = coeffscoloc[4] * fring * fring * fring * fring
        
        
    elseif IntPhaseVersion == 105
        CollocationPointsPhaseInt0 = gpoints5[1]* deltax + xmin
        CollocationPointsPhaseInt1 = gpoints5[2]* deltax + xmin
        CollocationPointsPhaseInt2 = gpoints5[3]* deltax + xmin
        CollocationPointsPhaseInt3 = gpoints5[4]* deltax + xmin
        CollocationPointsPhaseInt4 = gpoints5[5]* deltax + xmin

        CollocationValuesPhaseInt0 = phaseIN
        CollocationValuesPhaseInt4 = phaseRD

        v2IMmRDv4 = (((eta*(0.9951733419499662 + 101.21991715215253*eta + 632.4731389009143*eta2))/(0.00016803066316882238 + 0.11412314719189287*eta + 1.8413983770369362*eta2 + 1. *eta2*eta)) + ((totchi*(18.694178521101332 + 16.89845522539974*totchi + 4941.31613710257*eta2*totchi + eta*(-697.6773920613674 - 147.53381808989846*totchi2) + 0.3612417066833153*totchi2 + eta2*eta*(3531.552143264721 - 14302.70838220423*totchi + 178.85850322465944*totchi2)))/(2.965640445745779 - 2.7706595614504725*totchi + 1. *totchi2)) + (dchi*delta*eta2*(356.74395864902294 + 1693.326644293169*eta2*totchi)))
        v3IMmRDv4 = (((eta*(-5.126358906504587 - 227.46830225846668*eta + 688.3609087244353*eta2 - 751.4184178636324*eta2*eta))/(-0.004551938711031158 - 0.7811680872741462*eta + 1. *eta2)) + ((totchi*(0.1549280856660919 - 0.9539250460041732*totchi - 539.4071941841604*eta2*totchi + eta*(73.79645135116367 - 8.13494176717772*totchi2) - 2.84311102369862*totchi2 + eta2*eta*(-936.3740515136005 + 1862.9097047992134*totchi + 224.77581754671272*totchi2)))/(-1.5308507364054487 + 1. *totchi)) + (2993.3598520496153*dchi*delta*eta2*eta2*eta2))
        v2IM      = (((-82.54500000000004 - 5.58197349185435e6*eta - 3.5225742421184325e8*eta2 + 1.4667258334378073e9*eta2*eta)/(1. + 66757.12830903867*eta + 5.385164380400193e6*eta2 + 2.5176585751772933e6*eta2*eta)) + ((totchi*(19.416719811164853 - 36.066611959079935*totchi - 0.8612656616290079*totchi2 + eta2*(170.97203068800542 - 107.41099349364234*totchi - 647.8103976942541*totchi2*totchi) + 5.95010003393006*totchi2*totchi + eta2*eta*(-1365.1499998427248 + 1152.425940764218*totchi + 415.7134909564443*totchi2 + 1897.5444343138167*totchi2*totchi - 866.283566780576*totchi2*totchi2) + 4.984750041013893*totchi2*totchi2 + eta*(207.69898051583655 - 132.88417400679026*totchi - 17.671713040498304*totchi2 + 29.071788188638315*totchi2*totchi + 37.462217031512786*totchi2*totchi2)))/(-1.1492259468169692 + 1. *totchi)) + (dchi*delta*eta2*eta*(7343.130973149263 - 20486.813161100774*eta + 515.9898508588834*totchi)))

        CollocationValuesPhaseInt1 = 0.75*(v2IMmRDv4 + CollocationValuesPhaseRD3) + 0.25*v2IM
        CollocationValuesPhaseInt2 = v3IMmRDv4 + CollocationValuesPhaseRD3
        CollocationValuesPhaseInt3 = (((0.4248820426833804 - 906.746595921514*eta - 282820.39946006844*eta2 - 967049.2793750163*eta2*eta + 670077.5414916876*eta2*eta2)/(1. + 1670.9440812294847*eta + 19783.077247023448*eta2)) + ((totchi*(0.22814271667259703 + 1.1366593671801855*totchi + eta2*eta*(3499.432393555856 - 877.8811492839261*totchi - 4974.189172654984*totchi2) + eta*(12.840649528989287 - 61.17248283184154*totchi2) + 0.4818323187946999*totchi2 + eta2*(-711.8532052499075 + 269.9234918621958*totchi + 941.6974723887743*totchi2) + eta2*eta2*(-4939.642457025497 - 227.7672020783411*totchi + 8745.201037897836*totchi2)))/(-1.2442293719740283 + 1. *totchi)) + (dchi*delta*(-514.8494071830514 + 1493.3851099678195*eta)*eta2*eta))
        CollocationValuesPhaseInt3 = CollocationValuesPhaseInt3 + CollocationValuesPhaseInt2

        A0i = [1., fring/CollocationPointsPhaseInt0, (fring/CollocationPointsPhaseInt0)*(fring/CollocationPointsPhaseInt0), (fring/CollocationPointsPhaseInt0)^3, (fring/CollocationPointsPhaseInt0)^4]
        A1i = [1., fring/CollocationPointsPhaseInt1, (fring/CollocationPointsPhaseInt1)*(fring/CollocationPointsPhaseInt1), (fring/CollocationPointsPhaseInt1)^3, (fring/CollocationPointsPhaseInt1)^4]
        A2i = [1., fring/CollocationPointsPhaseInt2, (fring/CollocationPointsPhaseInt2)*(fring/CollocationPointsPhaseInt2), (fring/CollocationPointsPhaseInt2)^3, (fring/CollocationPointsPhaseInt2)^4]
        A3i = [1., fring/CollocationPointsPhaseInt3, (fring/CollocationPointsPhaseInt3)*(fring/CollocationPointsPhaseInt3), (fring/CollocationPointsPhaseInt3)^3, (fring/CollocationPointsPhaseInt3)^4]
        A4i = [1., fring/CollocationPointsPhaseInt4, (fring/CollocationPointsPhaseInt4)*(fring/CollocationPointsPhaseInt4), (fring/CollocationPointsPhaseInt4)^3, (fring/CollocationPointsPhaseInt4)^4]

        Acoloc = permutedims(cat(A0i, A1i, A2i, A3i, A4i; dims=2), (2,1))
        bcoloc = [CollocationValuesPhaseInt0 - ((4. * cLcoloc) / ((2. *fdamp)*(2. *fdamp) + (CollocationPointsPhaseInt0 - fring)*(CollocationPointsPhaseInt0 - fring))), CollocationValuesPhaseInt1 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt1 - fring)*(CollocationPointsPhaseInt1 - fring))), CollocationValuesPhaseInt2 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt2 - fring)*(CollocationPointsPhaseInt2 - fring))), CollocationValuesPhaseInt3 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt3 - fring)*(CollocationPointsPhaseInt3 - fring))), CollocationValuesPhaseInt4 - ((4. * cLcoloc) / ((2. *fdamp)*(2.0*fdamp) + (CollocationPointsPhaseInt4 - fring)*(CollocationPointsPhaseInt4 - fring)))]

        coeffscoloc = Acoloc \ bcoloc

        b0coloc = coeffscoloc[1]
        b1coloc = coeffscoloc[2] * fring
        b2coloc = coeffscoloc[3] * fring * fring
        b3coloc = coeffscoloc[4] * fring * fring * fring
        b4coloc = coeffscoloc[5] * fring * fring * fring * fring
    
    else
        error("Intermediate phase version not implemented.")
    end

    c4ov3 = c4coloc / 3.0
    cLovfda = cLcoloc / fdamp

    # Intermediate phase derivative at fPhaseMatchIN
    DPhiInt =
        b0coloc +
        b1coloc / fPhaseMatchIN +
        b2coloc / (fPhaseMatchIN * fPhaseMatchIN) +
        b3coloc / (fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN) +
        b4coloc / (fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN) +
        (4.0 * cLcoloc) /
        ((4.0 * fdamp * fdamp) + (fPhaseMatchIN - fring) * (fPhaseMatchIN - fring))

    C2Int = phaseIN - DPhiInt
    # Inspiral phase at fPhaseMatchIN
    phiIN =
        phi0 +
        phi1 * (fPhaseMatchIN^(1.0 / 3.0)) +
        phi2 * (fPhaseMatchIN^(2.0 / 3.0)) +
        phi3 * fPhaseMatchIN +
        phi4 * (fPhaseMatchIN^(4.0 / 3.0)) +
        (phi5 + phi5L * log(fPhaseMatchIN)) * (fPhaseMatchIN^(5.0 / 3.0)) +
        (phi6 + phi6L * log(fPhaseMatchIN)) * fPhaseMatchIN * fPhaseMatchIN +
        phi7 * (fPhaseMatchIN^(7.0 / 3.0)) +
        (phi8 + phi8L * log(fPhaseMatchIN)) * (fPhaseMatchIN^(8.0 / 3.0)) +
        (phi9 + phi9L * log(fPhaseMatchIN)) * fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN
    # add the pseudo-PN Coefficients
    phiIN =
        phiIN +
        sigma1 * (fPhaseMatchIN^(8.0 / 3.0)) +
        sigma2 * (fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN) +
        sigma3 * (fPhaseMatchIN^(10.0 / 3.0)) +
        sigma4 * (fPhaseMatchIN^(11.0 / 3.0)) +
        sigma5 * (fPhaseMatchIN^4)
    # finally the normalisation
    phiIN = phiIN * phiNorm * (fPhaseMatchIN^(-5.0 / 3.0))

    # Intermediate phase at fPhaseMatchIN

    phiIM =
        b0coloc * fPhaseMatchIN + b1coloc * log(fPhaseMatchIN) - b2coloc / fPhaseMatchIN -
        b3coloc / (fPhaseMatchIN * fPhaseMatchIN) / 2.0 -
        (b4coloc / (fPhaseMatchIN * fPhaseMatchIN * fPhaseMatchIN) / 3.0) +
        (2.0 * cLcoloc * atan((fPhaseMatchIN - fring) / (2.0 * fdamp))) / fdamp

    C1Int = phiIN - phiIM - (C2Int * fPhaseMatchIN)

    # Intermediate phase at fPhaseMatchIM + connections
    phiIMC =
        b0coloc * fPhaseMatchIM + b1coloc * log(fPhaseMatchIM) - b2coloc / fPhaseMatchIM -
        b3coloc / (fPhaseMatchIM * fPhaseMatchIM) / 2.0 -
        (b4coloc / (fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM) / 3.0) +
        (2.0 * cLcoloc * atan((fPhaseMatchIM - fring) / (2.0 * fdamp))) / fdamp +
        C1Int +
        C2Int * fPhaseMatchIM
    # Ringdown phase at fPhaseMatchIM
    phiRD =
        c0coloc * fPhaseMatchIM + 1.5 * c1coloc * (fPhaseMatchIM^(2.0 / 3.0)) -
        c2coloc / fPhaseMatchIM - c4ov3 / (fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM) +
        (cLovfda * atan((fPhaseMatchIM - fring) / fdamp))
    # Intermediate phase derivative at fPhaseMatchIM + connection
    DPhiIntC =
        b0coloc +
        b1coloc / fPhaseMatchIM +
        b2coloc / (fPhaseMatchIM * fPhaseMatchIM) +
        b3coloc / (fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM) +
        b4coloc / (fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM) +
        (4.0 * cLcoloc) /
        ((4.0 * fdamp * fdamp) + (fPhaseMatchIM - fring) * (fPhaseMatchIM - fring)) +
        C2Int
    # Ringdown phase derivative at fPhaseMatchIM
    DPhiRD = (
        c0coloc +
        c1coloc * (fPhaseMatchIM^(-1.0 / 3.0)) +
        c2coloc / (fPhaseMatchIM * fPhaseMatchIM) +
        c4coloc / (fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM * fPhaseMatchIM) +
        (cLcoloc / (fdamp * fdamp + (fPhaseMatchIM - fring) * (fPhaseMatchIM - fring)))
    )

    C2MRD = DPhiIntC - DPhiRD
    C1MRD = phiIMC - phiRD - C2MRD*fPhaseMatchIM



    Phase22 = Phase22Struct(
        fPhaseMatchIN,
        fPhaseMatchIM,
        phi0,
        phi1,
        phi2,
        phi3,
        phi4,
        phi5,
        phi5L,
        phi6,
        phi6L,
        phi7,
        phi8,
        phi8L,
        phi9,
        phi9L,
        dphase0,
        dphi0,
        dphi1,
        dphi2,
        dphi3,
        dphi4,
        dphi5,
        dphi6,
        dphi6L,
        dphi7,
        dphi8,
        dphi8L,
        dphi9,
        dphi9L,
        a0coloc,
        a1coloc,
        a2coloc,
        a3coloc,
        a4coloc,
        b0coloc,
        b1coloc,
        b2coloc,
        b3coloc,
        b4coloc,
        cLcoloc,
        C2Int,
        c0coloc,
        c1coloc,
        c2coloc,
        c4coloc,
        C2MRD,
        sigma1,
        sigma2,
        sigma3,
        sigma4,
        sigma5,
        C1Int,
        C1MRD,
    )
    
    return Phase22

end

# """
# Compute the phase of the GW as a function of frequency, given the connection coefficients. Used only in PhenomXHM.
# """
# function Phi(PhenomXAS(), f, eta, Phase22::Phase22Structure)

#     # Unpack the connection coefficients
#     fPhaseMatchIN = Phase22.fPhaseMatchIN
#     fPhaseMatchIM = Phase22.fPhaseMatchIM
#     phi0 = Phase22.phi0
#     phi1 = Phase22.phi1
#     phi2 = Phase22.phi2
#     phi3 = Phase22.phi3
#     phi4 = Phase22.phi4
#     phi5 = Phase22.phi5
#     phi5L = Phase22.phi5L
#     phi6 = Phase22.phi6
#     phi6L = Phase22.phi6L
#     phi7 = Phase22.phi7
#     phi8 = Phase22.phi8
#     phi8L = Phase22.phi8L
#     phi9 = Phase22.phi9
#     phi9L = Phase22.phi9L
#     dphase0 = Phase22.dphase0
#     dphi0 = Phase22.dphi0
#     dphi1 = Phase22.dphi1
#     dphi2 = Phase22.dphi2
#     dphi3 = Phase22.dphi3
#     a0coloc = Phase22.a0coloc
#     a1coloc = Phase22.a1coloc
#     a2coloc = Phase22.a2coloc
#     a3coloc = Phase22.a3coloc
#     a4coloc = Phase22.a4coloc
#     b0coloc = Phase22.b0coloc
#     b1coloc = Phase22.b1coloc
#     b2coloc = Phase22.b2coloc
#     b3coloc = Phase22.b3coloc
#     b4coloc = Phase22.b4coloc
#     cLcoloc = Phase22.cLcoloc
#     C2Int = Phase22.C2Int
#     c0coloc = Phase22.c0coloc
#     c1coloc = Phase22.c1coloc
#     c2coloc = Phase22.c2coloc
#     c4coloc = Phase22.c4coloc
#     C2MRD = Phase22.C2MRD
#     sigma1 = Phase22.sigma1
#     sigma2 = Phase22.sigma2
#     sigma3 = Phase22.sigma3
#     sigma4 = Phase22.sigma4
#     sigma5 = Phase22.sigma5
#     C1Int = Phase22.C1Int
#     C1MRD = Phase22.C1MRD

    
#     c4ov3   = c4coloc / 3.
#     cLovfda = cLcoloc / fdamp
        
#     # Intermediate phase derivative at fPhaseMatchIN
#     DPhiInt = b0coloc + b1coloc/fPhaseMatchIN + b2coloc/(fPhaseMatchIN*fPhaseMatchIN) + b3coloc/(fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN) + b4coloc/(fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN) + (4. *cLcoloc) / ((4. *fdamp*fdamp) + (fPhaseMatchIN - fring)*(fPhaseMatchIN - fring))
    
#     C2Int = phaseIN - DPhiInt
#     # Inspiral phase at fPhaseMatchIN
#     phiIN = phi0 + phi1*(fPhaseMatchIN^(1. /3.)) + phi2*(fPhaseMatchIN^(2. /3.)) + phi3*fPhaseMatchIN + phi4*(fPhaseMatchIN^(4. /3.)) + (phi5 + phi5L*log(fPhaseMatchIN))*(fPhaseMatchIN^(5. /3.)) + (phi6 + phi6L*log(fPhaseMatchIN))*fPhaseMatchIN*fPhaseMatchIN + phi7*(fPhaseMatchIN^(7. /3.)) + (phi8 + phi8L*log(fPhaseMatchIN))*(fPhaseMatchIN^(8. /3.)) + (phi9  + phi9L*log(fPhaseMatchIN))*fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN
#     # add the pseudo-PN Coefficients
#     phiIN = phiIN + sigma1*(fPhaseMatchIN^(8. /3.)) + sigma2*(fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN) + sigma3*(fPhaseMatchIN^(10. /3.)) + sigma4*(fPhaseMatchIN^(11. /3.)) + sigma5*(fPhaseMatchIN^4)
#     # finally the normalisation
#     phiIN = phiIN * phiNorm*(fPhaseMatchIN^(-5. /3.))
    
#     # Intermediate phase at fPhaseMatchIN
    
#     phiIM = b0coloc*fPhaseMatchIN + b1coloc*log(fPhaseMatchIN) - b2coloc/fPhaseMatchIN - b3coloc/(fPhaseMatchIN*fPhaseMatchIN)/2. - (b4coloc/(fPhaseMatchIN*fPhaseMatchIN*fPhaseMatchIN)/3.) + (2. * cLcoloc * atan((fPhaseMatchIN - fring) / (2. * fdamp)))/fdamp
        
#     C1Int = phiIN - phiIM - (C2Int * fPhaseMatchIN)
    
#     # Intermediate phase at fPhaseMatchIM + connections
#     phiIMC   = b0coloc*fPhaseMatchIM + b1coloc*log(fPhaseMatchIM) - b2coloc/fPhaseMatchIM - b3coloc/(fPhaseMatchIM*fPhaseMatchIM)/2. - (b4coloc/(fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM)/3.) + (2. * cLcoloc * atan((fPhaseMatchIM - fring) / (2. * fdamp)))/fdamp + C1Int + C2Int*fPhaseMatchIM
#     # Ringdown phase at fPhaseMatchIM
#     phiRD    = c0coloc*fPhaseMatchIM + 1.5*c1coloc*(fPhaseMatchIM^(2. /3.)) - c2coloc/fPhaseMatchIM - c4ov3/(fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM) + (cLovfda * atan((fPhaseMatchIM - fring)/fdamp))
#     # Intermediate phase derivative at fPhaseMatchIM + connection
#     DPhiIntC = b0coloc + b1coloc/fPhaseMatchIM + b2coloc/(fPhaseMatchIM*fPhaseMatchIM) + b3coloc/(fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM) + b4coloc/(fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM) + (4. *cLcoloc) / ((4. *fdamp*fdamp) + (fPhaseMatchIM - fring)*(fPhaseMatchIM - fring)) + C2Int
#     # Ringdown phase derivative at fPhaseMatchIM
#     DPhiRD   = (c0coloc + c1coloc*(fPhaseMatchIM^(-1. /3.)) + c2coloc/(fPhaseMatchIM*fPhaseMatchIM) + c4coloc/(fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM*fPhaseMatchIM) + (cLcoloc / (fdamp*fdamp + (fPhaseMatchIM - fring)*(fPhaseMatchIM - fring))))
    
#     C2MRD = DPhiIntC - DPhiRD
#     C1MRD = phiIMC - phiRD - C2MRD*fPhaseMatchIM
    

#     # Linear time and phase shifts so that model peaks near t ~ 0
#     lina = 0.
#     linb         = ((3155.1635543201924 + 1257.9949740608242*eta - 32243.28428870599*eta2 + 347213.65466875216*eta2*eta - 1.9223851649491738e6*eta2*eta2 + 5.3035911346921865e6*eta2*eta2*eta - 5.789128656876938e6*eta2*eta2*eta2) + ((-24.181508118588667 + 115.49264174560281*eta - 380.19778216022763*eta2)*totchi + (24.72585609641552 - 328.3762360751952*eta + 725.6024119989094*eta2)*totchi2 + (23.404604124552 - 646.3410199799737*eta + 1941.8836639529036*eta2)*totchi2*totchi + (-12.814828278938885 - 325.92980012408367*eta + 1320.102640190539*eta2)*totchi2*totchi2) + (-148.17317525117338*dchi*delta*eta2))
#     _completePhaseDer = infreqs -> @. ifelse(infreqs <= fPhaseMatchIN, (infreqs^(-8. /3.))*dphase0*(dphi0 + dphi1*(infreqs^(1. /3.)) + dphi2*(infreqs^(2. /3.)) + dphi3*infreqs + dphi4*(infreqs^(4. /3.)) + dphi5*(infreqs^(5. /3.)) + (dphi6 + dphi6L*log(infreqs))*infreqs*infreqs + dphi7*(infreqs^(7. /3.)) + (dphi8 + dphi8L*log(infreqs))*(infreqs^(8. /3.)) + (dphi9  + dphi9L*log(infreqs))*infreqs*infreqs*infreqs + a0coloc*(infreqs^(8. /3.)) + a1coloc*infreqs*infreqs*infreqs + a2coloc*(infreqs^(10. /3.)) + a3coloc*(infreqs^(11. /3.)) + a4coloc*(infreqs^4)), ifelse(infreqs <= fPhaseMatchIM, b0coloc + b1coloc/infreqs + b2coloc/(infreqs*infreqs) + b3coloc/(infreqs*infreqs*infreqs) + b4coloc/(infreqs*infreqs*infreqs*infreqs) + (4. *cLcoloc) / ((4. *fdamp*fdamp) + (infreqs - fring)*(infreqs - fring)) + C2Int, (c0coloc + c1coloc*(infreqs^(-1. /3.)) + c2coloc/(infreqs*infreqs) + c4coloc/(infreqs*infreqs*infreqs*infreqs) + (cLcoloc / (fdamp*fdamp + (infreqs - fring)*(infreqs - fring)))) + C2MRD))
#     _completePhase = infreqs ->  @. ifelse(infreqs <= fPhaseMatchIN, phiNorm*(infreqs^(-5. /3.))*(phi0 + phi1*(infreqs^(1. /3.)) + phi2*(infreqs^(2. /3.)) + phi3*infreqs + phi4*(infreqs^(4. /3.)) + (phi5 + phi5L*log(infreqs))*(infreqs^(5. /3.)) + (phi6 + phi6L*log(infreqs))*infreqs*infreqs + phi7*(infreqs^(7. /3.)) + (phi8 + phi8L*log(infreqs))*(infreqs^(8. /3.)) + (phi9  + phi9L*log(infreqs))*infreqs*infreqs*infreqs + sigma1*(infreqs^(8. /3.)) + sigma2*(infreqs*infreqs*infreqs) + sigma3*(infreqs^(10. /3.)) + sigma4*(infreqs^(11. /3.)) + sigma5*(infreqs^4)), ifelse(infreqs <= fPhaseMatchIM, b0coloc*infreqs + b1coloc*log(infreqs) - b2coloc/infreqs - b3coloc/(infreqs*infreqs)/2. - (b4coloc/(infreqs*infreqs*infreqs)/3.) + (2. * cLcoloc * atan((infreqs - fring) / (2. * fdamp)))/fdamp + C1Int + C2Int*infreqs, ifelse(infreqs < fcutPar, (c0coloc*infreqs + 1.5*c1coloc*(infreqs^(2. /3.)) - c2coloc/infreqs - c4ov3/(infreqs*infreqs*infreqs) + (cLovfda * atan((infreqs - fring)/fdamp))) + C1MRD + C2MRD*infreqs, 0.)))

#     dphi22Ref    = etaInv * _completePhaseDer(fring-fdamp)
#     psi4tostrain = ((13.39320482758057 - 175.42481512989315*eta + 2097.425116152503*eta2 - 9862.84178637907*eta2*eta + 16026.897939722587*eta2*eta2) + ((4.7895602776763 - 163.04871764530466*eta + 609.5575850476959*eta2)*totchi + (1.3934428041390161 - 97.51812681228478*eta + 376.9200932531847*eta2)*totchi2 + (15.649521097877374 + 137.33317057388916*eta - 755.9566456906406*eta2)*totchi2*totchi + (13.097315867845788 + 149.30405703643288*eta - 764.5242164872267*eta2)*totchi2*totchi2) + (105.37711654943146*dchi*Seta*eta2))
#     linb         = linb - dphi22Ref - 2. *pi*(500. +psi4tostrain)
    
#     phifRef = -(etaInv*_completePhase(fRef) + linb*fRef + lina) + pi/4. + pi
#     phis    = @. etaInv*_completePhase(fgrid) + ifelse(fgrid <= fcutPar, linb*fgrid + lina + phifRef, 0.)
    
#     return phis

# end

function Ampl_22_ConnectionCoefficients(
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar=0.3,
    IntAmpVersion=104,
    GMsun_over_c3= uc.GMsun_over_c3,
    GMsun_over_c2_Gpc= uc.GMsun_over_c2_Gpc,
    final_spin_override=nothing,

)
    """
    Compute the amplitude of the GW as a function of frequency, given the events parameters.

    :param numpy.ndarray f: Frequency grid on which the phase will be computed, in :math:`\\rm Hz`.
    :param dict(numpy.ndarray, numpy.ndarray, ...) kwargs: Dictionary with arrays containing the parameters of the events to compute the amplitude of, as in :py:data:`events`.
    :return: GW amplitude for the chosen events evaluated on the frequency grid.
    :rtype: numpy.ndarray

    """
    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)

    # We work in dimensionless frequency M*f, not f

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    q = 0.5 * (1.0 + Seta - 2.0 * eta) / eta
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    m1ByMSq = m1ByM * m1ByM
    m2ByMSq = m2ByM * m2ByM
    # PN symmetry coefficient
    delta = Seta
    totchi = ((m1ByMSq * chi1 + m2ByMSq * chi2) / (m1ByMSq + m2ByMSq))
    totchi2 = totchi * totchi
    dchi = chi1 - chi2
    # Normalised PN reduced spin parameter
    chi_eff = (m1ByM * chi1 + m2ByM * chi2)
    chiPN = (chi_eff - (38.0 / 113.0) * eta * (chi1 + chi2)) / (1.0 - (76.0 * eta / 113.0))
    chiPN2 = chiPN * chiPN

    amp0 = 2.0 * sqrt(5.0 / (64.0 * pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL
    ampNorm = sqrt(2.0 * eta / 3.0) * (pi^(-1.0 / 6.0))

    # Compute final spin and radiated energy
    aeff = (
        (
            (
                3.4641016151377544 * eta + 20.0830030082033 * eta2 -
                12.333573402277912 * eta2 * eta
            ) / (1 + 7.2388440419467335 * eta)
        ) +
        (
            (m1ByMSq + m2ByMSq) * totchi +
            (
                (
                    -0.8561951310209386 * eta - 0.09939065676370885 * eta2 +
                    1.668810429851045 * eta2 * eta
                ) * totchi +
                (
                    0.5881660363307388 * eta - 2.149269067519131 * eta2 +
                    3.4768263932898678 * eta2 * eta
                ) * totchi2 +
                (
                    0.142443244743048 * eta - 0.9598353840147513 * eta2 +
                    1.9595643107593743 * eta2 * eta
                ) *
                totchi2 *
                totchi
            ) / (
                1 +
                (
                    -0.9142232693081653 + 2.3191363426522633 * eta -
                    9.710576749140989 * eta2 * eta
                ) * totchi
            )
        ) +
        (
            0.3223660562764661 * dchi * Seta * (1 + 9.332575956437443 * eta) * eta2 -
            0.059808322561702126 * dchi * dchi * eta2 * eta +
            2.3170397514509933 *
            dchi *
            Seta *
            (1 - 3.2624649875884852 * eta) *
            eta2 *
            eta *
            totchi
        )
    )
    aeff_nonprec = aeff
    aeff = isnothing(final_spin_override) ? aeff : final_spin_override
    Erad = (
        (
            (
                (
                    0.057190958417936644 * eta + 0.5609904135313374 * eta2 -
                    0.84667563764404 * eta2 * eta + 3.145145224278187 * eta2 * eta2
                ) * (
                    1.0 +
                    (
                        -0.13084389181783257 - 1.1387311580238488 * eta +
                        5.49074464410971 * eta2
                    ) * totchi +
                    (-0.17762802148331427 + 2.176667900182948 * eta2) * totchi2 +
                    (
                        -0.6320191645391563 + 4.952698546796005 * eta -
                        10.023747993978121 * eta2
                    ) *
                    totchi *
                    totchi2
                )
            ) / (
                1.0 +
                (-0.9919475346968611 + 0.367620218664352 * eta + 4.274567337924067 * eta2) *
                totchi
            )
        ) + (
            -0.09803730445895877 * dchi * Seta * (1.0 - 3.2283713377939134 * eta) * eta2 +
            0.01118530335431078 * dchi * dchi * eta2 * eta -
            0.01978238971523653 *
            dchi *
            Seta *
            (1.0 - 4.91667749015812 * eta) *
            eta *
            totchi
        )
    )
    # Compute ringdown and damping frequencies from fits
    fring =
        (
            (
                0.05947169566573468 - 0.14989771215394762 * aeff +
                0.09535606290986028 * aeff * aeff +
                0.02260924869042963 * aeff * aeff * aeff -
                0.02501704155363241 * aeff * aeff * aeff * aeff -
                0.005852438240997211 * (aeff^5) +
                0.0027489038393367993 * (aeff^6) +
                0.0005821983163192694 * (aeff^7)
            ) / (
                1 - 2.8570126619966296 * aeff + 2.373335413978394 * aeff * aeff -
                0.6036964688511505 * aeff * aeff * aeff * aeff +
                0.0873798215084077 * (aeff^6)
            )
        ) / (1.0 - Erad)
    fdamp =
        (
            (
                0.014158792290965177 - 0.036989395871554566 * aeff +
                0.026822526296575368 * aeff * aeff +
                0.0008490933750566702 * aeff * aeff * aeff -
                0.004843996907020524 * aeff * aeff * aeff * aeff -
                0.00014745235759327472 * (aeff^5) + 0.0001504546201236794 * (aeff^6)
            ) / (
                1 - 2.5900842798681376 * aeff + 1.8952576220623967 * aeff * aeff -
                0.31416610693042507 * aeff * aeff * aeff * aeff +
                0.009002719412204133 * (aeff^6)
            )
        ) / (1.0 - Erad)

    aeff = aeff_nonprec
    # Fitting function for hybrid minimum energy circular orbit (MECO) function and computation of ISCO frequency
    Z1tmp = 1.0 + cbrt((1.0 - aeff * aeff)) * (cbrt(1.0 + aeff) + cbrt(1.0 - aeff))
    Z1tmp = ifelse(Z1tmp > 3.0, 3.0, Z1tmp)
    Z2tmp = sqrt(3.0 * aeff * aeff + Z1tmp * Z1tmp)
    fISCO =
        (
            1.0 / (
                (
                    3.0 + Z2tmp -
                    sign(aeff) * sqrt((3.0 - Z1tmp) * (3.0 + Z1tmp + 2.0 * Z2tmp))
                )^(3.0 / 2.0) + aeff
            )
        ) / pi

    fMECO = (
        (
            (
                0.018744340279608845 +
                0.0077903147004616865 * eta +
                0.003940354686136861 * eta2 - 0.00006693930988501673 * eta2 * eta
            ) / (1.0 - 0.10423384680638834 * eta)
        ) +
        (
            (
                chiPN * (
                    0.00027180386951683135 - 0.00002585252361022052 * chiPN +
                    eta2 *
                    eta2 *
                    (
                        -0.0006807631931297156 + 0.022386313074011715 * chiPN -
                        0.0230825153005985 * chiPN2
                    ) +
                    eta2 * (
                        0.00036556167661117023 - 0.000010021140796150737 * chiPN -
                        0.00038216081981505285 * chiPN2
                    ) +
                    eta * (
                        0.00024422562796266645 - 0.00001049013062611254 * chiPN -
                        0.00035182990586857726 * chiPN2
                    ) +
                    eta2 *
                    eta *
                    (
                        -0.0005418851224505745 +
                        0.000030679548774047616 * chiPN +
                        4.038390455349854e-6 * chiPN2
                    ) - 0.00007547517256664526 * chiPN2
                )
            ) / (
                0.026666543809890402 +
                (
                    -0.014590539285641243 - 0.012429476486138982 * eta +
                    1.4861197211952053 * eta2 * eta2 +
                    0.025066696514373803 * eta2 +
                    0.005146809717492324 * eta2 * eta
                ) * chiPN +
                (
                    -0.0058684526275074025 - 0.02876774751921441 * eta -
                    2.551566872093786 * eta2 * eta2 - 0.019641378027236502 * eta2 -
                    0.001956646166089053 * eta2 * eta
                ) * chiPN2 +
                (
                    0.003507640638496499 +
                    0.014176504653145768 * eta +
                    1.0 * eta2 * eta2 +
                    0.012622225233586283 * eta2 - 0.00767768214056772 * eta2 * eta
                ) *
                chiPN2 *
                chiPN
            )
        ) +
        (
            dchi * dchi * (0.00034375176678815234 + 0.000016343732281057392 * eta) * eta2 +
            dchi *
            Seta *
            eta *
            (
                0.08064665214195679 * eta2 +
                eta * (-0.028476219509487793 - 0.005746537021035632 * chiPN) -
                0.0011713735642446144 * chiPN
            )
        )
    )

    gamma2 = (
        (
            (0.8312293675316895 + 7.480371544268765 * eta - 18.256121237800397 * eta2) /
            (1.0 + 10.915453595496611 * eta - 30.578409433912874 * eta2)
        ) +
        (
            (
                totchi * (
                    0.5869408584532747 +
                    eta * (-0.1467158405070222 - 2.8489481072076472 * totchi) +
                    0.031852563636196894 * totchi +
                    eta2 * (0.25295441250444334 + 4.6849496672664594 * totchi)
                )
            ) / (3.8775263105069953 - 3.41755361841226 * totchi + 1.0 * totchi2)
        ) +
        (-0.00548054788508203 * dchi * delta * eta)
    )
    gamma3 = (
        (
            (
                1.3666000000000007 - 4.091333144596439 * eta + 2.109081209912545 * eta2 -
                4.222259944408823 * eta2 * eta
            ) / (1.0 - 2.7440263888207594 * eta)
        ) +
        (
            (
                0.07179105336478316 +
                eta2 * (2.331724812782498 - 0.6330998412809531 * totchi) +
                eta * (-0.8752427297525086 + 0.4168560229353532 * totchi) -
                0.05633734476062242 * totchi
            ) * totchi
        ) +
        (0.0 * delta * dchi)
    )
    fAmpRDMin = ifelse(
        gamma2 <= 1.0,
        abs(fring + fdamp * gamma3 * (sqrt(1.0 - gamma2 * gamma2) - 1.0) / gamma2),
        abs(fring + fdamp * (-1.0) * gamma3 / gamma2),
    )
    v1RD = (
        (
            (0.03689164742964719 + 25.417967754401182 * eta + 162.52904393600332 * eta2) /
            (1.0 + 61.19874463331437 * eta - 29.628854485544874 * eta2)
        ) +
        (
            (
                totchi * (
                    -0.14352506969368556 +
                    0.026356911108320547 * totchi +
                    0.19967405175523437 * totchi2 - 0.05292913111731128 * totchi2 * totchi +
                    eta2 *
                    eta *
                    (
                        -48.31945248941757 - 3.751501972663298 * totchi +
                        81.9290740950083 * totchi2 +
                        30.491948143930266 * totchi2 * totchi -
                        132.77982622925845 * totchi2 * totchi2
                    ) +
                    eta * (
                        -4.805034453745424 +
                        1.11147906765112 * totchi +
                        6.176053843938542 * totchi2 -
                        0.2874540719094058 * totchi2 * totchi -
                        8.990840289951514 * totchi2 * totchi2
                    ) - 0.18147275151697131 * totchi2 * totchi2 +
                    eta2 * (
                        27.675454081988036 - 2.398327419614959 * totchi -
                        47.99096500250743 * totchi2 - 5.104257870393138 * totchi2 * totchi +
                        72.08174136362386 * totchi2 * totchi2
                    )
                )
            ) / (-1.4160870461211452 + 1.0 * totchi)
        ) +
        (-0.04426571511345366 * dchi * delta * eta2)
    )
    F1 = fAmpRDMin
    gamma1 =
        (v1RD / (fdamp * gamma3)) *
        (F1 * F1 - 2.0 * F1 * fring + fring * fring + fdamp * fdamp * gamma3 * gamma3) *
        exp(((F1 - fring) * gamma2) / (fdamp * gamma3))
    
    Ampl22 = Ampl22Struct(gamma1, gamma2, gamma3)

    return Ampl22
end
