#! format: off
#######################################################
# IMRPhenomXE WAVEFORM
#######################################################
#
# Frequency-domain (2,2)-mode eccentric aligned-spin model IMRPhenomXE (arXiv:2601.03340), ported from the
# phenomxpy implementation (phenomxpy/phenomxe, branch phenomxehm_dev). The model is built on top of IMRPhenomXAS:
#
#   * the quasi-Keplerian (3PN EOB) secular dynamics (x, e, l, lambda) are integrated numerically, with the
#     quasi-circular part of dx/dt taken from the IMRPhenomT (2,2) frequency evolution (see PhenomXE_PhenomT22.jl
#     and PhenomXE_dynamics.jl);
#   * the inspiral phase of the j = 0 (circular) harmonic is obtained from the numerical stationary phase
#     approximation (SPA) of the dynamics and attached to the IMRPhenomXAS intermediate and ringdown phase with
#     newly computed connection coefficients;
#   * the j = 0 amplitude is the IMRPhenomXAS amplitude plus a windowed eccentric correction built from the 3PN
#     eccentric time-domain amplitude of the 22 mode and the NR-calibrated IMRPhenomT amplitude;
#   * the j != 0 mean-anomaly harmonics (|j| <= n_harmonics) are added with their own numerical SPA
#     (PhenomXE_amplitudes.jl contains the 3PN amplitude coefficients, expanded to O(e^12)).
#
# Only the "dense" reference algorithm of phenomxpy is ported (no multibanding or sparse assembly), evaluated
# directly on the user-supplied (arbitrary, increasing) frequency grid.
#
# Additional parameters with respect to the quasi-circular models: `ecc` (eccentricity) and `meanAno`
# (mean anomaly, in [0, 2 pi]) defined at the reference frequency `fRef` (default: the first frequency of the grid).

"""
Compute the plus and cross polarizations of IMRPhenomXE on the frequency grid `f` (in Hz).

    hphc(PhenomXE(), f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Increasing frequency grid, in Hz.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the binary, in Gpc.
-  `iota`: Inclination angle, in radians.
-  `ecc`: Eccentricity at the reference frequency.
-  `meanAno`: Mean anomaly at the reference frequency, in radians.
#### Optional arguments:
-  `fRef`: Reference frequency in Hz at which `ecc` and `meanAno` are defined. Default is the first frequency of the grid.
-  `phiRef`: Reference orbital phase, in radians. Default is 0.
-  `n_harmonics`: Number of eccentric mean-anomaly harmonics (j = -n_harmonics, ..., n_harmonics). Default is `model.n_harmonics` (6).
-  `fcutPar`: Dimensionless frequency (Mf) at which the waveform ends. Default is 0.3.
-  `rtol`, `atol`: Tolerances of the ODE integration of the eccentric dynamics. Default is 1e-12.
-  `dt_fine`, `dtheta_min`, `nn_max`: Parameters of the PhenomT time grid used to build dx/dt (defaults 0.5, 0.002, 10000).
-  `e_merger_max`: Maximum eccentricity allowed at the end of the dynamics. Default is 0.2.
#### Return:
-  `hp`, `hc`: complex arrays with the plus and cross polarizations (the full phase is included).
"""
function hphc(model::PhenomXE, f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno;
              fRef = nothing, phiRef = 0.0, n_harmonics = model.n_harmonics, fcutPar = 0.3,
              GMsun_over_c3 = uc.GMsun_over_c3, GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
              rtol = 1e-12, atol = 1e-12, dt_fine = 0.5, dtheta_min = 0.002, nn_max = 10000, e_merger_max = 0.2,
              container = nothing, call_number = 2, optimization = false, debug = false)

    h22, active, _ = _phenomXE_h22(model, f, mc, eta, chi1, chi2, dL, ecc, meanAno;
                                    fRef = fRef, phiRef = phiRef, n_harmonics = n_harmonics, fcutPar = fcutPar,
                                    GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc,
                                    rtol = rtol, atol = atol, dt_fine = dt_fine, dtheta_min = dtheta_min, nn_max = nn_max,
                                    e_merger_max = e_merger_max, debug = debug)

    # Project the 22 mode on the polarizations (same convention as phenomxpy / LAL PhenomX)
    cos_inclination = cos(iota)
    ylm_factor = sqrt(5.0 / (64.0 * pi))
    hp_factor = -ylm_factor * (1.0 + cos_inclination * cos_inclination)
    hc_factor = 2.0im * ylm_factor * cos_inclination
    hp = hp_factor .* h22
    hc = hc_factor .* h22

    if !isnothing(container) && eltype(real(hp[1])) !== Float64
        for i in eachindex(f)
            container[i] = real(hp[i]).value + 1im * imag(hp[i]).value
            container[i + length(f)] = real(hc[i]).value + 1im * imag(hc[i]).value
        end
    end

    if optimization == true
        if call_number == 1
            return [real(hp); imag(hp); real(hc); imag(hc)]
        else
            return [hp; hc]
        end
    end
    return hp, hc
end

"""
Compute the complex IMRPhenomXE (2,2) mode on the frequency grid `f`, together with the mask of active
frequencies (Mf <= Mf_cut) and a NamedTuple with intermediate quantities (for debugging and tests).
See `hphc(model::PhenomXE, ...)` for the arguments.
"""
function _phenomXE_h22(model::PhenomXE, f, mc, eta, chi1, chi2, dL, ecc, meanAno;
                       fRef = nothing, phiRef = 0.0, n_harmonics = model.n_harmonics, fcutPar = 0.3,
                       GMsun_over_c3 = uc.GMsun_over_c3, GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
                       rtol = 1e-12, atol = 1e-12, dt_fine = 0.5, dtheta_min = 0.002, nn_max = 10000, e_merger_max = 0.2,
                       debug = false)

    # this is the type of the variables when ForwardDiff is used (keeps track of the Dual tag)
    typeofFD = promote_type(typeof(mc), typeof(eta), typeof(chi1), typeof(chi2), typeof(dL), typeof(ecc), typeof(meanAno), typeof(phiRef))

    eta = ifelse(eta > 0.25, 0.25 * one(eta), eta)
    M = mc / (eta^(0.6))
    eta2 = eta * eta
    etaInv = 1.0 / eta
    Seta = _xe_delta(eta)
    delta = Seta
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    chiS = 0.5 * (chi1 + chi2)
    chiA = 0.5 * (chi1 - chi2)
    chi_eff = m1ByM * chi1 + m2ByM * chi2
    totchi = (m1ByM * m1ByM * chi1 + m2ByM * m2ByM * chi2) / (m1ByM * m1ByM + m2ByM * m2ByM)
    dchi = chi1 - chi2
    mm = 2

    # Frequency grid in dimensionless units
    issorted(f) || error("PhenomXE: the frequency grid must be increasing")
    fgrid = M * GMsun_over_c3 .* f
    len = length(fgrid)
    f_min = f[1]
    Mfmin = fgrid[1]
    fRef_ = isnothing(fRef) ? f_min : fRef
    if fRef_ < f_min
        # as in phenomxpy: f_ref < f_min changes the starting frequency to f_ref
        Mfmin = fRef_ * M * GMsun_over_c3
    end
    MfRef = fRef_ * M * GMsun_over_c3

    # End of the waveform. As in phenomxpy, Mf = 0.33 for very high effective spins
    fCutDef = ifelse(chi_eff > 0.99, 0.33, fcutPar)
    Mf_max_prime = min(fCutDef, fgrid[end])
    if fCutDef <= Mfmin
        error("PhenomXE: fCut = $(fCutDef / (M * GMsun_over_c3)) Hz <= f_min = $(f_min) Hz")
    end
    active = fgrid .<= Mf_max_prime
    Mfs_xas = fgrid[active]
    nact = length(Mfs_xas)

    # Amplitude normalisations (phenomxpy conventions)
    amp0 = M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL
    amp0_xas = 2.0 * sqrt(5.0 / (64.0 * pi)) * amp0

    phi0 = phiRef

    # ------------------------------------------------------------------ IMRPhenomXAS quantities
    rem = _xas_remnant_frequencies(eta, chi1, chi2)
    fring = rem.fring
    fdamp = rem.fdamp
    fMECO = rem.fMECO
    fISCO = rem.fISCO

    Phase22 = Phase_22_ConnectionCoefficients(mc, eta, chi1, chi2; fcutPar = fCutDef, GMsun_over_c3 = GMsun_over_c3)
    Amp22 = Ampl_22_ConnectionCoefficients(mc, eta, chi1, chi2, dL; GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc)
    gamma2 = Amp22.gamma2
    gamma3 = Amp22.gamma3
    fAmpMatchIN = fMECO + 0.25 * (fISCO - fMECO)
    fAmpRDMin = ifelse(gamma2 <= 1.0,
                       abs(fring + fdamp * gamma3 * (sqrt(1.0 - gamma2 * gamma2) - 1.0) / gamma2),
                       abs(fring + fdamp * (-1.0) * gamma3 / gamma2))
    fAmpIntMax = fAmpRDMin

    fInsp = Phase22.fPhaseMatchIN
    fInt = Phase22.fPhaseMatchIM
    phase_int_list = (Phase22.b0coloc, Phase22.b1coloc, Phase22.b2coloc, Phase22.b3coloc, Phase22.b4coloc, Phase22.cLcoloc)
    c0 = Phase22.c0coloc; c1 = Phase22.c1coloc; c2 = Phase22.c2coloc; c4 = Phase22.c4coloc; cL = Phase22.cLcoloc
    c4ov3 = c4 / 3.0
    cLovfda = cL / fdamp

    # Time shift and reference phase of IMRPhenomXAS (used by the XAS fallback for very short signals)
    lina = 0.0
    linb = TimeShift_22(model, eta, Seta, totchi, dchi, fring, fdamp, Phase22)
    phifRef_xas_anchor = -(etaInv * _completePhase(model, MfRef, Phase22, fdamp, fring, fCutDef) + linb * MfRef + lina) + pi / 4.0
    phifRef = phifRef_xas_anchor + 2.0 * phi0
    phifRef_xe = 0.0

    # IMRPhenomXAS IMR amplitude (without the overall amp0 normalisation, as in phenomxpy)
    amp_xas = Ampl(PhenomXAS(), f[active], mc, eta, chi1, chi2, dL; fcutPar = fCutDef, GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc) ./ amp0_xas

    # ------------------------------------------------------------------ Eccentric dynamics and numerical SPA
    x_start = (pi * Mfmin)^(2.0 / 3.0)
    x_ref = (pi * MfRef)^(2.0 / 3.0)
    dynT = _compute_qkp_dynamics(eta, delta, chi1, chi2, chiA, chiS, x_start, x_ref, ecc, meanAno, Mfmin;
                                 rtol = rtol, atol = atol, dt_fine = dt_fine, dtheta_min = dtheta_min, nn_max = nn_max, e_merger_max = e_merger_max)
    spa = _compute_num_spa_j0_phase(dynT; mm = mm)
    phT = dynT.phT
    Mfs_spa = spa.Mfs_spa
    phase_j_spa = spa.phase_j_spa
    tc = spa.tc

    # ------------------------------------------------------------------ j = 0 phase
    phase_xas = false
    amp_xas_only = false
    Mf_min_spa = Mfs_spa[1]
    f_spa_cut = fInsp
    if Mf_min_spa < fInsp
        f_spa_cut = fInsp
    elseif Mf_min_spa > fInsp && Mf_min_spa < 0.9 * fInt
        f_spa_cut = 0.9 * fInt
    else
        # the dynamics starts above the XAS intermediate region: fall back to IMRPhenomXAS
        phase_xas = true
        amp_xas_only = true
    end

    phi_j0 = Vector{typeofFD}(undef, nact)
    C1Int = Phase22.C1Int; C2Int = Phase22.C2Int; C1MRD = Phase22.C1MRD; C2MRD = Phase22.C2MRD
    if phase_xas
        phi_j0 .= etaInv .* _completePhase(model, Mfs_xas, Phase22, fdamp, fring, fCutDef) .+ linb .* Mfs_xas .+ lina .+ phifRef
    else
        # Inspiral: numerical SPA phase (rescaled by -eta so that it can be added to the XAS phase convention)
        phase_spa_rescaled = -eta .* (phase_j_spa .+ linb .* Mfs_spa .+ lina .+ phifRef_xe)
        cph, _ = _notaknot_spline(Mfs_spa, phase_spa_rescaled)
        nseg = length(Mfs_spa) - 1
        segment = clamp(searchsortedlast(Mfs_spa, f_spa_cut), 1, nseg)
        offset = f_spa_cut - Mfs_spa[segment]
        phase_spa_last = _spline_eval_segment(cph, segment, offset)
        derivative_spa_last = _spline_deriv_segment(cph, segment, offset)

        # New inspiral-intermediate connection coefficients
        DPhiInt = _xe_intermediate_dphase(f_spa_cut, fdamp, fring, phase_int_list)
        C2Int = derivative_spa_last - DPhiInt
        phiIM = _xe_intermediate_phase(f_spa_cut, fring, fdamp, phase_int_list)
        C1Int = phase_spa_last - phiIM - C2Int * f_spa_cut

        # New intermediate-ringdown connection coefficients
        phase_end = _xe_intermediate_phase(fInt, fring, fdamp, phase_int_list) + C1Int + C2Int * fInt
        derivative_end = _xe_intermediate_dphase(fInt, fdamp, fring, phase_int_list) + C2Int
        phiRD = _xe_ringdown_phase(fInt, fdamp, fring, c0, c1, c2, c4ov3, cLovfda)
        DPhiRD = _xe_ringdown_dphase(fInt, fdamp, fring, c0, c1, c2, c4, cL)
        C2MRD = derivative_end - DPhiRD
        C1MRD = phase_end - phiRD - C2MRD * fInt

        # Evaluate the IMR phase on the grid
        inspiral_end = searchsortedlast(Mfs_xas, f_spa_cut)
        intermediate_end = searchsortedlast(Mfs_xas, fInt)
        seg = 1
        for idx in 1:inspiral_end
            Mf = Mfs_xas[idx]
            while seg < nseg && Mf >= Mfs_spa[seg + 1]
                seg += 1
            end
            value = _spline_eval_segment(cph, seg, Mf - Mfs_spa[seg])
            phi_j0[idx] = value * etaInv + linb * Mf + lina + phifRef_xe
        end
        for idx in (inspiral_end + 1):intermediate_end
            Mf = Mfs_xas[idx]
            value = _xe_intermediate_phase(Mf, fring, fdamp, phase_int_list) + C1Int + C2Int * Mf
            phi_j0[idx] = value * etaInv + linb * Mf + lina + phifRef_xe
        end
        for idx in (intermediate_end + 1):nact
            Mf = Mfs_xas[idx]
            value = _xe_ringdown_phase(Mf, fdamp, fring, c0, c1, c2, c4ov3, cLovfda) + C1MRD + C2MRD * Mf
            phi_j0[idx] = value * etaInv + linb * Mf + lina + phifRef_xe
        end
    end

    # ------------------------------------------------------------------ j = 0 amplitude
    amp_j0 = Vector{Complex{typeofFD}}(undef, nact)
    if amp_xas_only
        amp_j0 .= amp_xas
    else
        t_insp = spa.tt; xt_insp = spa.xt; et_insp = spa.et
        dxdt = spa.dxdt; v_spa_qc = spa.v_spa_qc
        Mfs_spa_amp = Mfs_spa
        if Mfs_spa[end] > f_spa_cut
            # cut the dynamics up to f_spa_cut, keeping one extra point at the high-frequency end to improve the
            # interpolation (phenomxpy resize_dynamics). NOTE: phenomxpy shifts all indices by one and therefore also
            # drops the first SPA sample, which zeroes the j = 0 amplitude in the grid bins below the second dynamics
            # point (in practice the first bin at f_min); here the first sample is kept.
            idx_insp = 1:min(searchsortedlast(Mfs_spa, f_spa_cut) + 1, length(Mfs_spa))
            if length(idx_insp) > 2
                Mfs_spa_amp = Mfs_spa[idx_insp]
                t_insp = t_insp[idx_insp]; xt_insp = xt_insp[idx_insp]; et_insp = et_insp[idx_insp]
                dxdt = dxdt[idx_insp]; v_spa_qc = v_spa_qc[idx_insp]
            end
        end
        v_spa = sqrt.(xt_insp)
        # time derivative of the PhenomT (quasi-circular) and eccentric orbital frequencies
        om_dot_qc = 1.5 .* [_spline_eval(dynT.sc, dynT.sx, x) for x in xt_insp] .* v_spa_qc
        om_dot = 1.5 .* dxdt .* v_spa

        # Inspiral-merger IMRPhenomT amplitude evaluated along the dynamics
        omega_cut = _phT_imr_omega(phT, phT.amp_inspiral_cut)
        xt_cut = (omega_cut / 2.0)^(2.0 / 3.0)
        amp22_qc = _xe_insp_merg_PhenomT_amplitude(t_insp, xt_insp, xt_cut, phT)

        # Eccentric corrections to the SPA amplitude of the j = 0 harmonic
        amp_ecc = _xe_amp_ecc_corrections_j0(eta, delta, chiA, chiS, et_insp, v_spa, om_dot, om_dot_qc, amp22_qc, mm)

        # Window the eccentric corrections and add them to the IMRPhenomXAS amplitude
        _xe_windowed_amplitude!(amp_j0, Mfs_xas, Mfs_spa_amp, f_spa_cut, amp_xas, amp_ecc, 1000.0)
    end

    # ------------------------------------------------------------------ j != 0 harmonics
    h22 = Vector{Complex{typeofFD}}(undef, len)
    fill!(h22, zero(Complex{typeofFD}))
    for (k, idx) in enumerate(findall(active))
        h22[idx] = amp0 * amp_j0[k] * exp(1im * phi_j0[k])
    end
    if n_harmonics > 0 && !amp_xas_only
        n_ecc = searchsortedlast(Mfs_xas, fAmpIntMax)
        Mfs_ecc = Mfs_xas[1:n_ecc]
        hj_sum = _xe_eccentric_harmonics(eta, delta, chiA, chiS, tc, n_harmonics, spa, Mfs_ecc, mm)
        act_idx = findall(active)
        for k in 1:n_ecc
            h22[act_idx[k]] += amp0 * hj_sum[k]
        end
    end

    # Rotate the mode by the reference phase
    if !phase_xas && phi0 != 0.0
        h22 .*= exp(2im * phi0)
    end

    if debug
        println("PhenomXE: Mtot = $(_val(M)), eta = $(_val(eta)), chi1 = $(_val(chi1)), chi2 = $(_val(chi2))")
        println("PhenomXE: fring = $(_val(fring)), fdamp = $(_val(fdamp)), fMECO = $(_val(fMECO)), fISCO = $(_val(fISCO))")
        println("PhenomXE: fPhaseMatchIN = $(_val(fInsp)), fPhaseMatchIM = $(_val(fInt)), f_spa_cut = $(_val(f_spa_cut))")
        println("PhenomXE: fAmpMatchIN = $(_val(fAmpMatchIN)), fAmpRDMin = $(_val(fAmpRDMin))")
        println("PhenomXE: dynamics points = $(length(spa.tt)), Mfs_spa in [$(_val(Mfs_spa[1])), $(_val(Mfs_spa[end]))], tc = $(_val(tc))")
        println("PhenomXE: phase_xas = $phase_xas, amp_xas = $amp_xas_only")
        println("PhenomXE: C1Int = $(_val(C1Int)), C2Int = $(_val(C2Int)), C1MRD = $(_val(C1MRD)), C2MRD = $(_val(C2MRD))")
        println("PhenomXE: linb = $(_val(linb)), phifRef = $(_val(phifRef))")
    end

    info = (Mfs_xas = Mfs_xas, phi_j0 = phi_j0, amp_j0 = amp_j0, amp_xas = amp_xas, Mfs_spa = Mfs_spa, phase_j_spa = phase_j_spa,
            spa = spa, dynT = dynT, Phase22 = Phase22, f_spa_cut = f_spa_cut, fInt = fInt, fInsp = fInsp, fring = fring, fdamp = fdamp,
            linb = linb, phifRef = phifRef, C1Int = C1Int, C2Int = C2Int, C1MRD = C1MRD, C2MRD = C2MRD, amp0 = amp0,
            phase_xas = phase_xas, amp_xas_only = amp_xas_only, fAmpIntMax = fAmpIntMax)
    return h22, active, info
end

##############################################################################
#   IMRPhenomXAS auxiliary quantities
##############################################################################

# Final spin, radiated energy, ringdown/damping frequencies, MECO and ISCO frequencies of IMRPhenomXAS
function _xas_remnant_frequencies(eta, chi1, chi2)
    eta2 = eta * eta
    Seta = _xe_delta(eta)
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    m1ByMSq = m1ByM * m1ByM
    m2ByMSq = m2ByM * m2ByM
    totchi = ((m1ByMSq * chi1 + m2ByMSq * chi2) / (m1ByMSq + m2ByMSq))
    totchi2 = totchi * totchi
    dchi = chi1 - chi2
    chi_eff = (m1ByM * chi1 + m2ByM * chi2)
    chiPN = (chi_eff - (38. /113.)*eta*(chi1 + chi2)) / (1. - (76. *eta/113.))
    chiPN2 = chiPN * chiPN

    aeff   = (((3.4641016151377544*eta + 20.0830030082033*eta2 - 12.333573402277912*eta2*eta)/(1 + 7.2388440419467335*eta)) + ((m1ByMSq + m2ByMSq)*totchi + ((-0.8561951310209386*eta - 0.09939065676370885*eta2 + 1.668810429851045*eta2*eta)*totchi + (0.5881660363307388*eta - 2.149269067519131*eta2 + 3.4768263932898678*eta2*eta)*totchi2 + (0.142443244743048*eta - 0.9598353840147513*eta2 + 1.9595643107593743*eta2*eta)*totchi2*totchi) / (1 + (-0.9142232693081653 + 2.3191363426522633*eta - 9.710576749140989*eta2*eta)*totchi)) + (0.3223660562764661*dchi*Seta*(1 + 9.332575956437443*eta)*eta2 - 0.059808322561702126*dchi*dchi*eta2*eta + 2.3170397514509933*dchi*Seta*(1 - 3.2624649875884852*eta)*eta2*eta*totchi))
    Erad   = ((((0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta2*eta + 3.145145224278187*eta2*eta2)*(1. + (-0.13084389181783257 - 1.1387311580238488*eta + 5.49074464410971*eta2)*totchi + (-0.17762802148331427 + 2.176667900182948*eta2)*totchi2 + (-0.6320191645391563 + 4.952698546796005*eta - 10.023747993978121*eta2)*totchi*totchi2)) / (1. + (-0.9919475346968611 + 0.367620218664352*eta + 4.274567337924067*eta2)*totchi)) + (- 0.09803730445895877*dchi*Seta*(1. - 3.2283713377939134*eta)*eta2 + 0.01118530335431078*dchi*dchi*eta2*eta - 0.01978238971523653*dchi*Seta*(1. - 4.91667749015812*eta)*eta*totchi))
    fring = ((0.05947169566573468 - 0.14989771215394762*aeff + 0.09535606290986028*aeff*aeff + 0.02260924869042963*aeff*aeff*aeff - 0.02501704155363241*aeff*aeff*aeff*aeff - 0.005852438240997211*(aeff^5) + 0.0027489038393367993*(aeff^6) + 0.0005821983163192694*(aeff^7))/(1 - 2.8570126619966296*aeff + 2.373335413978394*aeff*aeff - 0.6036964688511505*aeff*aeff*aeff*aeff + 0.0873798215084077*(aeff^6)))/(1. - Erad)
    fdamp = ((0.014158792290965177 - 0.036989395871554566*aeff + 0.026822526296575368*aeff*aeff + 0.0008490933750566702*aeff*aeff*aeff - 0.004843996907020524*aeff*aeff*aeff*aeff - 0.00014745235759327472*(aeff^5) + 0.0001504546201236794*(aeff^6))/(1 - 2.5900842798681376*aeff + 1.8952576220623967*aeff*aeff - 0.31416610693042507*aeff*aeff*aeff*aeff + 0.009002719412204133*(aeff^6)))/(1. - Erad)

    Z1tmp = 1. + cbrt((1. - aeff*aeff) ) * (cbrt(1. + aeff) + cbrt(1. - aeff))
    Z1tmp = ifelse(Z1tmp>3., 3., Z1tmp)
    Z2tmp = sqrt(3. *aeff*aeff + Z1tmp*Z1tmp)
    fISCO  = (1. / ((3. + Z2tmp - sign(aeff)*sqrt((3. - Z1tmp) * (3. + Z1tmp + 2. *Z2tmp)))^(3. /2.) + aeff))/pi
    fMECO    = (((0.018744340279608845 + 0.0077903147004616865*eta + 0.003940354686136861*eta2 - 0.00006693930988501673*eta2*eta)/(1. - 0.10423384680638834*eta)) + ((chiPN*(0.00027180386951683135 - 0.00002585252361022052*chiPN + eta2*eta2*(-0.0006807631931297156 + 0.022386313074011715*chiPN - 0.0230825153005985*chiPN2) + eta2*(0.00036556167661117023 - 0.000010021140796150737*chiPN - 0.00038216081981505285*chiPN2) + eta*(0.00024422562796266645 - 0.00001049013062611254*chiPN - 0.00035182990586857726*chiPN2) + eta2*eta*(-0.0005418851224505745 + 0.000030679548774047616*chiPN + 4.038390455349854e-6*chiPN2) - 0.00007547517256664526*chiPN2))/(0.026666543809890402 + (-0.014590539285641243 - 0.012429476486138982*eta + 1.4861197211952053*eta2*eta2 + 0.025066696514373803*eta2 + 0.005146809717492324*eta2*eta)*chiPN + (-0.0058684526275074025 - 0.02876774751921441*eta - 2.551566872093786*eta2*eta2 - 0.019641378027236502*eta2 - 0.001956646166089053*eta2*eta)*chiPN2 + (0.003507640638496499 + 0.014176504653145768*eta + 1. *eta2*eta2 + 0.012622225233586283*eta2 - 0.00767768214056772*eta2*eta)*chiPN2*chiPN)) + (dchi*dchi*(0.00034375176678815234 + 0.000016343732281057392*eta)*eta2 + dchi*Seta*eta*(0.08064665214195679*eta2 + eta*(-0.028476219509487793 - 0.005746537021035632*chiPN) - 0.0011713735642446144*chiPN)))

    return (aeff = aeff, Erad = Erad, fring = fring, fdamp = fdamp, fMECO = fMECO, fISCO = fISCO, totchi = totchi, dchi = dchi, chiPN = chiPN)
end

# Intermediate phase derivative ansatz of IMRPhenomXAS (Eq. 7.6 of arXiv:2001.11412), phase_int_list = (b0, b1, b2, b3, b4, cL)
function _xe_intermediate_dphase(ff, fda, frd, phase_int_list)
    b0, b1, b2, b3, b4, cL = phase_int_list
    invff1 = 1.0 / ff; invff2 = invff1 * invff1; invff3 = invff2 * invff1; invff4 = invff3 * invff1
    return b0 + b1 * invff1 + b2 * invff2 + b3 * invff3 + b4 * invff4 + (4.0 * cL) / ((4.0 * fda * fda) + (ff - frd) * (ff - frd))
end

# Integrated intermediate phase ansatz of IMRPhenomXAS
function _xe_intermediate_phase(ff, frd, fda, phase_int_list)
    b0, b1, b2, b3, b4, cL = phase_int_list
    invff1 = 1.0 / ff; invff2 = invff1 * invff1; invff3 = invff1 * invff2
    return b0 * ff + b1 * log(ff) - b2 * invff1 - b3 * invff2 / 2.0 - (b4 * invff3 / 3.0) + (2.0 * cL * atan((ff - frd) / (2.0 * fda))) / fda
end

# Ringdown phase derivative ansatz of IMRPhenomXAS (Eq. 7.11 of arXiv:2001.11412)
function _xe_ringdown_dphase(ff, fda, frd, c0, c1, c2, c4, cL)
    invf2 = ff^(-2); invf4 = invf2 * invf2; invf1o3 = ff^(-1.0 / 3.0)
    return c0 + c1 * invf1o3 + c2 * invf2 + c4 * invf4 + (cL / (fda * fda + (ff - frd) * (ff - frd)))
end

# Integrated ringdown phase ansatz of IMRPhenomXAS
function _xe_ringdown_phase(ff, fda, frd, c0, c1, c2, c4ov3, cLovfda)
    invf = 1.0 / ff; invf3 = invf * invf * invf; f2o3 = ff^(2.0 / 3.0)
    return c0 * ff + 1.5 * c1 * f2o3 - c2 * invf - c4ov3 * invf3 + (cLovfda * atan((ff - frd) / fda))
end

##############################################################################
#   Eccentric amplitudes
##############################################################################

# 3PN time-domain amplitude of the j-th mean-anomaly harmonic of the 22 mode, evaluated along the dynamics
# (spa_xe_num.compute_td_j_harmonic). Note the convention: harmonic j is obtained with the coefficient index -j.
function _xe_td_j_harmonic(et::AbstractVector, v::AbstractVector, eta, delta, chiA, chiS, jj::Integer)
    T = promote_type(eltype(et), eltype(v), typeof(eta), typeof(chiA), typeof(chiS))
    out = Vector{Complex{T}}(undef, length(et))
    fac = 8.0 * eta * sqrt(pi / 5.0)
    for i in eachindex(et)
        v1 = v[i]; v2 = v1 * v1; v3 = v2 * v1; v4 = v3 * v1; v5 = v4 * v1; v6 = v5 * v1
        logx = 2.0 * log(v1)
        B0, B1, B15, B2, B25, B3 = _amph22_TD_3PN_e12(et[i], logx, eta, delta, chiA, chiS, jj)
        out[i] = (B0 + v2 * B1 + v3 * B15 + v4 * B2 + v5 * B25 + v6 * B3) * fac * v2
    end
    return out
end

# Inspiral-merger IMRPhenomT (2,2) amplitude along the dynamics (compute_insp_merg_PhenomT_amplitude_opt)
function _xe_insp_merg_PhenomT_amplitude(t_insp::AbstractVector, x_spa_qc::AbstractVector, xt_cut, phT)
    T = promote_type(eltype(t_insp), eltype(x_spa_qc), typeof(phT.fac0))
    n = length(t_insp)
    amp = Vector{Complex{T}}(undef, n)
    ninsp = x_spa_qc[1] < xt_cut ? searchsortedfirst(x_spa_qc, xt_cut) - 1 : 0
    for i in 1:ninsp
        amp[i] = phT.fac0 * x_spa_qc[i] * _phT_inspiral_amp_bare(x_spa_qc[i], phT.pn_real, phT.pn_imag, phT.pseudo)
    end
    for i in (ninsp + 1):n
        amp[i] = _phT_intermediate_amp(t_insp[i], phT.alpha1RD, phT.mergerC1, phT.mergerC2, phT.mergerC3, phT.mergerC4, phT.tshift)
    end
    return amp
end

# Eccentric corrections to the SPA amplitude of the j = 0 harmonic (compute_amp_ecc_corrections_j0)
function _xe_amp_ecc_corrections_j0(eta, delta, chiA, chiS, et::AbstractVector, v_spa::AbstractVector, om_dot::AbstractVector,
                                    om_dot_qc::AbstractVector, amp22_qc::AbstractVector, mm)
    amp_j_td_ecc = _xe_td_j_harmonic(et, v_spa, eta, delta, chiA, chiS, 0)
    amp_j_td_e0 = _xe_td_j_harmonic(zero(et), v_spa, eta, delta, chiA, chiS, 0)
    amp_ecc = similar(amp_j_td_ecc)
    for i in eachindex(amp_ecc)
        amp_j0_td_imr = amp22_qc[i] + (amp_j_td_ecc[i] - amp_j_td_e0[i])
        amp_j0_spa_imr = amp_j0_td_imr * sqrt(2.0 * pi / abs(mm * om_dot[i]))
        amp_j0_spa_e0_imr = amp22_qc[i] * sqrt(2.0 * pi / abs(mm * om_dot_qc[i]))
        amp_ecc[i] = amp_j0_spa_imr - amp_j0_spa_e0_imr
    end
    return amp_ecc
end

# Add the windowed eccentric amplitude corrections to the IMRPhenomXAS amplitude (compute_windowed_ecc_amplitude_corrections_v5)
function _xe_windowed_amplitude!(amp_xe::AbstractVector, Mfs_xas::AbstractVector, Mfs_spa::AbstractVector, f_spa_cut,
                                 amp_xas::AbstractVector, amp_ecc::AbstractVector, beta)
    n = length(Mfs_xas)
    # number of grid points below the SPA start. The SPA starts at the minimum frequency of the grid up to roundoff,
    # so a small relative tolerance avoids zeroing the first bin because of a one-ulp difference
    i_spa_start = searchsortedfirst(Mfs_xas, Mfs_spa[1] * (1.0 - 1e-10)) - 1
    i_spa_cut = searchsortedfirst(Mfs_xas, f_spa_cut) - 1      # number of grid points below the SPA cut
    for idx in 1:i_spa_start
        amp_xe[idx] = zero(eltype(amp_xe))
    end
    if i_spa_cut <= i_spa_start
        for idx in (i_spa_start + 1):n
            amp_xe[idx] = amp_xas[idx]
        end
        return amp_xe
    end
    c, _ = _notaknot_spline(Mfs_spa, amp_ecc)
    segment = 1
    nseg = length(Mfs_spa) - 1
    last_correction = zero(eltype(amp_ecc))
    for idx in (i_spa_start + 1):i_spa_cut
        Mf = Mfs_xas[idx]
        while segment < nseg && Mf >= Mfs_spa[segment + 1]
            segment += 1
        end
        correction = _spline_eval_segment(c, segment, Mf - Mfs_spa[segment])
        last_correction = correction
        window = 1.0 / (1.0 + exp(beta * (Mf - f_spa_cut)))
        amp_xe[idx] = amp_xas[idx] + correction * window
    end
    for idx in (i_spa_cut + 1):n
        Mf = Mfs_xas[idx]
        window = 1.0 / (1.0 + exp(beta * (Mf - f_spa_cut)))
        amp_xe[idx] = amp_xas[idx] + last_correction * window
    end
    return amp_xe
end

##############################################################################
#   j != 0 mean-anomaly harmonics (numerical SPA)
##############################################################################

"""
Sum of the j != 0 eccentric harmonics of the 22 mode on the grid `Mfs_ecc` (compute_3pn_ecc_h22_noj0_opt_v2 with
exact spline evaluation on an arbitrary grid). Each harmonic is evaluated with the SPA along the dynamics, split into
absolute value and phase, interpolated with not-a-knot cubic splines and summed as |A| exp(-i phase).
"""
function _xe_eccentric_harmonics(eta, delta, chiA, chiS, tc, n_harmonics, spa, Mfs_ecc::AbstractVector, mm)
    T = promote_type(eltype(Mfs_ecc), eltype(spa.xt), typeof(eta))
    hj_sum = zeros(Complex{T}, length(Mfs_ecc))
    xt = spa.xt; et = spa.et; tt = spa.tt; lt = spa.lt; lamt = spa.lamt; om_dot = spa.omega_dot; ndot = spa.ndot
    nn_norm = [_compute_dldtnorm_3pn(xt[i], et[i], eta, delta, chiS, chiA) for i in eachindex(xt)]
    om = xt .^ 1.5
    v_spa = sqrt.(xt)
    for jj in -n_harmonics:n_harmonics
        jj == 0 && continue
        # SPA amplitude of the harmonic (the - sign in the coefficient index follows the phenomxpy convention)
        amp_td = _xe_td_j_harmonic(et, v_spa, eta, delta, chiA, chiS, -jj)
        jp = 0.5 .* (mm .+ nn_norm .* jj)
        Mfs_j = om .* jp ./ pi
        phase_j = similar(Mfs_j)
        amp_j = similar(Mfs_j)
        for i in eachindex(Mfs_j)
            denSPA = abs(mm * om_dot[i] + jj * ndot[i])
            a = sqrt(2.0 * pi / denSPA) * amp_td[i]
            # absolute value and argument, guarded against an exactly vanishing amplitude (e.g. where the
            # eccentricity has been clipped to zero), whose derivatives would otherwise be NaN with ForwardDiff
            if _val(real(a)) == 0.0 && _val(imag(a)) == 0.0
                amp_j[i] = zero(real(a))
                arg_a = zero(real(a))
            else
                amp_j[i] = abs(a)
                arg_a = angle(a)
            end
            phase_j[i] = 2.0 * pi * Mfs_j[i] * (tt[i] - tc) - (mm * lamt[i] + jj * lt[i] + pi / 4.0) - arg_a
        end
        # keep the part of the dynamics where the SPA frequency is positive
        idx = findall(jp .> 0)
        isempty(idx) && continue
        Mfs_j = Mfs_j[idx]; phase_j = phase_j[idx]; amp_j = amp_j[idx]
        # cut at the first decreasing frequency so that the interpolation is well defined
        idiff = findfirst(<(0), diff(Mfs_j))
        if idiff !== nothing
            Mfs_j = Mfs_j[1:idiff]; phase_j = phase_j[1:idiff]; amp_j = amp_j[1:idiff]
        end
        length(Mfs_j) > 2 || continue
        # interpolate and accumulate
        camp, _ = _notaknot_spline(Mfs_j, amp_j)
        cph, _ = _notaknot_spline(Mfs_j, phase_j)
        i0 = searchsortedfirst(Mfs_ecc, Mfs_j[1])
        i1 = searchsortedlast(Mfs_ecc, Mfs_j[end])
        i1 < i0 && continue
        segment = 1
        nseg = length(Mfs_j) - 1
        for k in i0:i1
            Mf = Mfs_ecc[k]
            while segment < nseg && Mf >= Mfs_j[segment + 1]
                segment += 1
            end
            u = Mf - Mfs_j[segment]
            amplitude = _spline_eval_segment(camp, segment, u)
            phase = _spline_eval_segment(cph, segment, u)
            hj_sum[k] += amplitude * (cos(phase) - 1im * sin(phase))
        end
    end
    return hj_sum
end

##############################################################################
#   GWJulia interface
##############################################################################

"""
Absolute value of the plus and cross polarizations of IMRPhenomXE, see `hphc(model::PhenomXE, ...)`.
"""
function PolAbs(model::PhenomXE, f::AbstractVector, mc, eta, chi1, chi2, dL, iota, ecc, meanAno; kwargs...)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno; kwargs...)
    return [abs.(hp), abs.(hc)]
end

"""
Plus and cross polarizations of IMRPhenomXE with their complete phase, see `hphc(model::PhenomXE, ...)`.
ATTENTION: as for PhenomXHM, the complete phase is included here and `Phi(model::PhenomXE, ...)` returns zero.
"""
function Pol(model::PhenomXE, f::AbstractVector, mc, eta, chi1, chi2, dL, iota, ecc, meanAno; kwargs...)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno; kwargs...)
    return [hp, hc]
end

"""
The phase of IMRPhenomXE is already included in `Pol`, this returns zero (same convention as PhenomXHM).
"""
function Phi(model::PhenomXE, f, mc, eta, chi1, chi2, optional_param...; kwargs...)
    return f .* 0.0
end
