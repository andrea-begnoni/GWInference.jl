include("../../bosonSR_repo_/chris/julia_wf/annihilation.jl")
include("../../bosonSR_repo_/chris/julia_wf/level_transition.jl")


function PolAbs(model::BosonSR_ann,
    f::AbstractVector,
    M_solar_ann, 
    mua_ann, 
    dL,
    iota;
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        M_solar_ann,
        mua_ann,
        dL
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
function Pol(model::BosonSR_ann,
    f::AbstractVector,
    M_solar_ann, 
    mua_ann,
    dL,
    iota
)

    hp, hc = PolAbs(
        model,
        f,
        M_solar_ann,
        mua_ann,
        dL,
        iota
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

function Phi(model::BosonSR_ann,
    f::AbstractVector,
    M_solar_ann, 
    mua_ann,
    GMsun_over_c3 = uc.GMsun_over_c3, # not used in this function but included for consistency
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc, # not used in this function but included for consistency
)

    # is omega in Hz?
    omega = 2pi * f

    n = 4
    dL_kpc = dL * 1e6

    phase = angle.(h_ann(omega, M, n, l, alpha, r))


    return phase

end


function Ampl(
    model::BosonSR_ann,
    f::AbstractVector,
    p1, # mass BH
    p2, # mua_ann
    dL
)

    # is omega in Hz?
    omega = f
    M = p1
    mua_ann = p2 # in GeV
    alpha_ann = G*M_ann_*mua_ann*M_sun
    n = 4
    m = l = n-1
    # r in kpc, dL in Gpc
    r = dL * 1e6
    ampl = abs.(h_ann(omega, M, n, l, alpha, r))


    return ampl
    

end


function _fcut(model::BosonSR_ann,
    p1, # mass ratio
    p2, # alpha
)    # dummy cutoff frequency
    # to be replaced with a real one
    return 10_000 # in Hz
end




