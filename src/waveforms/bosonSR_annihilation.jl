include("../../bosonSR_repo/Annihilation.jl")

function PolAbs(model::BosonSR_ann,
    f::AbstractVector,
    p1, # mass ratio
    p2, # alpha
    dL,
    iota;
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        p1,
        p2,
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
    p1,
    p2,
    dL,
    iota
)

    hp, hc = PolAbs(
        model,
        f,
        p1,
        p2,
        dL,
        iota
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

function Phi(model::BosonSR_ann,
    f::AbstractVector,
    p1, # mass ratio
    p2; # alpha
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # is omega in Hz?
    omega = 2pi * f
    M = p1
    mua_ann = p2 # in GeV
    alpha_ann = G*M_ann_*mua_ann*M_sun
    n = 4
    m = l = n-1
    r = 1
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




