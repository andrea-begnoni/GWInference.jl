function PolAbs(model::BosonSR,
    f::AbstractVector,
    p1, # mass ratio
    p2, # alpha
    p3, # primary spin
    dL,
    iota;
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        p1,
        p2,
        p3,
        dL,
        iota;
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
function Pol(model::BosonSR,
    f::AbstractVector,
    p1,
    p2,
    p3,
    iota
)

    hp, hc = PolAbs(
        model,
        f,
        p1,
        p2,
        p3,
        dL,
        iota
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

function Phi(
    model::BosonSR,
    f::AbstractVector,
    p1, # mass ratio
    p2, # alpha
    p3, # primary spin
    iota;
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # Unpack parameters
    q = p1
    alpha = p2

    # dummy phase model
    # to be replaced with a real one
    phi = -2 * (f .^(-5/3) .+ 0.1 * alpha .* q ./ (1 .+ q).^2 .* f .^(-11/3))


    return phi

end


function Ampl(
    model::BosonSR,
    f::AbstractVector,
    p1, # mass ratio
    p2, # alpha
    p3, # primary spin
    dL,
    iota;
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # Unpack parameters
    q = p1
    alpha = p2

    # dummy amplitude model
    # to be replaced with a real one
    ampl = f .^(-7/6) .* (1 .+ 0.5 * alpha .* q ./ (1 .+ q).^2 .* f .^(-8/3))


    return ampl
    

end


function _fcut(model::BosonSR,
    p1, # mass ratio
    p2, # alpha
    p3, # primary spin
)    # dummy cutoff frequency
    # to be replaced with a real one
    return 2048.0
end