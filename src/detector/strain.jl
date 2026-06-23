
"""
The function computes the pattern functions of the detector for a GW event, they are function of the detector position at the time of arrival of the GW, the sky position of the source, the polarisation of the GW.

    _patternFunction(theta, phi, psi, tRef, DetectorCoordinates)

#### Input arguments:
-  `theta` : float, in rad
-  `phi`   : float, in rad
-  `psi`   : float, in rad
-  `tRef`  : float, in LMST, it is the time at which the GW arrives at the center of the Earth, it is equal to `tcoal` if we are not considering the Earth motion during the measurement
-  `DetectorCoordinates` : structure, containing the coordinates of the detector

#### Optional arguments:
-  `alpha_grad` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry

#### Output:
- `Fp`  : float, plus pattern function
- `Fc`  : float, cross pattern function

#### Example:
```julia
Fp, Fc = _patternFunction(0.1, 0.2, 0.3, 0.4, CE1Id_coordinates)
```

ToDo: Not up to date
"""
function _patternFunction(
    model::Model,
    DetectorCoordinates::DetectorStructure,
    tRef,
    theta,
    phi,
    psi;
    alpha_grad = 0.0,
)

    alpha_rad = alpha_grad * pi / 180.0
    lambda_rad = DetectorCoordinates.latitude_rad                        # latitude of detector site
    zeta_rad = DetectorCoordinates.arm_aperture_rad                      # angle between intereferrometer arms
    phi_mv_rad = @. DetectorCoordinates.longitude_rad + 2.0 * pi * tRef  # angles between the local meridian and the vernal point
    gamma_rad = DetectorCoordinates.orientation_rad + alpha_rad          # orientation of the detector’s arms with respec to local geographical directions, s measured counterclockwise from East to the bisector of the interferometer arms
    
    # detector tensor components in the celestial sphere frame
    cd = -cos(2*gamma_rad - zeta_rad)/4 + cos(2*gamma_rad + zeta_rad)/4
    sd = -sin(2*gamma_rad - zeta_rad)/4 + sin(2*gamma_rad + zeta_rad)/4   
    dt_cs_1 = @. cd/4*(cos(2*phi_mv_rad)*(3-cos(2*lambda_rad))-cos(2*lambda_rad)-1) - sd*sin(2*phi_mv_rad)*sin(lambda_rad)
    dt_cs_2 = @. cd/4*(3-cos(2*lambda_rad))*sin(2*phi_mv_rad) + sd*sin(lambda_rad)*cos(2*phi_mv_rad)
    dt_cs_3 = @. -cd/2*(sin(2*lambda_rad)*cos(phi_mv_rad)) + sd*sin(phi_mv_rad)*cos(lambda_rad)
    dt_cs_4 = @. cd/4*(cos(2*phi_mv_rad)*(cos(2*lambda_rad)-3)-cos(2*lambda_rad)-1) + sd*sin(2*phi_mv_rad)*sin(lambda_rad)
    dt_cs_5 = @. -cd/2*sin(phi_mv_rad)*sin(2*lambda_rad) - sd*cos(phi_mv_rad)*cos(lambda_rad)

    # calculate transformation matrix from the GW frame to the celestial sphere frame
    ra, dec = uc._ra_dec_from_theta_phi_rad(theta, phi)
    mat_1 = [sin(ra)*cos(psi)-cos(ra)*sin(dec)*sin(psi)   -sin(ra)*sin(psi)-cos(ra)*sin(dec)*cos(psi)   -cos(ra)*cos(dec);
             -sin(ra)*sin(dec)*sin(psi)-cos(ra)*cos(psi)  -sin(ra)*sin(dec)*cos(psi)+sin(psi)*cos(ra)  -sin(ra)*cos(dec);
             sin(psi)*cos(dec)                            cos(dec)*cos(psi)                            -sin(dec)]
    
    # contract polarization tensor with detector tensor in celestial sphere frame
    pol_list = _list_polarizations(model)
    pattern_functions = []
    for pol_name in pol_list
        # poltensor in GW frame
        pol_tensor = uc.polarization_dict[pol_name]
        # transfrom to celestial sphere frame
        pol_tensor = mat_1*pol_tensor*transpose(mat_1)

        # contract detector and polarization tensor

        fp = pol_tensor[1,1] .* dt_cs_1 .+ 2 .*pol_tensor[1,2] .* dt_cs_2 .+ 2 .*pol_tensor[1,3] .* dt_cs_3 .+ pol_tensor[2,2] .* dt_cs_4 .+ 2 .*pol_tensor[2,3] .* dt_cs_5 .- pol_tensor[3,3] .* (dt_cs_1 + dt_cs_4)
             
        push!(pattern_functions, fp)
    end

    return pattern_functions
end

"""
Compute the amplitude of the GW signal projected on the detector tensor for each polarization, given a precomputed list of polarizations, such as outputed by Pol(model, DetectorCoordinates, ... ).

    PolarizationDet(model, DetectorCoordinates, pol, f, mc, eta, dL, theta, phi, iota, psi, tcoal, alpha = 0., useEarthMotion = false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `pol`: 
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days

    #### Optional arguments:
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular  geometry.
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement

    #### Output:
    - Returns a list of vectors. Each vector corresponds to a contribution of a polarization of the signal, projected onto the detector.

    Example: 
        ToDo: need example
"""
function PolarizationDet(model::CBC,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    f::AbstractArray,
    mc,
    eta,
    theta,
    phi,
    psi,
    tcoal;
    alpha = 0.0,
    useEarthMotion = false
)

    if useEarthMotion   # technique from GWFAST
        tcoalRescaled = tcoal .- waveform._tau_star(model, f, mc, eta) ./ (3600.0 * 24.0)
        tRef =
            tcoalRescaled .+
            _deltLoc(theta, phi, tcoalRescaled, DetectorCoordinates) ./ (3600.0 * 24.0)
    else
        tRef = tcoal + _deltLoc(theta, phi, tcoal, DetectorCoordinates) / (3600.0 * 24.0)
    end

    #Fp, Fc = _patternFunction(theta, phi, psi, tRef, DetectorCoordinates, alpha_grad=alpha)
    #Ap = @. Fp .* pol[1]
    #Ac = @. Fc .* pol[2]

    # project polarizations onto detector
    Fp = _patternFunction(model, DetectorCoordinates, tRef, theta, phi, psi, alpha_grad=alpha)

    n_fp = length(Fp)
    pol_det = Vector{typeof(pol[1])}(undef, n_fp)
    for idx in 1:n_fp
        pol_det[idx] = @. Fp[idx] .* pol[idx] 
    end

    return pol_det

end

"""
Need documentation 
"""
function PolarizationDet(model::GrModel,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    alpha = 0.0,
    useEarthMotion = false
)

    pol = Pol(
        model,
        f,
        mc,
        eta, 
        chi1,
        chi2,
        dL,
        iota,
        Lambda1,
        Lambda2
    )

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        pol,
        f,
        mc,
        eta,
        theta,
        phi,
        psi,
        tcoal;
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    return pol_det

end

"""
Need documentation 
"""
function PolarizationDet(model::CBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    optional_param... ;
    alpha = 0.0,
    useEarthMotion = false
)

    pol = Pol(
        model,
        f,
        mc,
        eta, 
        chi1,
        chi2,
        dL,
        iota,
        optional_param... 
    )

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        pol,
        f,
        mc,
        eta,
        theta,
        phi,
        psi,
        tcoal;
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    return pol_det

end

function PolarizationDet(model::NonCBC,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    f::AbstractArray,
    theta,
    phi,
    psi,
    tcoal;
    alpha = 0.0,
    useEarthMotion = false
)

    if useEarthMotion   # technique from GWFAST
        tcoalRescaled = tcoal #.- waveform._tau_star(model, f, intrinsic_param_tuple...) ./ (3600.0 * 24.0)
        tRef =
            tcoalRescaled .+
            _deltLoc(theta, phi, tcoalRescaled, DetectorCoordinates) ./ (3600.0 * 24.0)
    else
        tRef = tcoal + _deltLoc(theta, phi, tcoal, DetectorCoordinates) / (3600.0 * 24.0)
    end

    #Fp, Fc = _patternFunction(theta, phi, psi, tRef, DetectorCoordinates, alpha_grad=alpha)
    #Ap = @. Fp .* pol[1]
    #Ac = @. Fc .* pol[2]

    # project polarizations onto detector
    Fp = _patternFunction(model, DetectorCoordinates, tRef, theta, phi, psi, alpha_grad=alpha)

    n_fp = length(Fp)
    pol_det = Vector{typeof(pol[1])}(undef, n_fp)
    for idx in 1:n_fp
        pol_det[idx] = @. Fp[idx] .* pol[idx] 
    end

    return pol_det

end

function PolarizationDet(model::NonCBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    intrinsic_param_tuple, #::Tuple,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    optional_param... ;
    alpha = 0.0,
    useEarthMotion = false
)

    pol = Pol(
        model,
        f,
        intrinsic_param_tuple...,
        dL,
        iota,
        optional_param... 
    )

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        pol,
        f,
        #intrinsic_param_tuple,
        theta,
        phi,
        psi,
        tcoal;
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    return pol_det

end

"""
ToDo: New documentation for this function
This function computes the phase of the waveform seen by the detector, given a waveform model. It already includes the phase due to the Earth motion.

    PhaseDet(model, DetectorCoordinates, f, mc, eta, chi1, chi2, theta, phi, tcoal, phiCoal, DetectorCoordinates, Lambda1, Lambda2; useEarthMotion = false, phase_precomputation = nothing)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `phase_precomputation` : array, default nothing, precomputed phase

    #### Output:
    - `phase`  : array, phase of the GW signal, radians

    #### Example:
    ```julia
    phase = PhaseDet(PhenomD(), CE1Id_coordinates, 1:100, 10.0, 0.25, 0.5, 0.5, 0.1, 0.2, 0.3, 0.4)
    ```

"""
function PhaseDet(
    model::CBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    theta,
    phi,
    tcoal,
    phiCoal;
    useEarthMotion = false,
)

    # Phase of the GW signal

    # Full GW strain expression (complex)
    if useEarthMotion # technique from GWFAST
        tcoalRescaled = tcoal .- waveform._tau_star(model, f, mc, eta) ./ (3600.0 * 24.0)  #tcoal is in fraction of days, tau_star in seconds

    else
        tcoalRescaled = tcoal
    end

    phiD = (2.0 * pi .* f) .* _deltLoc(theta, phi, tcoalRescaled, DetectorCoordinates) # phase due to Earth motion

    return @. 2.0 * pi * (tcoal * 3600.0 * 24.0) .* f .- phiCoal  .+ phiD #.- Phi


end

function PhaseDet(
    model::NonCBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    theta,
    phi,
    tcoal,
    phiCoal;
    useEarthMotion = false,
)


    phiD = (2.0 * pi .* f) .* _deltLoc(theta, phi, tcoal, DetectorCoordinates) # phase due to Earth motion

    return @. 2.0 * pi * (tcoal * 3600.0 * 24.0) .* f .- phiCoal  .+ phiD #.- Phi


end

"""
ToDo: Documentation
"""
function Strain(model::CBC,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    phase::AbstractArray,
    f::AbstractArray,
    mc,
    eta,
    theta,
    phi,
    tcoal,
    phiCoal;
    useEarthMotion = false
)

    phase_det = PhaseDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        theta,
        phi,
        tcoal,
        phiCoal,
        useEarthMotion = useEarthMotion
    )
        
    return sum(pol) .* exp.(1im .* (phase_det .- phase))

end

function Strain(model::NonCBC,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    phase::AbstractArray,
    f::AbstractArray,
    theta,
    phi,
    tcoal,
    phiCoal;
    useEarthMotion = false
)

    phase_det = PhaseDet(
        model,
        DetectorCoordinates,
        f,
        theta,
        phi,
        tcoal,
        phiCoal,
        useEarthMotion = useEarthMotion
    )
        
    return sum(pol) .* exp.(1im .* (phase_det .- phase))

end

"""
ToDo: Documentation
This function computes the full strain (complex) as a function of the parameters, at given frequencies, at detector location, as measured by the detector, since it include the pattern functions, via PolarizationDet and PhaseDet.

    Strain(model, DetectorCoordinates,  f, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2, useEarthMotion=false, alpha=0.)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `ampl_precomputation` : array, default nothing, precomputed amplitude
    -  `phase_precomputation` : array, default nothing, precomputed phase

    #### Output:
    - `strain`  : array, complex strain as seen by the detector

    #### Example:
    ```julia
    strain = Strain(PhenomD(), CE1Id_coordinates, 1:100, 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
    ```
"""
function Strain(model::GrModel,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    useEarthMotion = false,
    alpha = 0.0
)

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        Lambda1,
        Lambda2,
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    phase_wave = Phi(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        Lambda1,
        Lambda2,
    )

    strain_det = Strain(
        model,
        DetectorCoordinates,
        pol_det,
        phase_wave,
        f,
        mc,
        eta,
        theta,
        phi,
        tcoal,
        phiCoal,
        useEarthMotion = useEarthMotion
    )
        
    return strain_det

end

function Strain(model::CBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    optional_param...;
    useEarthMotion = false,
    alpha = 0.0
)

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        optional_param...,
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    phase_wave = Phi(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        optional_param...,
    )

    strain_det = Strain(
        model,
        DetectorCoordinates,
        pol_det,
        phase_wave,
        f,
        mc,
        eta,
        theta,
        phi,
        tcoal,
        phiCoal,
        useEarthMotion = useEarthMotion
    )
        
    return strain_det

end


function Strain(model::NonCBC,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    intrinsic_param,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    optional_param...;
    useEarthMotion = false,
    alpha = 0.0
)

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        f,
        intrinsic_param,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        optional_param...,
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    phase_wave = Phi(
        model,
        f,
        intrinsic_param...,
        iota,
        optional_param...,
    )

    strain_det = Strain(
        model,
        DetectorCoordinates,
        pol_det,
        phase_wave,
        f,
        theta,
        phi,
        tcoal,
        phiCoal,
        useEarthMotion = useEarthMotion
    )
        
    return strain_det

end