
"""
This function computes the Fisher Matrix for a single detector, given the parameters of the event and the detector.
To do the computation it uses the function FisherMatrix_internal(...) for L-shaped detectors and FisherMatrix_Tdetector(...) for T-shaped detectors.
Thus it is a wrapper function that calls the correct function depending on the shape of the detector. More information on the Fisher Matrix computation can be found in the documentation of FisherMatrix_internal(...) and FisherMatrix_Tdetector(...).

    FisherMatrix(model, detector , mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1=0.0, Lambda2=0.0, res=1000, useEarthMotion=false, alpha=0.0, SNR_thres=12., fmin=2., fmax=nothing, coordinate_shift=true, return_SNR=false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `detector` : structure, containing the detector information
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
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `SNR_thres` : float, default 12., SNR threshold for the computation of the Fisher Matrix
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `coordinate_shift` : bool, default true, valid for T detectors, if true the codes shifts the coordinates of the detector from the center of the triangle to the center of the arms (more realistic scenario, recommended)
    -  `return_SNR` : bool, default false, if true the function returns the SNR of the event (skipping the need to call the SNR function)

    #### Output:
    - `FisherMatrix`  : matrix, Fisher Matrix

    #### Example:
    ```julia
    FisherMatrix = FisherMatrix(PhenomD(), CE1Id , 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    ```


"""
function FisherMatrix(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    optional_param...;
    res = 1000,
    useEarthMotion::Bool = false,
    alpha = 0.0,
    SNR_thres::Union{Nothing, Float64} = 12.,
    fmin::Float64=2.,
    fmax::Union{Nothing, Float64}=nothing,
    coordinate_shift::Bool = true,
    return_SNR::Bool = false,
    optimization::Bool = false, # not supported for 1 detector
    call_number = 1 # not supported for 1 detector
)
    #function that is used only to divide between L and T detectors

    if detector.shape =='L'

        return FisherMatrix_internal(
            model,
            detector,
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
            optional_param...,
            SNR_thres=SNR_thres,
            res = res,
            useEarthMotion = useEarthMotion,
            alpha = alpha,
            fmin=fmin,
            fmax=fmax,
            return_SNR = return_SNR,
        )
    elseif detector.shape =='T'

        return FisherMatrix_Tdetector(
            model,
            detector,
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
            optional_param...,
            SNR_thres=SNR_thres,
            res = res,
            useEarthMotion = useEarthMotion,
            alpha = alpha,
            fmin=fmin,
            fmax=fmax,
            coordinate_shift=coordinate_shift,
            return_SNR = return_SNR,
        )
    else
        error("detector shape not recognized")
    end
end

"""
This is the function that computes the Fisher Matrix for a single detector with L shape.
It computes the derivatives of the strain w.r.t. each parameter and then computes the Fisher Matrix. It is called by the function FisherMatrix and FisherMatrix_Tdetector.
There is no need to call this function in your computations since all the logic is implemented in the FisherMatrix function. 

"""

function FisherMatrix_internal(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    optional_param...;
    res = 1000,
    useEarthMotion::Bool = false,
    alpha = 0.0,
    SNR_thres::Union{Nothing, Float64} = 12.,
    fmin::Float64=2.,
    fmax::Union{Nothing, Float64}=nothing,
    return_SNR::Bool = false,
)

    #Define/extract tidal diformabilites
    if _event_type(model::Model) == "BBH"
        Lambda1 = 0.
        Lambda2 = 0.
    elseif _event_type(model::Model) == "BNS"
        Lambda1 = optional_param[1]
        Lambda2 = optional_param[2]
    elseif _event_type(model::Model) == "NSBH"
        Lambda1 = optional_param[1]
        Lambda2 = 0.
    else
        #ToDo: Print error
    end

    if model isa TaylorF2
        nPar = _npar(model, Lambda1, Lambda2)
    else
        nPar = _npar(model)
    end

    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))
    psdGrid = linear_interpolation(detector.fNoise, detector.psd, extrapolation_bc = 1.0)(fgrid)  
    
    # compute SNR and procede only if it is above the threshold
    SNRval = nothing
    if SNR_thres !==nothing
        SNRval = SNR(
            model,
            detector,
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
            fmin = fmin,
            fmax = fmax,
            res = res,
            #ampl_precomputation = ampl_precomputation,
        )
        if SNRval < SNR_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        end
    end

    detectorCoordinates = DetectorCoordinates(
        detector.latitude_rad,
        detector.longitude_rad,
        detector.orientation_rad,
        detector.arm_aperture_rad
    )
    
    ###########  Derivatives of the strain w.r.t. each parameter
    strainAutoDiff_real = Matrix{Float64}(undef, res, nPar)
    strainAutoDiff_imag = Matrix{Float64}(undef, res, nPar)
    
    event_parameter = [mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, optional_param...]
    if model isa TaylorF2 && nPar == 12 && Lambda1 == 0. # if Lambda1 = 0. and Lambda2 != 0. we need to add a zero Lambda1
        event_parameter = [event_parameter[1:11]; Lambda2]
        strain_param = x -> (x[1:11]..., Lambda1, x[12])
    else
        event_parameter = event_parameter[1:nPar] # cut off non-required parameter
        strain_param = x -> x
    end
    strainAutoDiff_real = ForwardDiff.jacobian(
        x -> real(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                strain_param(x)... ,
                alpha = alpha,
                useEarthMotion = useEarthMotion
            ),
        ),
        event_parameter,
    )
    strainAutoDiff_imag = ForwardDiff.jacobian(
        x -> imag(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                strain_param(x)... ,
                alpha = alpha,
                useEarthMotion = useEarthMotion,
            ),
        ),
        event_parameter,
    )

    # It can happen that a certain frequency gives a Nan value, in this case we set the derivative to zero,
    # this happens less than one time per event and usually at the end of the frequency grid.
    strainAutoDiff_real[isnan.(strainAutoDiff_real)] .= 0.0
    strainAutoDiff_imag[isnan.(strainAutoDiff_imag)] .= 0.0
    ######### end of derivatives
    jacobian = Matrix{ComplexF64}(undef, nPar, res)
    for ii in 1:nPar
        jacobian[ii, :] = strainAutoDiff_real[:, ii] + 1im * strainAutoDiff_imag[:, ii]
    end
    jacobian[10,:] /= (3600.0 * 24.0)   # Change the units of the tcoal derivative from days to seconds (this improves conditioning)


    # compute the Fisher matrix
    Fisher = Matrix{Float64}(undef, nPar, nPar)

    for alpha = 1:nPar
        for beta = alpha:nPar
            Fisher[alpha, beta] =
                4.0 *
                trapz(fgrid, real(jacobian[alpha, :] .* conj(jacobian[beta, :])) ./ psdGrid)
            Fisher[beta, alpha] = Fisher[alpha, beta]
        end
    end

    if return_SNR == true
        return Fisher, SNRval
    else
        return Fisher
    end

end

"""
Same as FisherMatrix_internal above, but for PhenomHM and PhenomXHM using the waveform, `waveform_values`, and its jacobian w.r.t. (mc, eta, chi1, chi2, dL, iota),
`waveform_jacobian`, precomputed on `fgrid` (see _hphc_values_jacobian). It is called by FisherMatrix when `optimization = true`.
"""
function FisherMatrix_internal(model::Union{PhenomHM, PhenomXHM},
    detector::Detector,
    fgrid::AbstractArray,
    waveform_values::AbstractArray,
    waveform_jacobian::AbstractArray,
    mc::Float64,
    eta::Float64,
    theta::Float64,
    phi::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64;
    alpha = 0.0,
    useEarthMotion::Bool = false,
)
    psdGrid = linear_interpolation(detector.fNoise, detector.psd, extrapolation_bc = 1.0)(fgrid)  

    detectorCoordinates = DetectorCoordinates(
        detector.latitude_rad,
        detector.longitude_rad,
        detector.orientation_rad,
        detector.arm_aperture_rad
    )

    nPar = _npar(model)
    len = length(fgrid)

    ###########  Derivatives of the strain w.r.t. each parameter
    # chi1, chi2, dL and iota enter only through the precomputed waveform jacobian, their values here are not used
    event_parameter = [mc, eta, 0., 0., 0., theta, phi, 0., psi, tcoal, phiCoal]

    strainAutoDiff_real = ForwardDiff.jacobian(
        x -> real(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                waveform_values,
                waveform_jacobian,
                x[1],x[2],x[6],x[7],x[9],x[10],x[11],
                alpha = alpha,
                useEarthMotion = useEarthMotion
            ),
        ),
        event_parameter,
    )

    strainAutoDiff_imag = ForwardDiff.jacobian(
        x -> imag(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                waveform_values,
                waveform_jacobian,
                x[1],x[2],x[6],x[7],x[9],x[10],x[11],
                alpha = alpha,
                useEarthMotion = useEarthMotion
            ),
        ),
        event_parameter,
    )

    # It can happen that a certain frequency gives a Nan value, in this case we set the derivative to zero,
    # this happens less than one time per event and usually at the end of the frequency grid.
    strainAutoDiff_real[isnan.(strainAutoDiff_real)] .= 0.0
    strainAutoDiff_imag[isnan.(strainAutoDiff_imag)] .= 0.0
    ######### end of derivatives
    jacobian = Matrix{ComplexF64}(undef, nPar, len)
    for ii in 1:nPar
        jacobian[ii, :] = strainAutoDiff_real[:, ii] + 1im * strainAutoDiff_imag[:, ii]
    end
    jacobian[10,:] /= (3600.0 * 24.0)   # Change the units of the tcoal derivative from days to seconds (this improves conditioning)


    # compute the Fisher matrix
    Fisher = Matrix{Float64}(undef, nPar, nPar)

    for alpha = 1:nPar
        for beta = alpha:nPar
            Fisher[alpha, beta] =
                4.0 *
                trapz(fgrid, real(jacobian[alpha, :] .* conj(jacobian[beta, :])) ./ psdGrid)
            Fisher[beta, alpha] = Fisher[alpha, beta]
        end
    end

    return Fisher

end

"""
Computes hp and hc of PhenomHM or PhenomXHM on `fgrid`, together with their jacobian w.r.t. (mc, eta, chi1, chi2, dL, iota).
The output is used by FisherMatrix when `optimization = true`, so that the waveform derivatives are computed only once for all the detectors.

    waveform_values, waveform_jacobian = _hphc_values_jacobian(model, fgrid, mc, eta, chi1, chi2, dL, iota; call_number = 1)

    #### Optional arguments:
    -  `call_number` : int, default 1, how the waveform is called inside ForwardDiff:
        1: the values are stored in a container and real and imaginary parts are differentiated in a single call (fastest),
        2: the values are stored in a container and real and imaginary parts are differentiated separately,
        3: standard, the values are computed with a separate call to the waveform.

    #### Output:
    - `waveform_values` : array, [hp; hc]
    - `waveform_jacobian` : matrix, derivatives of [hp; hc], one column per parameter
"""
function _hphc_values_jacobian(model::Union{PhenomHM, PhenomXHM},
    fgrid::AbstractArray,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    iota::Float64;
    call_number = 1,
)
    res = length(fgrid)
    waveform_parameter = [mc, eta, chi1, chi2, dL, iota]

    if call_number == 1 # use of the container and real + imag in single call
        waveform_values = zeros(ComplexF64, 2*res)
        waveform_jacobian_ = ForwardDiff.jacobian( x-> hphc(model, fgrid, x..., container=waveform_values, call_number=1, optimization=true), waveform_parameter)
        waveform_jacobian_hp_real = waveform_jacobian_[1:res,:]
        waveform_jacobian_hp_imag = waveform_jacobian_[res+1:2*res,:]
        waveform_jacobian_hc_real = waveform_jacobian_[2*res+1:3*res,:]
        waveform_jacobian_hc_imag = waveform_jacobian_[3*res+1:4*res,:]
        waveform_jacobian = [waveform_jacobian_hp_real + 1im .* waveform_jacobian_hp_imag ; waveform_jacobian_hc_real + 1im .* waveform_jacobian_hc_imag]

    elseif call_number == 2 # use of the container
        waveform_values = zeros(ComplexF64, 2*res)
        waveform_jacobian_real = ForwardDiff.jacobian( x-> real(hphc(model, fgrid, x..., container=waveform_values, call_number=2, optimization=true)), waveform_parameter)
        waveform_jacobian_imag = ForwardDiff.jacobian( x-> imag(hphc(model, fgrid, x..., call_number=2, optimization=true)), waveform_parameter)
        waveform_jacobian = waveform_jacobian_real + 1im .* waveform_jacobian_imag

    elseif call_number == 3 # standard
        waveform_jacobian_real = ForwardDiff.jacobian( x-> real(hphc(model, fgrid, x..., call_number=3, optimization=true)), waveform_parameter)
        waveform_jacobian_imag = ForwardDiff.jacobian( x-> imag(hphc(model, fgrid, x..., call_number=3, optimization=true)), waveform_parameter)
        waveform_jacobian = waveform_jacobian_real + 1im .* waveform_jacobian_imag
        waveform_values = hphc(model, fgrid, mc, eta, chi1, chi2, dL, iota, call_number=3, optimization=true)

    else
        error("call_number must be 1, 2 or 3")
    end

    return waveform_values, waveform_jacobian
end

"""
This function computes the *Fisher Matrix*, as a function of the parameters of the event, as measured by a NETWORK of detectors.
It relies on the function FisherMatrix(..., detector::Detector, ...), which computes the Fisher for a single detector.
where the dots indicate the parameters equal to the previous function call.

#### Optional arguments:
-  `optimization` : bool, default true, only for PhenomHM and PhenomXHM with a network of detectors, if true the derivatives of the waveform are computed only once and reused for each detector (faster)
-  `call_number` : int, default 1, only used if `optimization` is true, selects how the waveform derivatives are computed (1, 2 or 3, see _hphc_values_jacobian)

#### Example:
```julia
    FisherMatrix(... [CE1Id, CE2NM], ...)
```

"""
function FisherMatrix(model::Model,
    detector::Vector{Detector},
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    optional_param...;
    res = 1000,
    useEarthMotion::Bool = false,
    SNR_thres::Union{Nothing, Float64}=12.,
    alpha = 0.0,
    fmin::Float64=2.0,
    fmax::Union{Nothing, Float64} = nothing,
    coordinate_shift::Bool = true,
    return_SNR::Bool = false,
    optimization::Bool = true,
    call_number = 1,
)

    #Define/extract tidal diformabilites
    if _event_type(model::Model) == "BBH"
        Lambda1 = 0.
        Lambda2 = 0.
    elseif _event_type(model::Model) == "BNS"
        Lambda1 = optional_param[1]
        Lambda2 = optional_param[2]
    elseif _event_type(model::Model) == "NSBH"
        Lambda1 = optional_param[1]
        Lambda2 = 0.
    else
        #ToDo: Print error
    end

    # compute SNR and procede only if it is above the threshold
    if model isa TaylorF2
        nPar = _npar(model, Lambda1, Lambda2)
    else
        nPar = _npar(model)
    end

    if optimization == true
        # check call_number
        if call_number != 1 && call_number != 2 && call_number != 3
            error("call_number must be 1, 2 or 3")
        end
    end

    SNRval = nothing
    if SNR_thres !==nothing || return_SNR == true
        SNRval = SNR(
            model,
            detector,
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
            fmin = fmin,
            fmax = fmax,
            res = res,
            useEarthMotion = useEarthMotion,
        )
        if SNR_thres !==nothing && SNRval < SNR_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        
        end
    end

    # For PhenomHM and PhenomXHM the derivatives of the waveform are computed only once and reused for each detector
    use_optimization = optimization == true && model isa Union{PhenomHM, PhenomXHM}
    if use_optimization
        if isnothing(fmax)
            fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        else
            fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
            fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
        end

        fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))

        waveform_values, waveform_jacobian = _hphc_values_jacobian(model, fgrid, mc, eta, chi1, chi2, dL, iota, call_number = call_number)
    end

    fisherList = Vector{Matrix{Float64}}(undef, length(detector))
    for i in eachindex(detector)
        

        if detector[i].shape == 'L' && use_optimization
            F = FisherMatrix_internal(
                model,
                detector[i],
                fgrid,
                waveform_values,
                waveform_jacobian,
                mc,
                eta,
                theta,
                phi,
                psi,
                tcoal,
                phiCoal,
                alpha = alpha,
                useEarthMotion = useEarthMotion,
            )
        elseif detector[i].shape == 'T' && use_optimization
            F = FisherMatrix_Tdetector(
                model,
                detector[i],
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
                optional_param...,
                SNR_thres=nothing,
                res = res,
                useEarthMotion = useEarthMotion,
                alpha = alpha,
                fmin=fmin,
                fmax=fmax,
                coordinate_shift = coordinate_shift,
                return_SNR=false,
                fgrid = fgrid,
                waveform_values = waveform_values,
                waveform_jacobian = waveform_jacobian,
            )
        elseif detector[i].shape == 'L'
            F = FisherMatrix_internal(
                model,
                detector[i],
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
                optional_param...,
                SNR_thres=nothing,
                res = res,
                useEarthMotion = useEarthMotion,
                alpha = alpha,
                fmin=fmin,
                fmax=fmax,
                return_SNR=false,
            )
        elseif detector[i].shape == 'T'
            F = FisherMatrix_Tdetector(
                model,
                detector[i],
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
                optional_param...,
                SNR_thres=nothing,
                res = res,
                useEarthMotion = useEarthMotion,
                alpha = alpha,
                fmin=fmin,
                fmax=fmax,
                coordinate_shift = coordinate_shift,
                return_SNR=false,
            )
        end
        fisherList[i] = F

    end
    if return_SNR == true
        return sum(fisherList, dims = 1)[1], SNRval
    else
        return sum(fisherList, dims = 1)[1]
    end

end

"""
This is a helper function that computes the Fisher Matrix for a single detector with T shape. It is called by the function FisherMatrix and calls FisherMatrix_internal.
"""
function FisherMatrix_Tdetector(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    optional_param...;
    res = 1000,
    useEarthMotion::Bool = false,
    SNR_thres::Union{Nothing, Float64}=12.,
    alpha = 0.0,
    fmin::Float64=2.0,
    fmax::Union{Nothing, Float64} = nothing,
    REarth_km = uc.REarth_km,
    coordinate_shift::Bool = true,
    return_SNR::Bool = false, 
    fgrid = nothing,
    waveform_values = nothing,
    waveform_jacobian = nothing,
)
    # if fgrid, waveform_values and waveform_jacobian are given (PhenomHM and PhenomXHM with optimization = true, see FisherMatrix)
    # the precomputed waveform derivatives are used for the three arms
    use_optimization = model isa Union{PhenomHM, PhenomXHM} && !isnothing(waveform_jacobian)

    SNRval = nothing

    if SNR_thres !== nothing

        #Define/extract tidal diformabilites
        if _event_type(model::Model) == "BBH"
            Lambda1 = 0.
            Lambda2 = 0.
        elseif _event_type(model::Model) == "BNS"
            Lambda1 = optional_param[1]
            Lambda2 = optional_param[2]
        elseif _event_type(model::Model) == "NSBH"
            Lambda1 = optional_param[1]
            Lambda2 = 0.
        else
            #ToDo: Print error
        end

        if model isa TaylorF2
            nPar = _npar(model, Lambda1, Lambda2)
        else
            nPar = _npar(model)
        end

        SNRval = SNR(
            model,
            detector,
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
            fmin = fmin,
            fmax = fmax,
            res = res,
            #ampl_precomputation = ampl_precomputation
        )
        if SNRval < SNR_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        end
    end

    # We write the T-detector as three detectors in slightly different positions
    # Detector has to be in the plane orthogonal to the radius
    # write coordinates on a tangent plane (https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates)
    # Note the minus sign because theta goes from north to south

    if coordinate_shift == true
        lat = detector.latitude_rad
        long = detector.longitude_rad
        first = [- cos(lat)*cos(long), 
                - cos(lat)*sin(long), 
                sin(lat)]
        #No sign needed here because phi is counterclosk wise (so goes towards east)
        second = [-sin(long), cos(long), 0.0]
        ETarm           = 10e3                      ## ET arms in m
        ET_cartesian   =      REarth_km * 1e3 .* [sin(lat)*cos(long), sin(lat)*sin(long), cos(lat)]  ## ET center in m
        #These are the positions of the 3 detectors
        ET1_cartesian = @. ET_cartesian + ETarm * (-.5 * first - 0.28867513 * second)
        ET2_cartesian = @. ET_cartesian + ETarm * (+.5 * first - 0.28867513 * second)
        ET3_cartesian = @. ET_cartesian + ETarm * (+0.57735027 * second)

        # go back to spherical coordinates and discard radius (it is REarth at 1e-7)

        ET1_coo = [acos(ET1_cartesian[3]/sqrt(sum(ET1_cartesian.^2))), atan(ET1_cartesian[2], ET1_cartesian[1])]
        ET2_coo = [acos(ET2_cartesian[3]/sqrt(sum(ET2_cartesian.^2))), atan(ET2_cartesian[2], ET2_cartesian[1])]
        ET3_coo = [acos(ET3_cartesian[3]/sqrt(sum(ET3_cartesian.^2))), atan(ET3_cartesian[2], ET3_cartesian[1])]

        ET1 = Detector(ET1_coo[1], ET1_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET2 = Detector(ET2_coo[1], ET2_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET3 = Detector(ET3_coo[1], ET3_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
    else
        ET1 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET2 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET3 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
    end

    if use_optimization
        parameters = [mc, eta, theta, phi, psi, tcoal, phiCoal]
        F1 = FisherMatrix_internal(
            model,
            ET1,
            fgrid,
            waveform_values,
            waveform_jacobian,
            parameters...,
            alpha = 0.0,
            useEarthMotion = useEarthMotion,
        )
        F2 = FisherMatrix_internal(
            model,
            ET2,
            fgrid,
            waveform_values,
            waveform_jacobian,
            parameters...,
            alpha = 60.0,
            useEarthMotion = useEarthMotion,
        )
        F3 = FisherMatrix_internal(
            model,
            ET3,
            fgrid,
            waveform_values,
            waveform_jacobian,
            parameters...,
            alpha = 120.0,
            useEarthMotion = useEarthMotion,
        )

    else
        F1 = FisherMatrix_internal(
            model,
            ET1,
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
            optional_param...,
            res = res,
            useEarthMotion = useEarthMotion,
            SNR_thres = nothing,
            alpha = 0.0,
            fmin=fmin,
            fmax=fmax,
            return_SNR=false,
        )
        F2 = FisherMatrix_internal(
            model,
            ET2,
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
            optional_param...,
            res = res,
            useEarthMotion = useEarthMotion,
            SNR_thres = nothing,
            alpha = 60.0,
            fmin=fmin,
            fmax=fmax,
            return_SNR=false,
        )
        F3 = FisherMatrix_internal(
            model,
            ET3,
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
            optional_param...,
            res = res,
            useEarthMotion = useEarthMotion,
            SNR_thres = nothing,
            alpha = 120.0,
            fmin=fmin,
            fmax=fmax,
            return_SNR=false,

        )
    end
    if return_SNR == true
        return F1 + F2 + F3, SNRval
    else
        return F1 + F2 + F3
    end
end

"""
The main function of the code, it computes the Fisher Matrix for an array of events, as measured by a single detector or a network of detectors. It also calculates the 
SNRs if requested and saves the results (Fisher matrices and SNRs) in a file if the optional argument `auto_save` is set to true. The file is saved in the folder `output/name_folder/Fishers_SNRs.h5`


    FisherMatrix(model, detector, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal,  Lambda1=0.0, Lambda2=0.0, res=1000, useEarthMotion=false, alpha=0.0, SNR_thres=12., fmin=2., fmax=nothing, coordinate_shift=true, return_SNR=false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `detector` : structure, containing the detector information
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
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `SNR_thres` : float, default 12., SNR threshold for the computation of the Fisher Matrix
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `coordinate_shift` : bool, default true, valid for T detectors, if true the codes shifts the coordinates of the detector from the center of the triangle to the center of the arms (more realistic scenario, recommended)
    -  `return_SNR` : bool, default false, if true the function returns the SNR of the event (skipping the need to call the SNR function)
    -  `auto_save` : bool, default false, if true the function saves the results in a file
    -  `name_folder` : string, name of the folder where the results are saved, if the default is left, it saves BBH in the folder "output/BBH" and so on for each source type
    -  `optimization` : bool, default true, only for PhenomHM and PhenomXHM with a network of detectors, if true the derivatives of the waveform are computed only once and reused for each detector (faster)
    -  `call_number` : int, default 1, only used if `optimization` is true, selects how the waveform derivatives are computed (1, 2 or 3, see _hphc_values_jacobian)

    #### Output:
    - `FisherMatrix`  : matrix, Fisher Matrix

    #### Example:
    ```julia
    FisherMatrix = FisherMatrix(PhenomD(), [10.0, 20.], [0.25, 0.25], [0.5, 1.], [0.5, -1.], [1.0, 3.], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.6], [0.6, 0.7], CE1Id)
    ```


"""
function FisherMatrix(model::Model,
    detector::Union{Detector, Vector{Detector}},
    mc::AbstractArray,
    eta::AbstractArray,
    chi1::AbstractArray,
    chi2::AbstractArray,
    dL::AbstractArray,
    theta::AbstractArray,
    phi::AbstractArray,
    iota::AbstractArray,
    psi::AbstractArray,
    tcoal::AbstractArray,
    phiCoal::AbstractArray,
    optional_param...;
    fmin::Union{Float64, AbstractArray}=2.0,
    fmax::Union{Nothing, Float64, AbstractArray} = nothing,
    res = 1000,
    useEarthMotion::Bool = false,
    SNR_thres::Union{Nothing, Float64} =12.,
    alpha = 0.0,
    coordinate_shift::Bool = true,
    return_SNR::Bool = false,
    auto_save::Bool =false,
    name_folder = nothing,
    save_catalog::Bool = false,
    optimization::Bool = true,
    call_number = 1,
)
    nEvents = length(mc)    

    if(fmin isa AbstractArray && length(fmin) != nEvents)
        throw(ArgumentError("fmin must be an array of the same length as the number of events (or otherwise a single scalar Float64)"))
    end

    if(fmax isa AbstractArray && length(fmax) != nEvents)
        throw(ArgumentError("fmax must be an array of the same length as the number of events (or otherwise a single scalar Float64 or Nothing)"))
    end

    if name_folder === nothing
        name_folder = _event_type(model) 
    end

    #Define/extract tidal diformabilites
    if _event_type(model::Model) == "BBH"
        Lambda1 = zeros(nEvents)
        Lambda2 = zeros(nEvents)
    elseif _event_type(model::Model) == "BNS"
        Lambda1 = optional_param[1]
        Lambda2 = optional_param[2]
        if name_folder == "BBH"
            name_folder = "BNS"
        end
    elseif _event_type(model::Model) == "NSBH"
        Lambda1 = optional_param[1]
        Lambda2 = zeros(nEvents)
        if name_folder == "BBH"
            name_folder = "NSBH"
        end
    else
        #ToDo: Print error
    end

    if model isa TaylorF2
        nPar = _npar(model, Lambda1[1], Lambda2[1])
    else
        nPar = _npar(model)
    end

    # restructure the optional parameter
    optional_param_reshaped = Array{Vector{Float64}}(undef, nEvents)
    for ii in 1:nEvents  
        optional_param_ii = []
        for op in optional_param
            append!(optional_param_ii, op[ii])
        end
        optional_param_reshaped[ii] = optional_param_ii
    end

    Fishers = Array{Float64}(undef, nEvents, nPar, nPar)
    if return_SNR == true
        SNRs = Array{Float64}(undef, nEvents)
        elapsed_time = @elapsed @showprogress desc="Computing Fishers and SNRs..." @threads for ii in 1:nEvents  

            Fishers[ii,:,:], SNRs[ii] = FisherMatrix(
                model,
                detector, 
                mc[ii], 
                eta[ii], 
                chi1[ii], 
                chi2[ii], 
                dL[ii], 
                theta[ii], 
                phi[ii], 
                iota[ii], 
                psi[ii], 
                tcoal[ii], 
                phiCoal[ii], 
                optional_param_reshaped[ii]..., 
                fmin= (fmin isa AbstractArray ? fmin[ii] : fmin), 
                fmax= (fmax isa AbstractArray ? fmax[ii] : fmax), 
                res = res, 
                useEarthMotion = useEarthMotion, 
                SNR_thres=SNR_thres, 
                alpha = alpha, 
                coordinate_shift = coordinate_shift, 
                return_SNR=true,
                optimization = optimization,
                call_number = call_number
            )
        end 

        println("Fisher matrices and SNRs computed!")
        if elapsed_time > 60.0
            elapsed_time = elapsed_time/60
            println("The evaluation took: ", elapsed_time, " minutes.")
        else
            println("The evaluation took: ", elapsed_time, " seconds.")
        end

        if auto_save == true  
            path = pwd()
            mkpath("output/"*name_folder)
            path = pwd()*"/output/"*name_folder*"/"
            date = Dates.now()
            date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
            h5open(path*"Fishers_SNRs.h5", "w") do file
                write(file, "Fishers", Fishers)
                attributes(file)["number_events"] = nEvents
                if typeof(detector) == Vector{Detector}
                    label = [detector[i].label for i in eachindex(detector)]
                    attributes(file)["Detectors"] = label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers and the SNRs for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*join(label, ",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format

                else
                    attributes(file)["Detectors"] = detector.label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers and the SNRs for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*detector.label*" detector and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format

                end     
                write(file, "SNRs", SNRs) 
                if save_catalog 
                    write(file, "mc", mc)
                    write(file, "eta", eta)
                    write(file, "chi1", chi1)
                    write(file, "chi2", chi2)
                    write(file, "dL", dL)
                    write(file, "theta", theta)
                    write(file, "phi", phi)
                    write(file, "iota", iota)
                    write(file, "psi", psi)
                    write(file, "tcoal", tcoal)
                    # Lambda is well defined, see above
                    write(file, "Lambda1", Lambda1)
                    write(file, "Lambda2", Lambda2)
                    # It would be nice to save fmin and fmax as well
                end
            end
        end
        return Fishers, SNRs
    else
        elapsed_time = @elapsed  @showprogress desc="Computing Fishers..."  @threads for ii in 1:nEvents  

                    Fishers[ii,:,:]=FisherMatrix(
                        model,
                        detector,
                        mc[ii],
                        eta[ii], 
                        chi1[ii],
                        chi2[ii],
                        dL[ii],
                        theta[ii],
                        phi[ii], 
                        iota[ii], 
                        psi[ii], 
                        tcoal[ii], 
                        phiCoal[ii],
                        optional_param_reshaped[ii]...,
                        fmin= (fmin isa AbstractArray ? fmin[ii] : fmin), 
                        fmax= (fmax isa AbstractArray ? fmax[ii] : fmax), 
                        res = res, 
                        useEarthMotion = useEarthMotion, 
                        SNR_thres=SNR_thres, 
                        alpha = alpha, 
                        coordinate_shift = coordinate_shift,
                        return_SNR=false,
                        optimization = optimization,
                        call_number = call_number
                    )
                end 
        println("Fisher matrices computed!")
        if elapsed_time > 60.0
            elapsed_time = elapsed_time/60
            println("The evaluation took: ", elapsed_time, " minutes.")
        else
            println("The evaluation took: ", elapsed_time, " seconds.")
        end
        if auto_save == true  
            path = pwd()
            mkpath("output/"*name_folder)
            path = pwd()*"/output/"*name_folder*"/"
            date = Dates.now()
            date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
            h5open(path*"Fishers.h5", "w") do file
                write(file, "Fishers", Fishers)
                attributes(file)["number_events"] = nEvents
                if typeof(detector) == Vector{Detector}
                    label = [detector[i].label for i in eachindex(detector)]
                    attributes(file)["Detectors"] = label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*join(label,",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format
                else
                    attributes(file)["Detectors"] = detector.label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*detector.label*" detector and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format
                end    
                if save_catalog
                    write(file, "mc", mc)
                    write(file, "eta", eta)
                    write(file, "chi1", chi1)
                    write(file, "chi2", chi2)
                    write(file, "dL", dL)
                    write(file, "theta", theta)
                    write(file, "phi", phi)
                    write(file, "iota", iota)
                    write(file, "psi", psi)
                    write(file, "tcoal", tcoal)
                    # Lambda is well defined, see above.
                    write(file, "Lambda1", Lambda1)
                    write(file, "Lambda2", Lambda2)
                    # It would be nice to save fmin and fmax as well
                end
            end
        end
        return Fishers
    end
        
end


"""
This function reads the Fisher matrices and/or the SNRs
    _read_Fishers_SNRs(path, SNR=true)

    #### Input arguments:
    -  `path` : string, path to the file
    -  `SNR` : bool, default true, if true the function reads also the SNRs

    #### Output:
    - `Fishers`  : matrix, Fisher Matrix
    - `SNRs`  : vector, SNRs (if SNR=true)

    #### Example:
    ```julia
    Fishers, SNRs = _read_Fishers_SNRs("output/BBH/Fishers_SNRs.h5")
    ```

"""
function _read_Fishers_SNRs(path; SNR::Bool=true)
    if SNR == true
        Fishers, SNRs = h5open(path, "r") do file
            println("Attributes: ", keys(attributes(file)))
            attributes_keys = keys(attributes(file))
            for i in 1:length(keys(attributes(file)))
                println(attributes_keys[i], ": ", read(attributes(file)[attributes_keys[i]]))
            end
            println("Keys: ", keys(file))
            Fishers = read(file, "Fishers")
            SNRs = read(file, "SNRs")
            return Fishers, SNRs
        end
    else
        Fishers = h5open(path, "r") do file
            println("Attributes: ", keys(attributes(file)))
            attributes_keys = keys(attributes(file))
            for i in 1:length(keys(attributes(file)))
                println(attributes_keys[i], ": ", read(attributes(file)[attributes_keys[i]]))
            end
            println("Keys: ", keys(file))
            Fishers = read(file, "Fishers")
            return Fishers
        end
    end
end

