"""
Compute the *signal-to-noise-ratio*, SNR, as a function of the parameters of the event, as measured a single detector. 
The SNR is computed as the square root of the integral of the signal-to-noise ratio squared.
The integral is computed using the trapezoidal rule.
This function is the main function to compute the SNR, then using multiple-dispatch this function is called when considering a network of detectors and more than a single event.

SNR(model, detector, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, Lambda1=0.0, Lambda2=0.0; fmin=2.0, fmax=nothing, res=1000, useEarthMotion=false, ampl_precomputation=nothing)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
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
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `ampl_precomputation` : array, default nothing, precomputed amplitude

    #### Output:
    - `SNR`  : float, signal-to-noise ratio

    #### Example:
    ```julia
    SNR = SNR(PhenomD(), CE1Id , 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4)
    ```
"""
function SNR(model::Model,
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
    optional_param...;
    fmin::Float64=2.0,
    fmax::Union{Nothing, Float64}=nothing,
    res = 1000,
    useEarthMotion::Bool = false,
    ampl_precomputation = nothing,
)

    #extract tidal diformabilites
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
    
    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))
    # Out of the PSD range, we use a constant value of 1, which results in completely negligible contributions
    psdGrid = linear_interpolation(detector.fNoise, detector.psd, extrapolation_bc = 1.0)(fgrid)
    detectorCoordinates = DetectorCoordinates(
        detector.latitude_rad,
        detector.longitude_rad,
        detector.orientation_rad,
        detector.arm_aperture_rad
    )

        if ampl_precomputation === nothing
           
            ampl_precomputation = PolAbs(
                model,
                fgrid,
                mc,
                eta,
                chi1,
                chi2,
                dL,
                iota,
                optional_param... 
            )
        else
            ampl_precomputation = ampl_precomputation
        end

    if detector.shape == 'L' # in Julia '' indicates a char, while "" a string 
        Aps, Acs = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            theta,
            phi,
            psi,
            tcoal,
            useEarthMotion = useEarthMotion,
        )
        Atot = Aps .* Aps .+ Acs .* Acs
        SNR = 2.0 * sqrt(trapz(fgrid, Atot ./ psdGrid))

    elseif detector.shape == 'T'
        SNRsq_Tshape = zeros(3)
        
        # The signal in 3 arms sums to zero for geometrical reasons, so we can use this to skip some calculations

        Aps1, Acs1 = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            theta,
            phi,
            psi,
            tcoal,
            alpha = 0.0,
            useEarthMotion = useEarthMotion,
        )
        Atot1 = Aps1 .* Aps1 .+ Acs1 .* Acs1
        Aps2, Acs2 = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            theta,
            phi,
            psi,
            tcoal,
            alpha = 60.0,
            useEarthMotion = useEarthMotion,
        )
        Atot2 = Aps2 .* Aps2 .+ Acs2 .* Acs2
        Aps3, Acs3 = -(Aps1 .+ Aps2), -(Acs1 .+ Acs2)
        Atot3 = Aps3 .* Aps3 .+ Acs3 .* Acs3
        SNRsq_Tshape[1] = trapz(fgrid, Atot1 ./ psdGrid)
        SNRsq_Tshape[2] = trapz(fgrid, Atot2 ./ psdGrid)
        SNRsq_Tshape[3] = trapz(fgrid, Atot3 ./ psdGrid)


        SNR = 2.0 * sqrt(sum(SNRsq_Tshape)) # The factor of two arises by cutting the integral from 0 to infinity
    
    end

    return SNR

end

"""
This function computes the *signal-to-noise-ratio*, SNR, as a function of the parameters of the event, as measured by a NETWORK of detectors.
It relies on the function SNR(..., detector::Detector, ...), which computes the SNR for a single detector.
where the dots indicate the parameters equal to the previous function call.

#### Optional arguments:
-  `precomputation` : bool, default true, if true the waveform is called just once and then reused for each detector.
-  `ampl_precomputation` : array, default nothing, precomputed amplitude of the polarizations (as returned by PolAbs), used if `precomputation` is true.

#### Example:
```julia
    SNR(... [CE1Id, CE2NM], ...)
```

"""
function SNR(model::Model,
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
    optional_param...;
    fmin::Float64=2.0,
    fmax::Union{Nothing, Float64} = nothing,
    res = 1000,
    useEarthMotion::Bool = false,
    precomputation::Bool = true,
    ampl_precomputation = nothing,
)
##########################
## This part is to precompute the amplitude of the waveform which is the longest part of the computation
    
    #extract tidal diformabilites
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

    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))

    if precomputation == true && isnothing(ampl_precomputation)
        ampl_precomputation = waveform.PolAbs(
            model,
            fgrid,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            iota,
            optional_param...;
        )
    elseif precomputation == false
        ampl_precomputation = nothing
    end

########################
    SNRList = Vector{Float64}(undef,length(detector))
    for i in 1:length(detector)
        SNRList[i] = SNR(
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
            optional_param...,
            fmin=fmin,
            fmax=fmax,
            res=res,
            useEarthMotion = useEarthMotion,
            ampl_precomputation = ampl_precomputation
        )
    end

    return sqrt(sum(SNRList.^2))
end

"""
This function computes the *signal-to-noise-ratio*, SNR, as a function of the parameters of an ARRAY of event, as measured by a SINGLE detector or a NETWORK of detectors.
It relies on the function SNR(..., detector::Detector, ...), which computes the SNR for a single detector.
where the dots indicate the parameters equal to the previous function call.

It is possible to save the SNRs in a file, if the optional argument `auto_save` is set to true. Use `name_folder` : string, to decide the name of the folder where the results are saved,
if the default is left it saves BBH in the folder "output/BBH" and so on for each source type. The file is saved in the folder `output/name_folder/SNRs.h5` 
and contains the SNRs for the events in the catalog. It contains also the parameters of the events if the optional argument `save_catalog` is set to true.
"""
function SNR(model::Model,
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
    optional_param...;
    fmin::Union{Float64, AbstractArray}=2.0,
    fmax::Union{Nothing, Float64, AbstractArray} = nothing,
    res = 1000,
    auto_save::Bool = false,
    name_folder = "BBH",
    save_catalog::Bool = false,
    useEarthMotion::Bool = false,
    precomputation::Bool = true,
)
    nEvents = length(mc)
    SNRs = Vector{Float64}(undef, nEvents)

    if(fmin isa AbstractArray && length(fmin) != nEvents)
        throw(ArgumentError("fmin must be an array of the same length as the number of events (or otherwise a single scalar Float64)"))
    end

    if(fmax isa AbstractArray && length(fmax) != nEvents)
        throw(ArgumentError("fmax must be an array of the same length as the number of events (or otherwise a single scalar Float64 or Nothing)"))
    end

    #check correct length of optional_parameters
    for op in optional_param
        if length(op) != nEvents
            trow(DomainError(op, "Lenght of each optional parameter must be a vectors with the same length as the number of events!"))
        end
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

    elapsed_time = @elapsed @showprogress desc="Computing SNRs..."  @threads for ii in 1:nEvents  
                    
        optional_param_ii = []
        for op in optional_param
            append!(optional_param_ii, op[ii])
        end

        SNRs[ii]=SNR(
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
            optional_param_ii..., 
            fmin= (fmin isa AbstractArray ? fmin[ii] : fmin), 
            fmax= (fmax isa AbstractArray ? fmax[ii] : fmax), 
            res = res, 
            useEarthMotion = useEarthMotion,
            precomputation = precomputation)
    end 

    println("SNRs computed!")
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
        h5open(path*"SNRs.h5", "w") do file
            attributes(file)["number_events"] = nEvents
            if typeof(detector) == Vector{Detector}
                label = [detector[i].label for i in eachindex(detector)]
                attributes(file)["Detectors"] = label
                attributes(file)["What_this_file_contains"] = "This file contains the SNR for the "*string(nEvents)*" events of the catalog 
                                generated with the "*string(typeof(model))*" waveform model. The SNR is calculated for the "*join(label,",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                attributes(file)["date"] = date_format
            else
                attributes(file)["Detectors"] = detector.label
                attributes(file)["What_this_file_contains"] = "This file contains the SNR for the "*string(nEvents)*" events of the catalog 
                                generated with the "*string(typeof(model))*" waveform model. The SNR is calculated for the "*detector.label*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                attributes(file)["date"] = date_format
            end   
            
            write(file, "SNRs", SNRs) 
            if save_catalog && typeof(model) <: GrModel 
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
                # Lambdas are well defined!
                # See above
                write(file, "Lambda1", Lambda1)
                write(file, "Lambda2", Lambda2)
                # It would be nice to save fmin and fmax as well
            end
        end
    end

    return SNRs
end