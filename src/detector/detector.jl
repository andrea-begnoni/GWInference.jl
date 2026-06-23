module detector


using ..waveform
import ..UtilsAndConstants as uc     

##### JULIA PACKAGES
using Trapz
using DelimitedFiles
using Interpolations
using ForwardDiff
using LinearAlgebra
using Base.Threads
using BenchmarkTools
using HDF5
using ProgressMeter
using Base.Threads
using Dates


export DetectorStructure, DetectorCoordinates, Detector, _readASD, _readPSD, getCoords, CE1Id_coordinates, CE1Id, CE2NM_coordinates,
         CE2NM, CE2NSW_coordinates, CE2NSW, ETS_coodinates, ETS, ETLS_coodinates, ETLS, ETMR_coordinates, ETMR, ETLMR_coordinates, ETLMR, 
         LIGO_L_coordinates, LIGO_L, LIGO_H_coordinates, LIGO_H, VIRGO_coordinates, VIRGO, KAGRA_coordinates, KAGRA, _available_detectors,
          _define_events, _deltLoc, _patternFunction, PolarizationDet, PhaseDet, Strain, SNR, FisherMatrix, _read_Fishers_SNRs

# Here we define all the structures used inside this module
# define the detector structures
abstract type DetectorStructure end

# this one contains only the coordinates of the detector, the orientation and the arm aperture, all in radians
mutable struct DetectorCoordinates <: DetectorStructure
    latitude_rad::Float64
    longitude_rad::Float64
    orientation_rad::Float64
    arm_aperture_rad::Float64
end

# this one contains the detector structure plus the noise power spectrum and the frequency of the noise, the shape of the detector and a label to identify it when saving the results
mutable struct Detector <: DetectorStructure
    latitude_rad::Float64
    longitude_rad::Float64
    orientation_rad::Float64
    arm_aperture_rad::Float64
    shape::Char
    fNoise::AbstractArray
    psd::AbstractArray
    label::String
end

"""
function to read the amplitude spectral density of the noise
    _readASD("path/to/file.txt")

#### Input arguments:
-  `pathASD` : string, path to the file .txt containing the amplitude spectral density of the noise

#### Optional arguments:
-  `cols` : array, default [1,2], columns of the file containing the frequency and the amplitude spectral density

#### Output:
-  `fNoise` : array, frequency of the noise
-  `psd` : array, amplitude spectral density of the noise

#### Example:
```julia
fNoise, psd = _readASD("path/to/file.txt", [1,2])
```
"""
function _readASD(pathASD::String; cols = [1,2])
    fNoise, psd = open(pathASD) do file

        data = readdlm(file)
        fNoise = data[:, cols[1]]
        psd = data[:, cols[2]] .^ 2

        return fNoise, psd
    end
    return fNoise, psd
end

"""
function to read the power spectral density of the noise, which is the square of the amplitude spectral density
    _readPSD("path/to/file.txt")

#### Input arguments:
-  `pathPSD` : string, path to the file .txt containing the power spectral density of the noise

#### Optional arguments:
-  `cols` : array, default [1,2], columns of the file containing the frequency and the power spectral density

#### Output:
-  `fNoise` : array, frequency of the noise
-  `psd` : array, power spectral density of the noise

#### Example:
```julia
fNoise, psd = _readPSD("path/to/file.txt", [1,2])
```
"""
function _readPSD(pathPSD::String; cols = [1,2])
    fNoise, psd = open(pathPSD) do file

        data = readdlm(file)
        fNoise = data[:, cols[1]]
        psd = data[:, cols[2]]

        return fNoise, psd
    end
    return fNoise, psd
end


# function to extract coordinates of a detector structure
function getCoords(detector)
    return detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad
end

# Here we define the detectors present in the code, to create a new one, just follow the structure of the existing ones, for more information look at the tutorials

# Get the path to the directory of this file
PACKAGE_DIR = @__DIR__

# Go one step back in the path (from ""GW.jl/src" to "GW.jl")
PARENT_DIR = dirname(dirname(PACKAGE_DIR))

# Construct the path to the "useful_files" folder from the parent directory
PSDS_DIR = joinpath(PARENT_DIR, "useful_files/")

CE1Id_coordinates=detector.DetectorCoordinates(43.827 * pi/180, -112.825 * pi/180, -45 * pi/180, 90. * pi/180)
CE1Id = detector.Detector(43.827 * pi/180, -112.825 * pi/180, -45 * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer.txt")..., "CE1Id")

CE2NM_coordinates=detector.DetectorCoordinates(33.160 * pi/180, -106.480 * pi/180, -105. * pi/180, 90. * pi/180)
CE2NM = detector.Detector(33.160 * pi/180, -106.480 * pi/180, -105. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer_20km.txt")..., "CE2NM")

CE2NSW_coordinates=detector.DetectorCoordinates(-34. * pi/180, 145. * pi/180, 0. * pi/180, 90. * pi/180)
CE2NSW=detector.Detector(-34. * pi/180, 145. * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer_20km.txt")...,  "CE2NSW")

ETS_coodinates = detector.DetectorCoordinates((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 60. * pi/180)
ETS=detector.Detector((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 60. * pi/180, 'T', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETS")

ETLS_coodinates = detector.DetectorCoordinates((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 90. * pi/180)
ETLS = detector.Detector((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETLS")

ETMR_coordinates = detector.DetectorCoordinates((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 60. * pi/180)
ETMR = detector.Detector((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 60. * pi/180, 'T', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETMR")

ETLMR_coordinates = detector.DetectorCoordinates((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 90. * pi/180)
ETLMR = detector.Detector((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETLMR")

LIGO_L_coordinates = detector.DetectorCoordinates(30.563 * pi/180, -90.774 * pi/180, 242.71636956358617* pi/180, 90. * pi/180)
LIGO_L = detector.Detector(30.563 * pi/180, -90.774 * pi/180, 242.71636956358617* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/aligo_O3actual_L1.txt")...,  "LIGO_L")

LIGO_H_coordinates = detector.DetectorCoordinates(46.455 * pi/180, -119.408 * pi/180, 170.99924234706103* pi/180, 90. * pi/180)
LIGO_H = detector.Detector(46.455 * pi/180, -119.408 * pi/180, 170.99924234706103* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/aligo_O3actual_H1.txt")...,  "LIGO_H")

VIRGO_coordinates = detector.DetectorCoordinates(43.631 * pi/180, 10.504 * pi/180, 115.56756342034298* pi/180, 90. * pi/180)
VIRGO = detector.Detector(43.631 * pi/180, 10.504 * pi/180, 115.56756342034298* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/avirgo_O3actual.txt")..., "VIRGO")

KAGRA_coordinates = detector.DetectorCoordinates(36.412 * pi/180, 137.306 * pi/180, 15.396* pi/180, 90. * pi/180)
KAGRA = detector.Detector(36.412 * pi/180, 137.306 * pi/180, 15.396* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/kagra_128Mpc.txt")...,  "KAGRA")


# with this function you can check all the detectors available, it also prints the length of the arms for future interferomenters
function _available_detectors()
    println("The detectors available are: ")
    return ["CE1Id, 40km", "CE2NM, 20km", "CE2NSW, 20km", "ETS, 10km", "ETLS, 10km", "ETMR, 10km", "ETLMR, 10km", "LIGO_L", "LIGO_H", "VIRGO", "KAGRA"]
end

# this function returns the detector structure given the name of the detector
function _available_detectors(det::String)
    # check if det is among the tabulated detectors
    if det=="CE1Id"
        return CE1Id
    elseif det=="CE2NM"
        return CE2NM
    elseif det=="CE2NSW"
        return CE2NSW
    elseif det=="ETS"
        return ETS
    elseif det=="ETLS"
        return ETLS
    elseif det=="ETMR"
        return ETMR
    elseif det=="ETLMR"
        return ETLMR
    elseif det=="LIGO_L"
        return LIGO_L
    elseif det=="LIGO_H"
        return LIGO_H
    elseif det=="VIRGO"
        return VIRGO
    elseif det=="KAGRA"
        return KAGRA
    else
        println("Detector not available. The ones available are: CE1Id, CE2NM, CE2NSW, ETS (triangular), ETLS, ETMR (triangular), ETLMR, LIGO_L, LIGO_H, VIRGO, KAGRA")
    end
end

# function to check that the event have sense and the variable are ordered correnctly
function _define_events(mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal)

    if mc <= 0 || mc >= 200
        error("Chirp mass out of range, mc = $mc")
    end
    if eta <= 0 || eta > 0.25   
        error("Symmetric mass ratio out of range, eta = $eta")
    end
    if dL <= 0
        error("Luminosity distance out of range, dL = $dL")
    end
    if theta < 0 || theta > pi
        error("Sky position angle out of range, theta = $theta")
    end
    if phi < 0 || phi > 2 * pi
        error("Sky position angle out of range, phi = $phi")
    end
    if iota < 0 || iota > pi
        error("Inclination angle out of range, iota = $iota")
    end
    if psi < 0 || psi > 2 * pi
        error("Polarisation angle out of range, psi = $psi")
    end
    if tcoal < 0
        error("Time of coalescence out of range, tcoal = $tcoal")
    end
    if phiCoal < 0 || phiCoal > 2 * pi
        error("GW phase at coalescence out of range, phiCoal = $phiCoal")
    end
    if chi1 < -1 || chi1 > 1
        error("Spin component out of range, chi1 = $chi1")
    end
    if chi2 < -1 || chi2 > 1
        error("Spin component out of range, chi2 = $chi2")
    end

    println("Everything looks fine")
    return 0
end

"""
The function computes the time (in seconds) it takes for the GW to go from Earth center, which is the point to which all GWs are referred, to the location of the detector.
    _deltLoc(theta, phi, tRef, DetectorCoordinates) 

#### Input arguments:
-  `theta` : float, in rad
-  `phi`   : float, in rad
-  `tRef`  : float, in LMST, it is the time at which the GW arrives at the center of the Earth, it is equal to `tcoal` if we are not considering the Earth motion during the measurement
-  `DetectorCoordinates` : structure, containing the coordinates of the detector

#### Optional arguments:
-  `REarth_km` : float, default 6371.0, Earth radius in km
-  `clight_kms` : float, default 299792.458, speed of light in km/s

#### Output:
- `delt`  : float, in seconds

#### Example:
```julia
delt = _deltLoc(0.1, 0.2, 0.3, CE1Id_coordinates)
```
"""
function _deltLoc(
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    tRef,
    DetectorCoordinates::DetectorStructure;
    REarth_km = uc.REarth_km,
    clight_kms = uc.clight_kms,
)

    # Time needed to go from Earth center to detector location

    ra, dec = uc._ra_dec_from_theta_phi_rad(theta, phi)

    comp1 = 
        @. cos(dec) *
        cos(ra) *
        cos(DetectorCoordinates.latitude_rad) *
        cos(DetectorCoordinates.longitude_rad + 2.0 * pi * tRef)
    comp2 =
        @. cos(dec) *
        sin(ra) *
        cos(DetectorCoordinates.latitude_rad) *
        sin(DetectorCoordinates.longitude_rad + 2.0 * pi * tRef)
    comp3 = sin(dec) * sin(DetectorCoordinates.latitude_rad)
    # The minus sign arises from the definition of the unit vector pointing to the source
    delt = @.  -REarth_km * (comp1 + comp2 + comp3) / clight_kms

    return delt # in seconds
end

# note order is important
include("strain.jl")
include("snr.jl")
include("fisher.jl")

end# of module
