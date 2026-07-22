module GW


# using DelimitedFiles
# using Plots
# using Interpolations
# using HDF5
# using BenchmarkTools
# using Test
# using Trapz 
# using ForwardDiff
# using Distributions
# using LinearAlgebra
# using QuadGK
# using Integrals
# using Random
# using ProgressMeter
# using Base.Threads
# using Polynomials
# using LaTeXStrings

include("utils.jl")
#include("waveform.jl")
include("waveforms/waveforms_master.jl")
include("detector/detector.jl")
include("catalog.jl")

using .waveform
using .UtilsAndConstants
using .detector
using .catalog

# export from waveform.jl
export Ampl, Phi, Pol, _npar, _event_type, _available_waveforms, _fcut, _finalspin, _radiatednrg, _tau_star, _list_polarizations, hphc
export TaylorF2, PhenomD, PhenomD_NRTidal, PhenomHM, PhenomNSBH, PhenomXAS, PhenomXHM, PhenomXPHM, PhenomD_TIGER, PhenomHM_TIGER, PhenomD_TIGER_spinless, PhenomHM_TIGER_spinless
export Model, GrModel, BgrModel


# export from detector.jl
export DetectorStructure, DetectorCoordinates, Detector, _readASD, _readPSD, getCoords, CE1Id_coordinates, CE1Id, CE2NM_coordinates,
         CE2NM, CE2NSW_coordinates, CE2NSW, ETS_coodinates, ETS, ETLS_coodinates, ETLS, ETMR_coordinates, ETMR, ETLMR_coordinates, ETLMR, 
         LIGO_L_coordinates, LIGO_L, LIGO_H_coordinates, LIGO_H, VIRGO_coordinates, VIRGO, KAGRA_coordinates, KAGRA, _available_detectors,
          _define_events, _deltLoc, _patternFunction, PolarizationDet, PhaseDet, Strain, SNR, FisherMatrix, _read_Fishers_SNRs

# export from catalog.jl
export GenerateCatalog, ReadCatalog, get_dL

# export from utils.jl
export GMsun_over_c3, GMsun_over_c2, uGpc, GMsun_over_c2_Gpc, REarth_km, clight_kms, clightGpc, Lamt_delLam_from_Lam12,
     _ra_dec_from_theta_phi_rad, _theta_phi_from_ra_dec_rad, CovMatrix, Errors, SkyArea, _orientationBigCircle, CovMatrix_Lamt_delLam




end
