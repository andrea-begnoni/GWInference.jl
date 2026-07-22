#! format: off
module waveform
### This is the first module in the pipeline (after the catalog generation), it computes the waveform

# The waveforms presented here are adapted and modified from LALSimulation and GWFAST (https://github.com/CosmoStatGW/gwfast)


import ..UtilsAndConstants as uc    

### Import Julia packages relevant for this module
using DelimitedFiles
using Interpolations
using ForwardDiff
using Roots
using Elliptic
using LinearAlgebra


export TaylorF2, PhenomD, PhenomD_NRTidal, PhenomHM, PhenomNSBH, PhenomXAS, PhenomXHM, PhenomXPHM, PhenomD_TIGER, PhenomHM_TIGER, PhenomD_TIGER_spinless, PhenomHM_TIGER_spinless
export Model, GrModel, BgrModel
export Ampl, Phi, PolAbs, Pol, _npar, _event_type, _available_waveforms, _fcut, _finalspin, _radiatednrg, _tau_star, _list_polarizations, hphc

# Define an abstract type for the models
abstract type Model end
abstract type GrModel <: Model end    # Model for GrWaveforms
abstract type BgrModel <: Model end   # Model for Wavefrom including beyond GR polarizations


# Define a function to give an error message if the model is not implemented
function Ampl(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    # Implementation specific to each model
    error("Ampl not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function Phi(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_PHI = 0.018,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    # Implementation specific to each model
    error("Phi not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function Pol(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    itoa,
)
    # Implementation specific to each model
    error("Pol not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function PolAbs(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
)
    # Implementation specific to each model
    error("PolAbs not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

function _npar(model::Model)
    # Implementation specific to each model
    error("_npar not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function _list_polarizations(model::Model) 
    # Implementation specific to each model
    error("_list_polarizations not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
 end

# All standard GR polarizations have the same polarizations
"""
Returns a list of the names of the polarizations of this waveform.

    _list_polarizations(model)

    The order of polarizations handed in this list is the same as the signals returned by Pol(model, ...). 
"""
function _list_polarizations(model::GrModel) 
    return ["plus", "cross"]
 end

# Define concrete types for each model
# GR Waveforms
struct PhenomD <: GrModel 
    event_type::String 
    PhenomD() = new("BBH")
end

struct PhenomHM <: GrModel
    event_type::String 
    PhenomHM() = new("BBH")
end

struct PhenomD_NRTidal <: GrModel
    event_type::String 
    PhenomD_NRTidal() = new("BNS")
end

struct PhenomNSBH <: GrModel 
    event_type::String 
    PhenomNSBH() = new("NSBH")
end

struct PhenomXAS <: GrModel 
    event_type::String 
    PhenomXAS() = new("BBH")
end

struct PhenomXHM <: GrModel 
    event_type::String 
    PhenomXHM() = new("BBH")
end

struct PhenomXPHM <: GrModel
    event_type::String
    PhenomXPHM() = new("BBH")
end

struct TaylorF2 <: GrModel 
    event_type::String 
    TaylorF2(event_type::String = "BBH") = new(event_type)
end

#Beyond GR Waveforms
struct PhenomD_TIGER <: BgrModel 
    PNorder::Float64
    event_type::String 
    PhenomD_TIGER(PNorder::Float64) = new(PNorder, "BBH")
end

struct PhenomHM_TIGER <: BgrModel 
    PNorder::Float64
    event_type::String 
    PhenomHM_TIGER(PNorder::Float64) = new(PNorder, "BBH")
end

struct PhenomD_TIGER_spinless <: BgrModel 
    PNorder::Float64
    event_type::String 
    PhenomD_TIGER_spinless(PNorder::Float64) = new(PNorder, "BBH")
end

struct PhenomHM_TIGER_spinless <: BgrModel 
    PNorder::Float64
    event_type::String 
    PhenomHM_TIGER_spinless(PNorder::Float64) = new(PNorder, "BBH")
end

"""
Returns the event_type of a struct<:Model as a string.
"""
function _event_type(model::Model)
    return model.event_type
end

function _available_waveforms()
    return ["TaylorF2", "PhenomD", "PhenomHM", "PhenomD_NRTidal", "PhenomNSBH", "PhenomXAS", "PhenomXHM", "PhenomXPHM", "PhenomD_TIGER", "PhenomHM_TIGER", "PhenomD_TIGER_spinless", "PhenomHM_TIGER_spinless"]
end

#@doc "Function to check the available waveforms and return the corresponding model."
"""
Function to check the available waveforms and return the corresponding model.
"""
function _available_waveforms(waveform::String)

    if waveform == "TaylorF2"
        return TaylorF2()
    elseif waveform == "PhenomD"
        return PhenomD()
    elseif waveform == "PhenomHM"
        return PhenomHM()
    elseif waveform == "PhenomD_NRTidal"
        return PhenomD_NRTidal()
    elseif waveform == "PhenomNSBH"
        return PhenomNSBH()
    elseif waveform == "PhenomXAS"
        return PhenomXAS()
    elseif waveform == "PhenomXHM"
        return PhenomXHM()
    elseif waveform == "PhenomXPHM"
        return PhenomXPHM()
    else
        error("Waveform not available. Choose between: $(_available_waveforms())")
    end
end
##############################################################################
#   INCLUDE ALL THE WAVEFORMS
##############################################################################

# GR waveforms
include("TaylorF2.jl")
include("PhenomD.jl")
include("PhenomHM.jl")
include("PhenomD_NRTidalv2.jl")
include("PhenomNSBH.jl")
include("PhenomXAS.jl")
include("PhenomXHM.jl")
include("PhenomXPHM.jl")
include("ConnectionFunctionsXAS.jl") # This is needed for PhenomXHM

#Beyond GR waveforms
include("PhenomD_TIGER.jl")
include("PhenomHM_TIGER.jl")
include("PhenomD_TIGER_spinless.jl")
include("PhenomHM_TIGER_spinless.jl")

##############################################################################
#   STRUCTURE USED IN THE MODULE
##############################################################################

struct TF2EccCoeffsStruct
    zero::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    twoV::Union{Float64,ForwardDiff.Dual}
    twoV0::Union{Float64,ForwardDiff.Dual}
    threeV::Union{Float64,ForwardDiff.Dual}
    threeV0::Union{Float64,ForwardDiff.Dual}
    fourV4::Union{Float64,ForwardDiff.Dual}
    fourV2V02::Union{Float64,ForwardDiff.Dual}
    fourV04::Union{Float64,ForwardDiff.Dual}
    fiveV5::Union{Float64,ForwardDiff.Dual}
    fiveV3V02::Union{Float64,ForwardDiff.Dual}
    fiveV2V03::Union{Float64,ForwardDiff.Dual}
    fiveV05::Union{Float64,ForwardDiff.Dual}
    sixV6::Union{Float64,ForwardDiff.Dual}
    sixV4V02::Union{Float64,ForwardDiff.Dual}
    sixV3V03::Union{Float64,ForwardDiff.Dual}
    sixV2V04::Union{Float64,ForwardDiff.Dual}
    sixV06::Union{Float64,ForwardDiff.Dual}
end

struct TF2coeffsStructure
    zero::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
    three::Union{Float64,ForwardDiff.Dual}
    four::Union{Float64,ForwardDiff.Dual}
    five::Union{Float64,ForwardDiff.Dual}
    five_log::Union{Float64,ForwardDiff.Dual}
    six::Union{Float64,ForwardDiff.Dual}
    six_log::Union{Float64,ForwardDiff.Dual}
    seven::Union{Float64,ForwardDiff.Dual}
end

struct TF2coeffsStructure_BGR
    minus_two::Union{Float64,ForwardDiff.Dual} # added for PhenomD_TIGER
    zero::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
    three::Union{Float64,ForwardDiff.Dual}
    four::Union{Float64,ForwardDiff.Dual}
    five::Union{Float64,ForwardDiff.Dual}
    five_log::Union{Float64,ForwardDiff.Dual}
    six::Union{Float64,ForwardDiff.Dual}
    six_log::Union{Float64,ForwardDiff.Dual}
    seven::Union{Float64,ForwardDiff.Dual}
end

struct PhiInspcoeffsStructure
    initial_phasing::Union{Float64,ForwardDiff.Dual}
    two_thirds::Union{Float64,ForwardDiff.Dual}
    third::Union{Float64,ForwardDiff.Dual}
    third_log::Union{Float64,ForwardDiff.Dual}
    log::Union{Float64,ForwardDiff.Dual}
    min_third::Union{Float64,ForwardDiff.Dual}
    min_two_thirds::Union{Float64,ForwardDiff.Dual}
    min_one::Union{Float64,ForwardDiff.Dual}
    min_four_thirds::Union{Float64,ForwardDiff.Dual}
    min_five_thirds::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    four_thirds::Union{Float64,ForwardDiff.Dual}
    five_thirds::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
end

struct PhiInspcoeffsStructure_BGR
    initial_phasing::Union{Float64,ForwardDiff.Dual}
    two_thirds::Union{Float64,ForwardDiff.Dual}
    third::Union{Float64,ForwardDiff.Dual}
    third_log::Union{Float64,ForwardDiff.Dual}
    log::Union{Float64,ForwardDiff.Dual}
    min_third::Union{Float64,ForwardDiff.Dual}
    min_two_thirds::Union{Float64,ForwardDiff.Dual}
    min_one::Union{Float64,ForwardDiff.Dual}
    min_four_thirds::Union{Float64,ForwardDiff.Dual}
    min_five_thirds::Union{Float64,ForwardDiff.Dual}
    min_seven_thirds::Union{Float64,ForwardDiff.Dual} # added for PhenomD_TIGER
    one::Union{Float64,ForwardDiff.Dual}
    four_thirds::Union{Float64,ForwardDiff.Dual}
    five_thirds::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
end

struct AcoeffsStructure
    two_thirds::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    four_thirds::Union{Float64,ForwardDiff.Dual}
    five_thirds::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
    seven_thirds::Union{Float64,ForwardDiff.Dual}
    eight_thirds::Union{Float64,ForwardDiff.Dual}
    three::Union{Float64,ForwardDiff.Dual}
end

struct PhislmpStructure
    two_one::Union{Float64, Vector{Float64}}
    two_two::Union{Float64, Vector{Float64}}
    three_two::Union{Float64, Vector{Float64}}
    three_three::Union{Float64, Vector{Float64}}
    four_three::Union{Float64, Vector{Float64}}
    four_four::Union{Float64, Vector{Float64}}
end

struct AmplslmpStructure
    two_one::Union{Float64, Vector{Float64}}
    two_two::Union{Float64, Vector{Float64}}
    three_two::Union{Float64, Vector{Float64}}
    three_three::Union{Float64, Vector{Float64}}
    four_three::Union{Float64, Vector{Float64}}
    four_four::Union{Float64, Vector{Float64}}
end

struct Phase22Struct
    fPhaseMatchIN::Union{Float64,ForwardDiff.Dual}
    fPhaseMatchIM::Union{Float64,ForwardDiff.Dual}
    phi0::Union{Float64,ForwardDiff.Dual}
    phi1::Union{Float64,ForwardDiff.Dual}
    phi2::Union{Float64,ForwardDiff.Dual}
    phi3::Union{Float64,ForwardDiff.Dual}
    phi4::Union{Float64,ForwardDiff.Dual}
    phi5::Union{Float64,ForwardDiff.Dual}
    phi5L::Union{Float64,ForwardDiff.Dual}
    phi6::Union{Float64,ForwardDiff.Dual}
    phi6L::Union{Float64,ForwardDiff.Dual}
    phi7::Union{Float64,ForwardDiff.Dual}
    phi8::Union{Float64,ForwardDiff.Dual}
    phi8L::Union{Float64,ForwardDiff.Dual}
    phi9::Union{Float64,ForwardDiff.Dual}
    phi9L::Union{Float64,ForwardDiff.Dual}
    dphase0::Union{Float64,ForwardDiff.Dual}
    dphi0::Union{Float64,ForwardDiff.Dual}
    dphi1::Union{Float64,ForwardDiff.Dual}
    dphi2::Union{Float64,ForwardDiff.Dual}
    dphi3::Union{Float64,ForwardDiff.Dual}
    dphi4::Union{Float64,ForwardDiff.Dual}
    dphi5::Union{Float64,ForwardDiff.Dual}
    dphi6::Union{Float64,ForwardDiff.Dual}
    dphi6L::Union{Float64,ForwardDiff.Dual}
    dphi7::Union{Float64,ForwardDiff.Dual}
    dphi8::Union{Float64,ForwardDiff.Dual}
    dphi8L::Union{Float64,ForwardDiff.Dual}
    dphi9::Union{Float64,ForwardDiff.Dual}
    dphi9L::Union{Float64,ForwardDiff.Dual}
    a0coloc::Union{Float64,ForwardDiff.Dual}
    a1coloc::Union{Float64,ForwardDiff.Dual}
    a2coloc::Union{Float64,ForwardDiff.Dual}
    a3coloc::Union{Float64,ForwardDiff.Dual}
    a4coloc::Union{Float64,ForwardDiff.Dual}
    b0coloc::Union{Float64,ForwardDiff.Dual}
    b1coloc::Union{Float64,ForwardDiff.Dual}
    b2coloc::Union{Float64,ForwardDiff.Dual}
    b3coloc::Union{Float64,ForwardDiff.Dual}
    b4coloc::Union{Float64,ForwardDiff.Dual}
    cLcoloc::Union{Float64,ForwardDiff.Dual}
    C2Int::Union{Float64,ForwardDiff.Dual}
    c0coloc::Union{Float64,ForwardDiff.Dual}
    c1coloc::Union{Float64,ForwardDiff.Dual}
    c2coloc::Union{Float64,ForwardDiff.Dual}
    c4coloc::Union{Float64,ForwardDiff.Dual}
    C2MRD::Union{Float64,ForwardDiff.Dual}
    sigma1::Union{Float64,ForwardDiff.Dual}
    sigma2::Union{Float64,ForwardDiff.Dual}
    sigma3::Union{Float64,ForwardDiff.Dual}
    sigma4::Union{Float64,ForwardDiff.Dual}
    sigma5::Union{Float64,ForwardDiff.Dual}
    C1Int::Union{Float64,ForwardDiff.Dual}
    C1MRD::Union{Float64,ForwardDiff.Dual}
end


struct Ampl22Struct
    gamma1::Union{Float64,ForwardDiff.Dual}
    gamma2::Union{Float64,ForwardDiff.Dual}
    gamma3::Union{Float64,ForwardDiff.Dual}
end

##############################################################################
#   HELPER FUNCTIONS, NEEDED BY MORE THAN ONE WAVEFORM
##############################################################################

function _readQNMgrid_a(pathWF::String)
    return readdlm(pathWF * "QNMData_a.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end

function _readQNMgrid_fring(pathWF::String)
    return readdlm(pathWF * "QNMData_fring.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end

function _readQNMgrid_fdamp(pathWF::String)
    return readdlm(pathWF * "QNMData_fdamp.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end


##############################################################################
#   HELPER FUNCTIONS, NEEDED ONLY BY ONE WAVEFORM
##############################################################################

##############################################################################
#
#                              TaylorF2
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::TaylorF2, Lambda1, Lambda2)
    more_par = 0
    if Lambda1 > 0.
        more_par += 1
    end

    if Lambda2 > 0.
        more_par += 1
    end

    return 11 + more_par
end

# function Phi(model::TaylorF2, f, mc, eta, chi1, chi2, Lambda1, Lambda2; GMsun_over_c3 = uc.GMsun_over_c3)
#     return Phi(model, f, mc, eta, chi1, chi2, GMsun_over_c3 = GMsun_over_c3)
# end

function Ampl(model::TaylorF2, f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    return Ampl(model, f, mc, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end

##############################################################################
#
#                              PhenomD
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomD)
    return 11
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1,
    Lambda2,
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    container = nothing,
)
    return Phi(model, f, mc, eta, chi1, chi2, fInsJoin_PHI=fInsJoin_PHI, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3, container=container)
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing,
)

    return Ampl(model, f, mc, eta, chi1, chi2, dL, fcutPar = fcutPar, fInsJoin_Ampl = fInsJoin_Ampl, GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc, container = container)
end

##############################################################################
#
#                              PhenomD_NRTidalv2
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomD_NRTidal)
    return 13
end
    
##############################################################################
#
#                              PhenomNSBH
#
##############################################################################
"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomNSBH)
    return 12
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, 
    Lambda,
    Lambda2;
    fInsJoin = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)

    return Phi(model, f, mc, eta, chi1, chi2, Lambda; fInsJoin = fInsJoin, fcutPar = fcutPar, GMsun_over_c3 = GMsun_over_c3)
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, #assumed to be chi2=0
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    return Ampl(model, f, mc, eta, chi1, chi2, dL, Lambda1, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3, GMsun_over_c2_Gpc=GMsun_over_c2_Gpc)
end

##############################################################################
#
#                              PhenomHM
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomHM)
    return 11
end

function Phi(model::PhenomHM, f, mc, eta, chi1, chi2, Lambda1, Lambda2; GMsun_over_c3 = uc.GMsun_over_c3)
    return Phi(model, f, mc, eta, chi1, chi2, GMsun_over_c3 = GMsun_over_c3)
end

function Ampl(model::PhenomHM, f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    return Ampl(model, f, mc, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end

##############################################################################
#
#                              PhenomHM_TIGER
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomHM_TIGER)
    return 12
end

function Phi(model::PhenomHM_TIGER, f, mc, eta, chi1, chi2, optional_param...; GMsun_over_c3 = uc.GMsun_over_c3)
    o1 = optional_param[1]
    return Phi(model, f, mc, eta, chi1, chi2, o1, GMsun_over_c3 = GMsun_over_c3)
end

function Ampl(model::PhenomHM_TIGER, f, mc, eta, chi1, chi2, dL, optional_param...; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    o1 = optional_param[1]
    return Ampl(model, f, mc, chi1, chi2, o1, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end

function _list_polarizations(model::PhenomHM_TIGER) 
    return ["plus", "cross"]
 end


##############################################################################
#
#                              PhenomHM_TIGER_spinless
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomHM_TIGER_spinless)
    return 12
end

function Phi(model::PhenomHM_TIGER_spinless, f, mc, eta, chi1, chi2, optional_param...; GMsun_over_c3 = uc.GMsun_over_c3)
    o1 = optional_param[1]
    return Phi(model, f, mc, eta, chi1, chi2, o1, GMsun_over_c3 = GMsun_over_c3)
end

function Ampl(model::PhenomHM_TIGER_spinless, f, mc, eta, chi1, chi2, dL, optional_param...; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    o1 = optional_param[1]
    return Ampl(model, f, mc, chi1, chi2, o1, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end

function _list_polarizations(model::PhenomHM_TIGER_spinless) 
    return ["plus", "cross"]
 end

##############################################################################
#
#                              PhenomXAS
#
##############################################################################


"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomXAS)
    return 11
end

"""
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomXAS,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1,
    Lambda2,
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)

    return Phi(model, f, mc, eta, chi1, chi2, fInsJoin_PHI=fInsJoin_PHI, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3)
end


"""
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomXAS,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.3,
    IntAmpVersion=104,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    return Ampl(model, f, mc, eta, chi1, chi2, dL, fcutPar=fcutPar, IntAmpVersion=IntAmpVersion, GMsun_over_c3=GMsun_over_c3, GMsun_over_c2_Gpc=GMsun_over_c2_Gpc)
end


##############################################################################
#
#                              PhenomXHM
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomXHM)
    return 11
end

function _npar(model::PhenomXPHM)
    return 15
end

function Phi(model::PhenomXHM, f, mc, eta, chi1, chi2, Lambda1, Lambda2; GMsun_over_c3 = uc.GMsun_over_c3)
    return Phi(model, f, mc, eta, chi1, chi2, GMsun_over_c3 = GMsun_over_c3)
end

function Ampl(model::PhenomXHM, f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    return Ampl(model, f, mc, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end


##############################################################################
#
#                              PhenomD_TIGER_spinless
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomD_TIGER_spinless)
    return 12 #13  # ANDREA: I think this should be 12, not 13
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomD_TIGER_spinless,
    f,
    mc,
    eta,
    chi1,
    chi2,
    optional_param... ;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)

    o1 = optional_param[1]
    return Phi(model, f, mc, eta, chi1, chi2, o1, fInsJoin_PHI=fInsJoin_PHI, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3)
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomD_TIGER_spinless,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    optional_param... ;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    o1 = optional_param[1]
    return Ampl(model, f, mc, eta, chi1, chi2, dL, o1, fcutPar = fcutPar, fInsJoin_Ampl = fInsJoin_Ampl, GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc)
end

"""
Returns polarizations of the PhenomD_TIGER_spinless waveform
"""
function _list_polarizations(model::PhenomD_TIGER_spinless) 
    return ["plus", "cross"]
 end


##############################################################################
#
#                              PhenomD_TIGER
#
##############################################################################

"""
Returns the number of parameter of a struct<:Model as integer number. 
"""
function _npar(model::PhenomD_TIGER)
    return 12 #13  # ANDREA: I think this should be 12, not 13
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomD_TIGER,
    f,
    mc,
    eta,
    chi1,
    chi2,
    optional_param... ;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)

    o1 = optional_param[1]
    return Phi(model, f, mc, eta, chi1, chi2, o1, fInsJoin_PHI=fInsJoin_PHI, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3)
end

""" 
helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomD_TIGER,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    optional_param... ;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    o1 = optional_param[1]
    return Ampl(model, f, mc, eta, chi1, chi2, dL, o1, fcutPar = fcutPar, fInsJoin_Ampl = fInsJoin_Ampl, GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc)
end

"""
Returns polarizations of the PhenomD_TIGER waveform
"""
function _list_polarizations(model::PhenomD_TIGER) 
    return ["plus", "cross"]
 end

##############################################################################
#   FUNCTIONS USED BY MANY WAVEFORMS
##############################################################################

"""
Compute the spin of the final object, as in LALSimIMRPhenomD_internals.c line 161 and 142, which is taken from `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_ eq. (3.6).
Valid for all the models considered in this code.


#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (float) spin of the final object.

"""
function _finalspin(model::Model, eta, chi1, chi2)

    eta2 = eta * eta
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1 * m1 * chi1 + m2 * m2 * chi2)

    af1 =
        eta * (
            3.4641016151377544 - 4.399247300629289 * eta + 9.397292189321194 * eta2 -
            13.180949901606242 * eta2 * eta
        )
    af2 =
        eta * (
            s * (
                (1.0 / eta - 0.0850917821418767 - 5.837029316602263 * eta) +
                (0.1014665242971878 - 2.0967746996832157 * eta) * s
            )
        )
    af3 =
        eta * (
            s * (
                (-1.3546806617824356 + 4.108962025369336 * eta) * s * s +
                (-0.8676969352555539 + 2.064046835273906 * eta) * s * s * s
            )
        )
    return af1 + af2 + af3
end

"""
Compute the total radiated energy, as in `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_ eq. (3.7) and (3.8).
Valid for all the models considered in this code.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (float) total radiated energy.

"""
function _radiatednrg(model::Model, eta, chi1, chi2)


    eta2 = eta * eta
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)

    EradNS =
        eta * (
            0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 +
            3.352411249771192 * eta2 * eta
        )

    return (
        EradNS * (
            1.0 +
            (
                -0.0030302335878845507 - 2.0066110851351073 * eta +
                7.7050567802399215 * eta2
            ) * s
        )
    ) / (
        1.0 +
        (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s
    )
end

"""
Compute the time to coalescence (in seconds) as a function of frequency in `Hz`. Used when including the Earth motion effect.

We use the expression in `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_ eq. (3.8b).

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency grid for a single event on which the time to coalescence will be computed, in `Hz`
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (array) time to coalescence for the chosen event evaluated on the frequency grid, in seconds.

"""
function _tau_star(model::Model, f, mc, eta; GMsun_over_c3 = uc.GMsun_over_c3)


    Mtot_sec = mc * GMsun_over_c3 / (eta^(0.6))
    v = (pi * Mtot_sec .* f) .^ (1.0 / 3.0)
    eta2 = eta * eta

    OverallFac = @. 5.0 / 256 * Mtot_sec / (eta * (v .^ 8.0))
    v2 = v .* v
    v3 = v2 .* v
    v4 = v3 .* v
    v5 = v4 .* v
    v6 = v5 .* v

    t05 =
        1.0 .+ (743.0 / 252.0 + 11.0 / 3.0 * eta) .* v2 .- 32.0 / 5.0 * pi .* v3 .+
        (3058673.0 / 508032.0 + 5429.0 / 504.0 * eta + 617.0 / 72.0 * eta2) .* v4 .-
        (7729.0 / 252.0 - 13.0 / 3.0 * eta) * pi .* v5
    t6 = 
        @. (
            -10052469856691.0 / 23471078400.0 +
            128.0 / 3.0 * pi^2 +
            6848.0 / 105.0 * MathConstants.eulergamma +
            (3147553127.0 / 3048192.0 - 451.0 / 12.0 * pi^2) * eta -
            15211.0 / 1728.0 * eta2 +
            25565.0 / 1296.0 * eta2 * eta +
            3424.0 / 105.0 * log(16.0 .* v2)
        ) * (v6)
    t7 = 
        @. (-15419335.0 / 127008.0 - 75703.0 / 756.0 * eta + 14809.0 / 378.0 * eta2) * pi .*
        v6 .* v

    return OverallFac .* (t05 .+ t6 .+ t7)
end

function _fcut(model::Model, mc, eta, Lambda1, Lambda2; fcutPar = 0.2, GMsun_over_c3 = uc.GMsun_over_c3)
    return _fcut(model, mc, eta, fcutPar = fcutPar, GMsun_over_c3 = GMsun_over_c3)
end

"""
Compute the cut frequency of the waveform as a function of the events parameters, in `Hz`.
Valid for TaylorF2, PhenomD, PhenomHM. Not for PhenomD_NRTidal.
We use the expression in `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_ eq. (3.8b).

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.

#### Optional arguments:
-  `fcutPar`: The cut frequency factor of the waveform as an adimensional frequency (Mf). Default is 0.2.

#### Return:
-  (float) Cut frequency of the waveform, in Hz.

"""
function _fcut(model::Model, mc, eta; fcutPar = 0.2, GMsun_over_c3 = uc.GMsun_over_c3)

    if typeof(model)== PhenomD_NRTidal
        error("You need to provide also Lambda1 and Lambda2 since fcut depends on them")
    end
    if typeof(model) == TaylorF2
        return 1. / ( 6 * pi * sqrt(6.) * GMsun_over_c3 )
    end

    return fcutPar / (mc * GMsun_over_c3 / (eta^(0.6)))

end

"""
useful functions for hphc()
"""
function _spinWeighted_SphericalHarmonic(theta, modes)
    # Taken from arXiv:0709.0093v3 eq. (II.7), (II.8) and LALSimulation for the s=-2 case and up to l=4.
    # We assume already phi=0 and s=-2 to simplify the function

    Ylm = @. ifelse(
        modes == 21,
        sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 + cos(theta)),
        ifelse(
            modes == 22,
            sqrt(5.0 / (64.0 * pi)) * (1.0 + cos(theta)) * (1.0 + cos(theta)),
            ifelse(
                modes == 32,
                sqrt(7.0 / pi) * (cos(theta * 0.5)^4) * (-2.0 + 3.0 * cos(theta)) * 0.5,
                ifelse(
                    modes == 33,
                    -sqrt(21.0 / (2.0 * pi)) * (cos(theta * 0.5)^5) * sin(theta * 0.5),
                    ifelse(
                        modes == 43,
                        -3.0 *
                        sqrt(7.0 / (2.0 * pi)) *
                        (cos(theta * 0.5)^5) *
                        (-1.0 + 2.0 * cos(theta)) *
                        sin(theta * 0.5),
                        3.0 * sqrt(7.0 / pi) * (cos(theta * 0.5)^6) * (sin(theta * 0.5)^2),
                    ),
                ),
            ),
        ),
    )
    Ylminm = @. ifelse(
        modes == 21,
        sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 - cos(theta)),
        ifelse(
            modes == 22,
            sqrt(5.0 / (64.0 * pi)) * (1.0 - cos(theta)) * (1.0 - cos(theta)),
            ifelse(
                modes == 32,
                sqrt(7.0 / (4.0 * pi)) *
                (2.0 + 3.0 * cos(theta)) *
                ((sin(theta * 0.5))^(4.0)),
                ifelse(
                    modes == 33,
                    sqrt(21.0 / (2.0 * pi)) * cos(theta * 0.5) * ((sin(theta * 0.5))^(5.0)),
                    ifelse(
                        modes == 43,
                        3.0 *
                        sqrt(7.0 / (2.0 * pi)) *
                        cos(theta * 0.5) *
                        (1.0 + 2.0 * cos(theta)) *
                        ((sin(theta * 0.5))^5.0),
                        3.0 *
                        sqrt(7.0 / pi) *
                        (cos(theta * 0.5) * cos(theta * 0.5)) *
                        ((sin(theta * 0.5))^6.0),
                    ),
                ),
            ),
        ),
    )

    return Ylm, Ylminm
end

end
