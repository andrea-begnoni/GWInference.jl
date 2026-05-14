####        GW POPULATION CATALOG       ####
module catalog # this is the name of this module
include("utils.jl") # include the utils module
import .UtilsAndConstants as uc # this is the module of utilities and constants

using Distributions # these are the imported Julia modules
using QuadGK
using Integrals
using HDF5
using Random
using Dates
using DelimitedFiles
using Interpolations
using Unitful
using UnitfulAstro      # Adds units of measurement commonly used in astrophysics
using Trapz


export GenerateCatalog, ReadCatalog, get_dL, ReadHyperparam # these are the functions that will be exported from this module

"""
Hubble constant as a function of redshift, results in km/s 1/Mpc

"""
function get_H_z(z, H0, Omega0_m, Omega0_Lambda)
    H_z = H0 * (Omega0_m * (1 .+ z) .^ 3 .+ Omega0_Lambda) .^ 0.5 # km/s 1/Mpc
    return H_z
end

"""
Luminosity distance as a function of redshift, results in Mpc. clight in km/s

"""
function get_dL(z, clight= uc.clight_kms, H0 = uc.H0, Omega0_m = uc.Omega0_m, Omega0_Lambda = uc.Omega0_Lambda)
    dL = zeros(length(z))
    jj = 1
    for ii in z
        dL[jj] =
            quadgk(k -> 1.0 ./ get_H_z(k, H0, Omega0_m, Omega0_Lambda), 0.0, ii)[1] *
            clight *
            (1.0 + ii)
        jj += 1
    end
    return dL
end

"""
Comoving distance as a function of redshift, results in Mpc. clight in km/s
"""
function get_comoving_distance(z, clight, H0, Omega0_m, Omega0_Lambda)
    chi = zeros(length(z))
    jj = 1
    for ii in z
        chi[jj] =
            quadgk(k -> 1.0 ./ get_H_z(k, H0, Omega0_m, Omega0_Lambda), 0.0, ii)[1] * clight
        jj += 1
    end
    return chi
end

"""
Comoving volume element as a function of redshift, results in Mpc^3. clight in km/s

"""
function get_dV_dz(z, clight, H0, Omega0_m, Omega0_Lambda)
    H_z = get_H_z(z, H0, Omega0_m, Omega0_Lambda)
    chi = quadgk(k -> 1.0 ./ get_H_z(k, H0, Omega0_m, Omega0_Lambda), 0.0, z)[1] * clight
    dV_dz = 4 * pi .* chi .^ 2 .* clight ./ H_z # Mpc^3
    return dV_dz
end

"""

Black hole mass function for NSBH, results in Mpc^-3

"""
function BHmass_NSBH(M_BH, a_1, b_1, a_2, b_2, a_3, b_3)
    f =
        (
            1 / (a_1 * exp(-b_1 * M_BH) + a_2 * exp(-b_2 * M_BH)) +
            1 / (a_3 * exp(b_3 * M_BH))
        )^(-1)
    return f
end

"""
Helper function
"""
function smoothing(m_prime, delta_m)
    f = exp(delta_m / m_prime + delta_m / (m_prime - delta_m))
    return f
end

"""
Helper function
"""
function S(m, mass_min, delta_m)
    if m < mass_min
        S_el = 0.0
    elseif m > mass_min && m < mass_min + delta_m
        
        S_el = (smoothing(m - mass_min, delta_m) + 1)^(-1)
    else
        S_el = 1.0
    end
    return S_el
end

"""
First mass of a BBH system, result in Msun, GWTC-3 (deprecated)
    
"""
function BHmass_BBH_first_mass(
    M_BH,
    lambda_peak,
    mass_min,
    mass_max,
    mu_m,
    sigma_m,
    alpha,
    delta_m,
)
    f =
        (
            (1 - lambda_peak) *
            pdf(truncated(Pareto(alpha, mass_min), mass_min, mass_max), M_BH) +
            lambda_peak * pdf(Normal(mu_m, sigma_m), M_BH)
        ) * S(M_BH, mass_min, delta_m)

    return f
end

"""
Helper function
"""
function BHmass_BBH_q(q, mass_min, delta_m, beta_q, m_1)
    f = q^beta_q * S(q * m_1, mass_min, delta_m)   # pdf of q given m_1
    return f
end

"""
Helper function
"""
function normalization_BHmass_BBH_q(
    mass_min,
    delta_m,
    beta_q,
    m_1,
    BHmass_BBH_q,
)
    f = quadgk(
        q -> BHmass_BBH_q(q, mass_min, delta_m, beta_q, m_1),
        1e-2,
        1.0,
    )[1] # normalization of the pdf of q given m_1
    return f
end

"""
Function for the Star Formation Rate (SFR) as a function of redshift, results needs to be normalized
"""
function get_Madau_and_friends(z, alpha_z, beta_z, z_p, a_z, b_z, c_z, d_z, friend)
    if friend == "Dickinson"
        f = (1 .+ z)^alpha_z / (1 .+ ((1 .+ z) / (1 .+ z_p))^(alpha_z + beta_z))    # Madau Dickinson
    elseif friend == "Fragos"
        f = a_z * (1 .+ z)^b_z / (1 .+ (c_z * (1 + z))^d_z)     # Madau Fragos
    end
    return f
end

# function get_Madau_and_Fragos(z, a_z, b_z, c_z, d_z)
#     f = a_z * (1 .+ z)^b_z / (1 .+ (c_z * (1 + z))^d_z)     # Madau Fragos
#     return f
# end
# function get_Madau_and_Dickinson(z, alpha_z, beta_z, z_p, amplitude)
#     f = amplitude * (1 .+ z)^alpha_z / (1 .+ ((1 .+ z) / (1 .+ z_p))^(alpha_z + beta_z))    # Madau Dickinson 
#     return f
# end


"""
Rate of BBH/BNS/NSBH mergers as a function of redshift, results need to be normalized, DEPRECATED

"""
function Redshift(
    z,
    population,
    alpha_z,
    beta_z,
    z_p,
    a_z_fit,
    b_z_fit,
    c_z_fit,
    d_z_fit,
    H0,
    Omega0_m,
    Omega0_Lambda,
    clight,
    get_dV_dz,
    friend,
)

    Rate =
        get_Madau_and_friends(
            z,
            alpha_z,
            beta_z,
            z_p,
            a_z_fit,
            b_z_fit,
            c_z_fit,
            d_z_fit,
            friend,
        ) / (1 .+ z) * get_dV_dz(z, clight, H0, Omega0_m, Omega0_Lambda)

    return Rate
end



"""
Rate of BBH/BNS/NSBH mergers as a function of redshift, results need to be normalized. Simplified for GWTC-4, since we only consider Madau Fragos and the time delay is arbitrary

"""
function Redshift(
    z,
    merger_rate_td_interp,
    t_delay_min,
    population,
    a_z,
    b_z,
    c_z,
    d_z,
    H0,
    Omega0_m,
    Omega0_Lambda,
    clight,
)

    # @assert population == "BBH" "Should not be here"

    if t_delay_min == 0.0
        Rate =
        madua_fragos(
            z,
            a_z,
            b_z,
            c_z,
            d_z
        ) / (1 .+ z) * get_dV_dz(z, clight, H0, Omega0_m, Omega0_Lambda)
    else
        vol = get_dV_dz(z, clight, H0, Omega0_m, Omega0_Lambda) * u"Mpc^3"
        
        Rate = merger_rate_td_interp(z) / (1 .+ z) * vol
    end
    

    return Rate
end
# Define the rejection sampling algorithm
# """
# Rejection sampling algorithm for the first mass of a BBH system, results in Msun
# """
# function rejection_sampling_first_mass(
#     BHmass_BBH_first_mass,
#     lambda_peak,
#     mass_min,
#     mass_max,
#     mu_m,
#     sigma_m,
#     alpha,
#     delta_m,
#     S,
#     smoothing,
#     normalization,
#     x_support,
#     n_samples,
# )   # x_support = [x_in, x_fin]
#     samples = zeros(n_samples)
#     i = 1
#     x_in = x_support[1]
#     x_fin = x_support[2]
#     while i <= n_samples
#         x = rand(Uniform(x_in, x_fin))
#         y = rand() # choose a y value between 0 and 1
#         if y <
#            BHmass_BBH_first_mass(
#             x,
#             lambda_peak,
#             mass_min,
#             mass_max,
#             mu_m,
#             sigma_m,
#             alpha,
#             delta_m,
#             S,
#             smoothing,
#         ) / normalization[1]
#             samples[i] = x
#             i += 1
#         end
#     end
#     return samples
# end


"""
Rejection sampling algorithm for the second mass of a BBH system, results in Msun
"""
function rejection_sampling_second_mass(
    BHmass_BBH_q,
    mass_min,
    delta_m,
    beta_q,
    x_support,
    n_samples,
    m_1,
    normalization_BHmass_BBH_q,
)   # x_support = [x_in, x_fin]
    samples = zeros(n_samples)
    i = 1
    x_in = x_support[1]

    x_fin = x_support[2]

    while i <= n_samples
        m = m_1[i]
        norm = normalization_BHmass_BBH_q(
            mass_min,
            delta_m,
            beta_q,
            m,
            BHmass_BBH_q,
        )
        x = rand(Uniform(x_in, x_fin))
        y = rand() * BHmass_BBH_q(1.01, mass_min, delta_m, beta_q, m) / norm # choose a y value between 0 and 1



        if y < BHmass_BBH_q(x, mass_min, delta_m, beta_q, m) / norm
            samples[i] = x * m

            i += 1

        end
    end
    return samples
end

"""
Rejection sampling algorithm for the redshift
"""
function rejection_sampling_redshift(
    Redshift,
    population,
    alpha_z,
    beta_z,
    z_p,
    a_z_fit,
    b_z_fit,
    c_z_fit,
    d_z_fit,
    H0,
    Omega0_m,
    Omega0_Lambda,
    clight,
    friend,
    normalization,
    x_support,
    n_samples,
)   # x_support = [x_in, x_fin]
    samples = zeros(n_samples)
    i = 1
    x_in = x_support[1]
    x_fin = x_support[2]
    while i <= n_samples
        x = rand(Uniform(x_in, x_fin))
        y = rand() # choose a y value between 0 and 1

        if y <
           Redshift(
            x,
            population,
            alpha_z,
            beta_z,
            z_p,
            a_z_fit,
            b_z_fit,
            c_z_fit,
            d_z_fit,
            H0,
            Omega0_m,
            Omega0_Lambda,
            clight,
            friend,
        ) / normalization[1]

            samples[i] = x
            i += 1
        end
    end
    return samples
end


"""
Rejection sampling algorithm for the redshift, optimized for GWTC-4
"""
function rejection_sampling_redshift_v2(z_grid, merger_rate_td, t_delay_min, population::String,
    a_z,
    b_z,
    c_z,
    d_z,
    H0,
    Omega0_m,
    Omega0_Lambda,
    clight,
    x_support,
    n_samples,
)   
    samples = zeros(n_samples)
    i = 1
    x_in = x_support[1]
    x_fin = x_support[2]
    merger_rate_td_interp = LinearInterpolation(z_grid, merger_rate_td, extrapolation_bc=0.)

    merger_rate_proxy(z) = Redshift(z, merger_rate_td_interp, t_delay_min, population, a_z, b_z, c_z, d_z, H0, Omega0_m, Omega0_Lambda, clight)
    
    max_pop = maximum(merger_rate_proxy.(x_support[1]:0.01:x_support[2])) * 1.1 # add a 10% safety margin to ensure the proposal distribution is always above the target distribution
    
    while i <= n_samples
        x = rand(Uniform(x_in, x_fin))
        y_ = rand() # choose a y value between 0 and 1
        y = y_ * max_pop # scale y to the maximum value of the population distribution in the support range

        if y < merger_rate_proxy(x)
            samples[i] = x
            i += 1
        end
    end
    return samples
end


"""
Rejection sampling algorithm for the BH mass of a NSBH system, results in Msun
"""
function rejection_sampling_BHmass_NSBH(
    BHmass_NSBH,
    a_1,
    b_1,
    a_2,
    b_2,
    a_3,
    b_3,
    normalization,
    x_support,
    n_samples,
)   # x_support = [x_in, x_fin]
    samples = zeros(n_samples)
    i = 1
    x_in = x_support[1]
    x_fin = x_support[2]
    while i <= n_samples
        x = rand(Uniform(x_in, x_fin)) # choose an x value between 0 and 10
        y = rand() # choose a y value between 0 and 1

        if y < BHmass_NSBH(x, a_1, b_1, a_2, b_2, a_3, b_3) / normalization[1]

            samples[i] = x
            i += 1
        end
    end
    return samples
end

"""
Main function of this module, generates a catalog of events.

#### Input arguments:
-  `nEvents`: int, number of events to generate.
-  `population` : string, The population of the events, can be "BBH", "BNS" or "NSBH".

#### Optional arguments:
-  `time_delay_in_Myr` : float, The time delay between the formation of the binary and the merger, in Myr. The ones availables are 0.0, 10.0, 20.0, 50.0.
-  `seed_par` : int, seed for the random number generator.
-  `name_catalog` : string, name of file containing the catalog.
-  `local_rate` : float, local rate of events in Gpc^-3 yr^-1. It is only used to give an estimate of the total number of events in a year.
-  `EoS` : string, equation of state, can be "AP3", "APR4_EPP", "ENG", "SLy", "WFF1".

#### Outputs :
-  `chirp_mass_detector_frame` : array of floats, chirp mass in the detector frame.
-  `eta` : array of floats, symmetric mass ratio.
-  `chi1` : array of floats, dimensionless spin of the most massive object. (chi1 = chi1z, since chi1x = chi1y = 0)
-  `chi2` : array of floats, dimensionless spin of the least massive object. (chi2 = chi2z, since chi2x = chi2y = 0)
-  `dL` : array of floats, luminosity distance in Gpc.
-  `theta` : array of floats, polar angle in radians.
-  `phi` : array of floats, azimuthal angle in radians.
-  `iota` : array of floats, inclination angle in radians.
-  `psi` : array of floats, polarization angle in radians.
- `tcoal` : array of floats, coalescence time in fraction of a day.
-  `phiCoal` : array of floats, coalescence phase in radians.
- `Lambda1` : array of floats, tidal deformability of object 1.
- `Lambda2` : array of floats, tidal deformability of object 2.

#### Example:
```julia
GenerateCatalog(100, "BBH", time_delay_in_Myr = 10., seed_par = 1234, name_catalog = "catalog.h5", local_rate = 17.0)
```

The code also saves the redshift in the catalog, in case you want to use it for other purposes.



"""
function GenerateCatalog(nEvents::Int, population::String; time_delay_in_Myr = 20., seed_par = nothing, name_catalog = nothing, local_rate = nothing, EoS = "AP3")

    if seed_par === nothing
        seed = rand(1:10000)
    else
        seed = seed_par
    end

    # time_delay_in_Myr = 10.0
    # population = "BBH"
    n_samples = nEvents    


    if name_catalog === nothing
        mkpath("catalogs/")
        name_file = string(
            "catalogs/catalog_",
            population,
            "_n_",
            n_samples,
            "_td_",
            time_delay_in_Myr,
            "Myrs.h5",
        )
    else
        mkpath("catalogs/")
        path = pwd()*"/catalogs/"
        name_file = path*name_catalog
    end


    # NSBH mass
    a_1 = 1.04e11
    b_1 = 2.15
    a_2 = 800
    b_2 = 0.29
    a_3 = 0.00285
    b_3 = 1.686

    # BBH mass

    # Get the path to the directory of this file
    PACKAGE_DIR = @__DIR__

    # Go one step back in the path (from ""GWInference.jl/src" to "GWInference.jl")
    PARENT_DIR = dirname(PACKAGE_DIR)

    # Construct the path to the "useful_files" folder from the parent directory
    path_hyper_param = joinpath(PARENT_DIR, "useful_files/extracted_hyperparameters.h5")


    lam0, lam1, alpha_1, alpha_2, beta, break_mass, delta_m_1, delta_m_2, m_peak_1, m_peak_2, mlow_1, mlow_2, rate_z0, sigma_peak_1, sigma_peak_2, mmax, mu_chi, sigma_chi, mu_spin, sigma_spin, xi_spin  = ReadHyperparam(path_hyper_param, verbose = false);
    mass_dist = BrokenPowerLawGaussianMixture(
        lam0,          # λ₀
        lam1,          # λ₁  
        alpha_1,       # α₁
        alpha_2,       # α₂
        break_mass,    # break_mass
        m_peak_1,      # m_peak_1
        sigma_peak_1,  # σ_peak_1
        m_peak_2,      # m_peak_2
        sigma_peak_2,  # σ_peak_2
        mlow_1,        # m_low
        mmax,          # m_high
        delta_m_1      # δ_m
    )

    if local_rate === nothing
        # source arXiv.2508.18083
        if population == "BBH"
            local_rate = round(rate_z0, digits=4)
        elseif population == "BNS"
            local_rate = 89.
        elseif population == "NSBH"
            local_rate = 23.
        end
        println("No local rate provided, using default value of ", local_rate, " Gpc^-3 yr^-1")
    else
        println("Local rate of ", local_rate, " Gpc^-3 yr^-1, the total number of events in a year will be estimated accordingly")
        
        #local_rate = rate_z0 * u"Gpc^-3*yr^-1" 
    end


    # Redshift

    # Madau Fragos
    a_z = 0.01
    b_z = 2.6
    c_z = 0.3125
    d_z = 6.2

    # Physics constant

    clight = uc.clight_kms  #km/s
    Omega0_m = uc.Omega0_m
    Omega0_Lambda = 1 - Omega0_m
    H0 = uc.H0 # km/s 1/Mpc
    H0_yr = 2.25e-18 * 3.154e7 # 1/yr

    if time_delay_in_Myr === nothing
        println("No time delay")
        t_delay_min = 0.0
    elseif time_delay_in_Myr !== nothing
        println("Time delay of ", time_delay_in_Myr, " Myr")
        t_delay_min = time_delay_in_Myr * u"Myr"
        
    end
        




    Random.seed!(seed)
    theta = acos.(rand(Uniform(-1, 1), n_samples))
    phi = rand(Uniform(0, 2 * pi), n_samples)
    iota = acos.(rand(Uniform(-1, 1), n_samples))
    psi = rand(Uniform(0, pi), n_samples)
    phiCoal = rand(Uniform(0, 2 * pi), n_samples)
    tcoal = rand(Uniform(0, 1), n_samples)  #fraction of a day, same as GWFAST. 
    #No need to consider overlaps now, all the events are independent and we only consider earth rotation around its axis (no needs to consider month and year)

    normalization_redshift = 0.
    m_NS_max = 2.5


    if population != "BBH"
        if EoS != "AP3" && EoS != "APR4_EPP" && EoS != "ENG" && EoS != "SLy" && EoS != "WFF1" && EoS != "SQM3" && EoS != "random"
            error("EoS not found, the ones available are AP3, APR4_EPP, ENG, SLy, WFF1, SQM3 and random (i.e., uniform distribution)")
        end

        # Get the path to the directory of this file
        PACKAGE_DIR = @__DIR__

        # Go one step back in the path (from ""GWInference.jl/src" to "GWInference.jl")
        PARENT_DIR = dirname(PACKAGE_DIR)

        # Construct the path to the "useful_files" folder from the parent directory
        EOS_DIR = joinpath(PARENT_DIR, "useful_files/EOS_table/")

        if EoS == "AP3"
            table_AP3 = readdlm(EOS_DIR*"eos_AP3_mass.dat")
            mass_AP3, Lambda_AP3 = table_AP3[:,1], table_AP3[:,2]
            Lambda_interp = linear_interpolation(mass_AP3, Lambda_AP3)
            m_NS_max = mass_AP3[end]
        elseif EoS == "APR4_EPP"
            table_APR4_EPP = readdlm(EOS_DIR*"eos_APR4_EPP_mass.dat")
            mass_APR4_EPP, Lambda_APR4_EPP = table_APR4_EPP[:,1], table_APR4_EPP[:,2]
            Lambda_interp = linear_interpolation(mass_APR4_EPP, Lambda_APR4_EPP)
            m_NS_max = mass_APR4_EPP[end]
        elseif EoS == "ENG"
            table_ENG = readdlm(EOS_DIR*"eos_ENG_mass.dat")
            mass_ENG, Lambda_ENG = table_ENG[:,1], table_ENG[:,2]
            Lambda_interp = linear_interpolation(mass_ENG, Lambda_ENG)
            m_NS_max = mass_ENG[end]
        elseif EoS == "SLy"
            table_SLy = readdlm(EOS_DIR*"eos_SLy_mass.dat")
            mass_SLy, Lambda_SLy = table_SLy[:,1], table_SLy[:,2]
            Lambda_interp = linear_interpolation(mass_SLy, Lambda_SLy)
            m_NS_max = mass_SLy[end]
        elseif EoS == "WFF1"
            table_WFF1 = readdlm(EOS_DIR*"eos_WFF1_mass.dat")
            mass_WFF1, Lambda_WFF1 = table_WFF1[:,1], table_WFF1[:,2]
            Lambda_interp = linear_interpolation(mass_WFF1, Lambda_WFF1)
            m_NS_max = mass_WFF1[end]
        elseif EoS == "SQM3"
            table_AP3 = readdlm(EOS_DIR*"eos_AP3_mass.dat")
            mass_AP3, Lambda_AP3 = table_AP3[:,1], table_AP3[:,2]
            m_NS_max = 2.
            mass_SQM3 = mass_AP3[mass_AP3 .<= m_NS_max]
            Lambda_SQM3 = Lambda_AP3[mass_AP3 .<= m_NS_max]
            Lambda_interp = linear_interpolation(mass_SQM3, Lambda_SQM3)
        elseif EoS == "random"
            Lambda_interp = x -> rand(Uniform(0.,2e3), length(x))
        end        
    end

    if population == "BBH"

        Lambda1 = zeros(n_samples)
        Lambda2 = zeros(n_samples)

        #spins

        dist_chi = truncated(Normal(mu_chi, sigma_chi), 0.0, 1.0)
        chi1 = rand(dist_chi,n_samples)
        chi2 = rand(dist_chi,n_samples)

        dist1 = truncated(Normal(mu_spin, sigma_spin), -1.0, 1.0)
        dist2 = Uniform(-1.0, 1.0)
        dist_theta_spin = MixtureModel([dist1, dist2], [xi_spin, 1 - xi_spin])

        cos_theta_spin_1 = rand(dist_theta_spin,n_samples)
        cos_theta_spin_2 = rand(dist_theta_spin,n_samples)

        chiz1 = chi1 .* cos_theta_spin_1
        chiz2 = chi2 .* cos_theta_spin_2
        

        m_1 = rand(mass_dist, n_samples)

        m_2 = rejection_sampling_second_mass(
            BHmass_BBH_q,
            mlow_2,
            delta_m_2,
            beta,
            [1e-2, 1.0],
            n_samples,
            10. .* ones(n_samples),
            normalization_BHmass_BBH_q,
        );

        k_GWTC4= 3.2 # did not find it in the hyperparameter file, but it is mentioned in the paper
        zref = 0.2
        rate_zref = rate_z0 * (1 + zref)^k_GWTC4
        LVK_rate_zref = rate_zref * u"Gpc^-3*yr^-1" 


    elseif population == "BNS"
        chiz1 = rand(Uniform(-0.05, 0.05), n_samples)
        chiz2 = rand(Uniform(-0.05, 0.05), n_samples)


        m_1 = rand(Uniform(1, m_NS_max), n_samples)
        m_2 = rand(Uniform(1, m_NS_max), n_samples)
        for i in 1:n_samples
            if m_1[i] < m_2[i]
                tmp = m_1[i]
                m_1[i] = m_2[i]
                m_2[i] = tmp
            end
        end

        Lambda1 = Lambda_interp.(m_1)
        Lambda2 = Lambda_interp.(m_2)

        zref = 0.0
        LVK_rate_zref = local_rate * u"Gpc^-3*yr^-1"

    elseif population == "NSBH"
        chiz1 = rand(Uniform(-0.05, 0.05), n_samples)
        chiz2 = rand(Normal(0.0, 0.15), n_samples)

        normalization_BHmass_NSBH =
            quadgk(M_BH -> BHmass_NSBH(M_BH, a_1, b_1, a_2, b_2, a_3, b_3), 2.0, 25.0)[1]
        m_1 = rejection_sampling_BHmass_NSBH(
            BHmass_NSBH,
            a_1,
            b_1,
            a_2,
            b_2,
            a_3,
            b_3,
            normalization_BHmass_NSBH,
            [2.0, 25.0],
            n_samples,
        )
        m_2 = rand(Normal(1.33, 0.09), n_samples)

        Lambda1 = Lambda_interp.(m_2)
        Lambda2 = zeros(n_samples)

        zref = 0.0
        LVK_rate_zref = local_rate * u"Gpc^-3*yr^-1"

    end

    z_grid, merger_rate_td = merger_rate_td_min(t_delay_min, zref, LVK_rate_zref, a_z, b_z, c_z, d_z);
    merger_rate_td_interp = LinearInterpolation(z_grid, merger_rate_td, extrapolation_bc=0.)

    z = rejection_sampling_redshift_v2(z_grid, merger_rate_td, t_delay_min, population, a_z, b_z, c_z, d_z, H0, Omega0_m, Omega0_Lambda, clight, [0.0, 10.0], n_samples)

    fproxy(z_grid) = Redshift(z_grid, merger_rate_td_interp, t_delay_min, population, a_z, b_z, c_z, d_z, H0, Omega0_m, Omega0_Lambda, clight)
    merger_rate_td_not_normalized = fproxy.(z_grid)
    normalization_redshift = trapz(z_grid, ustrip.(merger_rate_td_not_normalized))  # not used for normalization, just to check the total number of events in a year


    ## number of events in an year
    # integrate the rate over redshift
    total_number_sources_yr = normalization_redshift / 1e9 # to convert the normalization factor from Mpc^3 to Gpc^3
    # this number is not used in the code, it is just to have an estimate of the total number of sources in a year



    chirp_mass = (m_1 .* m_2) .^ (3 / 5) ./ (m_1 .+ m_2) .^ (1 / 5)
    chirp_mass_detector_frame = chirp_mass .* (1 .+ z)
    eta = (m_1 .* m_2) ./ (m_1 .+ m_2) .^ 2
    if unit(H0) == u"km/s/Mpc"
        H0_kms_Mpc = ustrip(H0) 
        dL = get_dL(z, clight, H0_kms_Mpc, Omega0_m, Omega0_Lambda) ./ 1e3 # Gpc
    else 
        error("H0 should be in km/s/Mpc")
    end
    date = Dates.now()
    date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
    println("Name of the catalog: ", name_file)
    h5open(name_file, "w") do file
        HDF5.attributes(file)["format"] = "GWJulia"
        HDF5.attributes(file)["number_events"] = nEvents
        HDF5.attributes(file)["seed"] = seed
        HDF5.attributes(file)["population"] = population
        HDF5.attributes(file)["time_delay_in_Myrs"] = time_delay_in_Myr
        HDF5.attributes(file)["SFR"] = "Madau&Fragos"
        HDF5.attributes(file)["total_number_sources_yr"] = Int(round(total_number_sources_yr, digits=0))
        HDF5.attributes(file)["local_rate"] = local_rate
        HDF5.attributes(file)["date"] = date_format
        if population != "BBH"
            HDF5.attributes(file)["EoS"] = EoS
        else
            HDF5.attributes(file)["EoS"] = "none"
        end


        write(file, "mc", chirp_mass_detector_frame)
        write(file, "eta", eta)
        write(file, "chi1", chiz1)  # chi1 = chi1z, since chi1x = chi1y = 0
        write(file, "chi2", chiz2)  # chi2 = chi2z, since chi2x = chi2y = 0
        write(file, "dL", dL)
        write(file, "theta", theta)
        write(file, "phi", phi)
        write(file, "iota", iota)
        write(file, "psi", psi)
        write(file, "phiCoal", phiCoal)
        write(file, "tcoal", tcoal)
        write(file, "Lambda1", Lambda1)
        write(file, "Lambda2", Lambda2)
        write(file, "z", z)
    end

    return chirp_mass_detector_frame,
    eta,
    chiz1,
    chiz2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    Lambda1,
    Lambda2
    

end


"""
Read the catalog generated with the function GenerateCatalog.

#### Input arguments:
-  `name_file` : string, name of the file containing the catalog.

#### Optional arguments:
-  `folder` : string, folder where the catalog is stored.
-  `redshift` : bool, if true, the redshift is read from the catalog.

#### Outputs:
-  `mc` : array of floats, chirp mass in the detector frame.
-  `eta` : array of floats, symmetric mass ratio.
-  `chi1` : array of floats, dimensionless spin of the most massive object. (chi1 = chi1z, since chi1x = chi1y = 0)
-  `chi2` : array of floats, dimensionless spin of the least massive object. (chi2 = chi2z, since chi2x = chi2y = 0)
-  `dL` : array of floats, luminosity distance in Gpc.  
-  `theta` : array of floats, polar angle in radians.
-  `phi` : array of floats, azimuthal angle in radians.
-  `iota` : array of floats, inclination angle in radians.
-  `psi` : array of floats, polarization angle in radians.
-  `tcoal` : array of floats, coalescence time in fraction of a day.
-  `phiCoal` : array of floats, coalescence phase in radians.
-  `Lambda1` : array of floats, tidal deformability of object 1.
-  `Lambda2` : array of floats, tidal deformability of object 2.
-  `z` : array of floats, redshift, returned only if redshift = true.

"""
function ReadCatalog(name_file; folder= "catalogs/", redshift = false)

    
    mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2, z = h5open(folder*name_file, "r") do file
        ### see HDF5.attributes
        println("attributes: ", keys(HDF5.attributes(file)))
        attributes_keys = keys(HDF5.attributes(file))
        for i in 1:length(keys(HDF5.attributes(file)))
            println(attributes_keys[i], ": ", read(HDF5.attributes(file)[attributes_keys[i]]))
        end
        println("Parameters: ", keys(file))
        nEvents = read(HDF5.attributes(file)["number_events"])
        population = read(HDF5.attributes(file)["population"])
        if population == "BBH"
            Lambda1 = zeros(nEvents)
            Lambda2 = Lambda1
        elseif population == "BNS"
            Lambda1 = read(file, "Lambda1")
            Lambda2 = read(file, "Lambda2")
        elseif population == "NSBH"
            Lambda1 = read(file, "Lambda1")
            Lambda2 = zeros(nEvents)
        else
            println("Population not recognized, the ones available are BBH, BNS, NSBH")
        end
        chi1 = read(file, "chi1")
        chi2 = read(file, "chi2")
        mc = read(file, "mc")
        eta = read(file, "eta")
        dL = read(file, "dL")
        if redshift == true
            z = read(file, "z")
        else
            z = zeros(nEvents)
        end
        theta = read(file, "theta")
        phi = read(file, "phi")
        iota = read(file, "iota")
        psi = read(file, "psi")
        phiCoal = read(file, "phiCoal")
        tcoal = read(file, "tcoal")
        return     mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2, z

    end
    if redshift == true
        return mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2, z
    else 
        return mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2
    end
end


function ReadHyperparam(name_file; verbose = true)

    params = h5open(name_file, "r") do file
        ### see HDF5.attributes
        if verbose == true
            println("HDF5.attributes: ", keys(HDF5.attributes(file)))
        end
        attributes_keys = keys(HDF5.attributes(file))
        for i in 1:length(keys(HDF5.attributes(file)))
            if verbose == true
                println(attributes_keys[i], ": ", read(HDF5.attributes(file)[attributes_keys[i]]))
            end
        end
        if verbose == true
            println("Groups: ", keys(file))
        end 
        groups_keys = keys(file)
        for i in 1:length(keys(file))
            if verbose == true
                println(groups_keys[i], ": ", keys(file[groups_keys[i]]))
            end
        end
        lam0 = read(file["individual_parameters/lambda_0"])
        lam1 = read(file["individual_parameters/lambda_1"])
        alpha_1 = read(file["individual_parameters/alpha_1"])
        alpha_2 = read(file["individual_parameters/alpha_2"])
        beta = read(file["individual_parameters/beta"])
        break_mass = read(file["individual_parameters/break_mass"])
        delta_m_1 = read(file["individual_parameters/delta_m_1"])
        delta_m_2 = read(file["individual_parameters/delta_m_2"])
        m_peak_1 = read(file["individual_parameters/m_peak_1"])
        m_peak_2 = read(file["individual_parameters/m_peak_2"])
        mlow_1 = read(file["individual_parameters/mlow_1"])
        mlow_2 = read(file["individual_parameters/mlow_2"])
        rate = read(file["individual_parameters/rate"])
        sigma_peak_1 = read(file["individual_parameters/sigma_peak_1"])
        sigma_peak_2 = read(file["individual_parameters/sigma_peak_2"])
        mmax = read(file["individual_parameters/mmax"]) 
        mu_chi = read(file["individual_parameters/mu_chi"])
        sigma_chi = read(file["individual_parameters/sigma_chi"])
        mu_spin = read(file["individual_parameters/mu_spin"])
        sigma_spin = read(file["individual_parameters/sigma_spin"])
        xi_spin = read(file["individual_parameters/xi_spin"])
        return [lam0, lam1, alpha_1, alpha_2, beta, break_mass, delta_m_1, delta_m_2, m_peak_1, m_peak_2, mlow_1, mlow_2, rate, sigma_peak_1, sigma_peak_2, mmax, mu_chi, sigma_chi, mu_spin, sigma_spin, xi_spin]
    end

end


# custom distribution type for the broken power law + Gaussian mixture

struct BrokenPowerLawGaussianMixture <: ContinuousUnivariateDistribution
    λ₀::Float64        # lambda_0 (broken power law weight)
    λ₁::Float64        # lambda_1 (first Gaussian weight) 
    α₁::Float64        # alpha_1 (first power law index)
    α₂::Float64        # alpha_2 (second power law index)
    break_mass::Float64 # break mass
    m_peak_1::Float64   # first Gaussian peak
    σ_peak_1::Float64   # first Gaussian std
    m_peak_2::Float64   # second Gaussian peak  
    σ_peak_2::Float64   # second Gaussian std
    m_low::Float64      # minimum mass
    m_high::Float64     # maximum mass
    δ_m::Float64        # smoothing parameter
    norm_constant::Float64  # normalization constant for the full PDF
    
    function BrokenPowerLawGaussianMixture(λ₀, λ₁, α₁, α₂, break_mass, m_peak_1, σ_peak_1, 
                                          m_peak_2, σ_peak_2, m_low, m_high, δ_m)
        # Validation
        @assert 0 ≤ λ₀ ≤ 1 "λ₀ must be between 0 and 1"
        @assert 0 ≤ λ₁ ≤ 1 "λ₁ must be between 0 and 1" 
        @assert λ₀ + λ₁ ≤ 1 "λ₀ + λ₁ must be ≤ 1"
        @assert m_low < m_high "m_low must be < m_high"
        @assert m_low < break_mass < m_high "break_mass must be between m_low and m_high"
        
        # Compute normalization constant
        norm_const = compute_normalization(λ₀, λ₁, α₁, α₂, break_mass, m_peak_1, σ_peak_1, 
                                           m_peak_2, σ_peak_2, m_low, m_high, δ_m)
        
        new(λ₀, λ₁, α₁, α₂, break_mass, m_peak_1, σ_peak_1, m_peak_2, σ_peak_2, m_low, m_high, δ_m, norm_const)
    end
end

# Define the support (domain) of the distribution
Distributions.support(d::BrokenPowerLawGaussianMixture) = Distributions.RealInterval(d.m_low, d.m_high)



# Helper function to compute normalization constant

function compute_normalization(λ₀, λ₁, α₁, α₂, break_mass, m_peak_1, σ_peak_1, 
                               m_peak_2, σ_peak_2, m_low, m_high, δ_m)
    # Unnormalized PDF function
    function unnormalized_pdf(x)
        if x < m_low || x > m_high
            return 0.0
        end
        
        # Broken power law component
        if x < break_mass
            p_BP = (x / break_mass)^(-α₁)
        else
            p_BP = (x / break_mass)^(-α₂)
        end
        
        # Normalization constant for broken power law only
        norm_const_BP = break_mass * (
            (1.0 - (m_low / break_mass)^(1.0 - α₁)) / (1.0 - α₁) + 
            ((m_high / break_mass)^(1.0 - α₂) - 1.0) / (1.0 - α₂)
        )
        p_BP_normalized = p_BP / norm_const_BP
        
        # Gaussian components  
        gauss1 = pdf(Normal(m_peak_1, σ_peak_1), x)
        gauss2 = pdf(Normal(m_peak_2, σ_peak_2), x)
        
        # Mixture
        mixture = λ₀ * p_BP_normalized + λ₁ * gauss1 + (1 - λ₀ - λ₁) * gauss2
        
        # Apply smoothing function
        return mixture * S(x, m_low, δ_m)
    end
    
    # Compute integral
    integral, _ = quadgk(unnormalized_pdf, m_low, m_high, rtol=1e-8)
    return integral
end

# Define the probability density function
function Distributions.pdf(d::BrokenPowerLawGaussianMixture, x::Real)
    if x < d.m_low || x > d.m_high
        return 0.0
    end
    
    # Broken power law component
    if x < d.break_mass
        p_BP = (x / d.break_mass)^(-d.α₁)
    else
        p_BP = (x / d.break_mass)^(-d.α₂)
    end
    
    # Normalization constant for broken power law only
    norm_const_BP = d.break_mass * (
        (1.0 - (d.m_low / d.break_mass)^(1.0 - d.α₁)) / (1.0 - d.α₁) + 
        ((d.m_high / d.break_mass)^(1.0 - d.α₂) - 1.0) / (1.0 - d.α₂)
    )
    p_BP_normalized = p_BP / norm_const_BP
    
    # Gaussian components  
    gauss1 = pdf(Normal(d.m_peak_1, d.σ_peak_1), x)
    gauss2 = pdf(Normal(d.m_peak_2, d.σ_peak_2), x)
    
    # Mixture
    mixture = d.λ₀ * p_BP_normalized + d.λ₁ * gauss1 + (1 - d.λ₀ - d.λ₁) * gauss2
    
    # Apply smoothing function and normalize
    return (mixture * S(x, d.m_low, d.δ_m)) / d.norm_constant
end

# Define sampling using rejection sampling (since this is a complex mixture)
function Distributions.rand(rng::AbstractRNG, d::BrokenPowerLawGaussianMixture)
    # Find approximate maximum of the pdf for rejection sampling
    x_test = range(d.m_low, d.m_high, length=1000)
    pdf_max = maximum([pdf(d, x) for x in x_test]) * 1.1  # Add small buffer
    
    # Rejection sampling
    while true
        x = rand(rng, Uniform(d.m_low, d.m_high))
        u = rand(rng) * pdf_max
        if u ≤ pdf(d, x)
            return x
        end
    end
end

# Convenience method for default RNG
Distributions.rand(d::BrokenPowerLawGaussianMixture) = rand(Random.GLOBAL_RNG, d)
    



### ALESSIO


#export merger_rate_powerlaw, hubble_parameter, delaytime, HubbleParameter, madau_fragos, HubbleParameter



# Define an object which stores the value of H(z) and the z at which it was computed (to perform checks)
# Assumes that time is expressed in megayears
struct HubbleParameter
    unit::Unitful.Unit
    value::AbstractVector   # Value of the hubble parameter 
    range::AbstractVector   # Redshift at which the above values of the Hubble parameter were computed
end



# Common approximation for the merger rate at low redshifts using a power law
#
# ---Parameters---
# z: redshift
# k: the power law index
# norm: a normalization factor

function merger_rate_powerlaw(z::AbstractVector, norm, k=uc.k_GWTC4)
    return norm.*(1 .+z).^k
end

# When one is approximating the merger rate with a powelaw, one may want to specify the desired value of the merger rate at a certain redshift different from 0.
#
# ---Parameters---
# z_eval: the value of redshift at which merger_rate is referred.
# k: the exponent to the powerlaw to be normalized
# merger_rate: value of the merger rate at z_eval which is to be imposed
# ---Returns---
# normalization: the normalization constant by which merger_rate_powerlaw needs to be multiplied

function get_normalization_constant(z_eval::Float64, k::Float64, merger_rate::typeof(1.0u"Gpc^-3*yr^-1"))
    normalization = merger_rate/(1 + z_eval)^k
    return normalization
end



# Computes the value of a primitive of the integral in function delaytime, so that one of the bounds of integration can be determined knowing t_delay and the other 
# bound. Default values of the cosmological parameters from Planck 2018 (table 2, arxiv: 1807.06209)
#
# ---Parameters---
# z: the redshift at which the value of the primitive is to be computed
# Omega0_m: the matter density in the local universe
# Omega0_Lambda: the cosmological constant density
# H0: the hubble parameter in the local universe
# ---Returns---
# primitive: the value of the primitive at z.

function primitive_time_delay(z, Omega0_m = uc.Omega0_m, Omega0_Lambda = uc.Omega0_Lambda, H0=uc.H0)
    k = sqrt(Omega0_m/Omega0_Lambda)*(1 .+z).^(3/2)    
    primitive = -2/(3*H0*sqrt(Omega0_Lambda)) * acsch.(k)
    return primitive
end

# Computes the formation redshift once the merger redshift and the time delay is known. Essentially, if one has a primitive of the integral which yields the time delay
# corresponding to the interval (zmerger, zformation). Default values of the cosmological parameters are from Planck 2018 (table 2, arxiv: 1807.06209)
#
# ---Parameters---
# z_merger: redshift of the merger
# t_delay: the time delay between the merger redshift and the formation redshift
# Omega0_m: the matter density in the local universe
# Omega0_Lambda: the cosmological constant density
# H0: the Hubble parameter in the local universe
# ---Returns---
# z_formation: the formation redshift

function get_formation_redshift(z_merger, t_delay, Omega0_m=uc.Omega0_m, Omega0_Lambda=uc.Omega0_Lambda, H0=uc.H0)
    var1 = csch.((-3/2)*H0*sqrt(Omega0_Lambda)*(t_delay .+ primitive_time_delay(z_merger, Omega0_m, Omega0_Lambda, H0)))
    zplusone32 = var1.*sqrt(Omega0_Lambda/Omega0_m)

    zplusone = zplusone32.^(2/3)

    z_formation = zplusone .- 1
    
    return z_formation
end


# Computes the hubble parameter in a flat ΛCDM universe. Default parameters are from the Planck 2018 Cosmology, table 2, arxiv: 1807.06209
#
# ---Parameters---
# H0: the hubble constant today
# Omega0_m: the matter density
# Omega0_Lambda: the cosmological constant density

function hubble_parameter(z, H0 = uc.H0, Omega0_m = uc.Omega0_m, Omega0_Lambda=uc.Omega0_Lambda)

    H = H0 .*sqrt.(Omega0_m .*(1 .+z).^3 .+ Omega0_Lambda)      # Friedmann's equation for a flat universe

    #output = HubbleParameter(hubble_parameter, z) # This line should hopefully build an object of type HubbleParameter to return (at least this is what i understand from
    # the tutorial)
    return H
end

# This is a parametrization of the star formation rate density [M_sun Mpc^-3 yr^-1]. This function is a transcription of equation 9 of 2303.10693v2
# Default values for the parameters are taken from Madau and Fragos (2017, 1606.07887). 
#
# ---Parameters---
# z: redshift
# a, b, c, d: parameters of the Madau-Fragos model
# ---Returns---
# ψ: the star formation rate density

function madau_fragos(z, a = 0.01u"Msun/Mpc^3/yr", b = 2.6, c = 1/3.2, d = 6.2)
    ψ = a* ((1 .+z).^b) ./ (1 .+(c.*(1 .+z)).^d )
    return ψ
end

# An alternative parametrization of the star formation rate density, which is found for instance in 2304.06368v2. Note that differently from the cited work, we 
# do not include a multiplication constant (so by default the returned ψ is ADIMENSIONAL)
# Default values for the parameters are taken from 2304.06368v2. 
#
# ---Parameters---
# z: redshift
# k, z_peak, r: parameters of the alternative Madau-Fragos model
# ---Returns---
# ψ: the star formation rate density

function madau_fragos_alternative(z, k = 2.6, z_peak=2.04, r=3.6)
    ψ = ((1 .+z).^k) ./ (1 .+ k/r.*((1 .+z)/(1+z_peak)).^(k+r))
    return ψ
end


# The delay time probability density function
#
# ---Parameters---
# t_delay: the delay time
# t_delay_max: the maximum delay time between formation and merger
# t_delay_min: the minumum delay time between formation and merger
# ---Returns---
# p: the probability of havin a certain t_delay


function time_delay_pdf(t_delay, t_delay_max, t_delay_min)
    p = (t_delay .* log(t_delay_max/t_delay_min))^(-1) 

    if t_delay < t_delay_min || t_delay > t_delay_max
        p = 0
    end
    
    return p
end

# Computes the age of the universe at redshift z for given cosmological parameters (Omega0_Lambda, H0, Omega0_m). Does not necessarily require that universe is flat. 
# Default values: Planck 2018 (table 2, arxiv: 1807.06209). Equation from M. Cignoni, lecture notes (2024). 
#
# ---Parameters---
# z: the redshift at which the age of the universe will be computed
# Omega0_Lambda: the cosmological constant density
# Omega0_m: the matter density at present time
# H0: the Hubble parameter in the local universe. 
# ---Returns---
# div: the age of the universe at redshift z expressed in Myr

function how_old_universe(z, Omega0_Lambda = uc.Omega0_Lambda, Omega0_m = uc.Omega0_m, H0=uc.H0)
    # recall a = 1/(1+z)
    domain=(0, 1/(1+z))
    Omega_0 = Omega0_m + Omega0_Lambda

    f(a, p) = sqrt(1/(Omega0_m*a^(-1) + Omega0_Lambda*a^2 + 1 - Omega_0))

    n = IntegralProblem(f, domain)
    solution = solve(n, QuadGKJL(), reltol=1e-8, abstol = 1e-8)    # Solve the integral problem. Arguments: the problem, the algorithm, tolerances in the numerical integration

    div = solution.u/H0
    div = uconvert(u"Myr", div) # Convert and express result in megayears
    return div

end

"""
Computes the merger rate density as a function of redshift, after marginalizing over the time delay between formation and merger. The merger rate density is normalized to match the LVK rate density at a reference redshift zref.
"""
function merger_rate_td_min(t_delay_min, zref, LVK_rate_zref, a, b, c, d; Omega0_Lambda = uc.Omega0_Lambda, Omega0_m = uc.Omega0_m, H0=uc.H0)

    r = range(0, stop = 10, length=401)  # redshifts, length = 101 is necessary so that there is a row of the dataframe whose redshift is exactly 0.2
    merger_rate_td = zeros(length(r)) # merger rate density as a function of redshift, after marginalizing over time delay

    for (i, z) in enumerate(r)

        age = how_old_universe(z, Omega0_Lambda, Omega0_m, H0) # Calculate the age of the universe at redshift z
        t_delay_max = min(age, 10u"Gyr")    # Maximum possible time delay is the age of the universe if it is less than 10 Gyr, else it is 10 Gyr

        if t_delay_max <= t_delay_min
            throw(ArgumentError("Minimum time delay cannot be bigger than maximum time delay"))
        end


        madau_fragos_proxy(z, t_delay) = madau_fragos(get_formation_redshift(z, t_delay, Omega0_m, Omega0_Lambda, H0), a, b, c, d)

        f(t_delay, z) = madau_fragos_proxy(z, t_delay)*time_delay_pdf(t_delay, 10u"Gyr", t_delay_min)

        solution = quadgk(t_delay -> f(t_delay, z), t_delay_min, t_delay_max)[1]

        merger_rate_td[i] = ustrip(solution)
    end

    ### normalize the merger rate density to match the LVK rate density at redshift 0.2
    zref_index = findfirst(r .== zref)
    merger_rate_td_normalized = merger_rate_td * (LVK_rate_zref / merger_rate_td[zref_index])

    return r, merger_rate_td_normalized
end



end
    

