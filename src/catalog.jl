####        GW POPULATION CATALOG       ####
module catalog # this is the name of this module

import ..UtilsAndConstants as uc # this is the module of utilities and constants

using Distributions # these are the imported Julia modules
using QuadGK
using Integrals
using HDF5
using Random
using Dates

export GenerateCatalog, ReadCatalog, get_dL # these are the functions that will be exported from this module

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
function get_comoving_distance(z, clight, get_H_z, H0, Omega0_m, Omega0_Lambda)
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
function get_dV_dz(z, clight, get_H_z, H0, Omega0_m, Omega0_Lambda)
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
function S(m, mass_min, delta_m, smoothing)
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
First mass of a BBH system, result in Msun
    
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
    S,
    smoothing,
)
    f =
        (
            (1 - lambda_peak) *
            pdf(truncated(Pareto(alpha, mass_min), mass_min, mass_max), M_BH) +
            lambda_peak * pdf(Normal(mu_m, sigma_m), M_BH)
        ) * S(M_BH, mass_min, delta_m, smoothing)

    return f
end

"""
Helper function
"""
function BHmass_BBH_q(q, mass_min, delta_m, beta_q, S, smoothing, m_1)
    f = q^beta_q * S(q * m_1, mass_min, delta_m, smoothing)   # pdf of q given m_1
    return f
end

"""
Helper function
"""
function normalization_BHmass_BBH_q(
    mass_min,
    delta_m,
    beta_q,
    S,
    smoothing,
    m_1,
    BHmass_BBH_q,
)
    f = quadgk(
        q -> BHmass_BBH_q(q, mass_min, delta_m, beta_q, S, smoothing, m_1),
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
Rate of BBH/BNS/NSBH mergers as a function of redshift, results need to be normalized

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
        ) / (1 .+ z) * get_dV_dz(z, clight, get_H_z, H0, Omega0_m, Omega0_Lambda)

    return Rate
end

# Define the rejection sampling algorithm
"""
Rejection sampling algorithm for the first mass of a BBH system, results in Msun
"""
function rejection_sampling_first_mass(
    BHmass_BBH_first_mass,
    lambda_peak,
    mass_min,
    mass_max,
    mu_m,
    sigma_m,
    alpha,
    delta_m,
    S,
    smoothing,
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
           BHmass_BBH_first_mass(
            x,
            lambda_peak,
            mass_min,
            mass_max,
            mu_m,
            sigma_m,
            alpha,
            delta_m,
            S,
            smoothing,
        ) / normalization[1]
            samples[i] = x
            i += 1
        end
    end
    return samples
end


"""
Rejection sampling algorithm for the second mass of a BBH system, results in Msun
"""
function rejection_sampling_second_mass(
    BHmass_BBH_q,
    mass_min,
    delta_m,
    beta_q,
    S,
    smoothing,
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
            S,
            smoothing,
            m,
            BHmass_BBH_q,
        )
        x = rand(Uniform(x_in, x_fin))
        y = rand() * BHmass_BBH_q(1.01, mass_min, delta_m, beta_q, S, smoothing, m) / norm # choose a y value between 0 and 1



        if y < BHmass_BBH_q(x, mass_min, delta_m, beta_q, S, smoothing, m) / norm
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
    get_dV_dz,
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
            get_dV_dz,
            friend,
        ) / normalization[1]

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
-  `SFR` : string, Star Formation Rate, can be "Madau&Dickinson" or "Madau&Fragos".
-  `name_catalog` : string, name of file containing the catalog.
-  `local_rate` : float, local rate of events in Gpc^-3 yr^-1. It is only used to give an estimate of the total number of events in a year.

#### Outputs (for the format "GWJulia"):
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
GenerateCatalog(100, "BBH", time_delay_in_Myr = 10., seed_par = 1234, SFR = "Madau&Dickinson", name_catalog = "catalog.h5", local_rate = 17.0)
```

The code also saves the redshift in the catalog, in case you want to use it for other purposes.



"""
function GenerateCatalog(nEvents::Int, population::String; time_delay_in_Myr = 10., seed_par = nothing, SFR = "Madau&Dickinson", name_catalog = nothing, local_rate = nothing, auto_save=true)

    if seed_par === nothing
        seed = rand(1:10000)
    else
        seed = seed_par
    end

    # time_delay_in_Myr = 10.0
    # population = "BBH"
    n_samples = nEvents    

    #friend = "Dickinson"
    if SFR != "Madau&Dickinson" && SFR != "Madau&Fragos"
        error("SFR not recognized, the ones available are Madau&Dickinson and Madau&Fragos")
    elseif SFR == "Madau&Dickinson"
        friend = "Dickinson"
    elseif SFR == "Madau&Fragos"
        friend = "Fragos"
    end

    if auto_save
        if name_catalog === nothing
            mkpath("catalogs/")
            name_file = string(
                "catalogs/catalog_",
                population,
                "_n_",
                n_samples,
                "_",
                SFR,
                "_td_",
                time_delay_in_Myr,
                "Myrs.h5",
            )
        else
            mkpath("catalogs/")
            path = pwd()*"/catalogs/"
            name_file = path*name_catalog
        end
    end


    # NSBH mass
    a_1 = 1.04e11
    b_1 = 2.15
    a_2 = 800
    b_2 = 0.29
    a_3 = 0.00285
    b_3 = 1.686

    # BBH mass
    lambda_peak = 0.039
    mass_min = 5.1
    mass_max = 87.0
    mu_m = 34.0
    sigma_m = 3.6
    alpha = 3.4
    delta_m = 4.8
    beta_q = 1.1

    # Redshift

    # Madau Fragos
    a_z = 0.01
    b_z = 2.6
    c_z = 0.3125
    d_z = 6.2

    if time_delay_in_Myr == nothing
        # Madau Dickinson
        println("No time delay")
        alpha_z = 2.7
        beta_z = 3.0
        z_p = 2.0

    elseif time_delay_in_Myr == 10.0
        # Madau Dickinson
        alpha_z, beta_z, z_p = (2.039337237230744, 3.1415896086773802, 1.9669095091968694)

    elseif time_delay_in_Myr == 20.0
        # Madau Dickinson
        alpha_z, beta_z, z_p = (1.9987169881669873, 3.163039104322259, 1.942109690188695)
    elseif time_delay_in_Myr == 50.0
        # Madau Dickinson
        alpha_z, beta_z, z_p = (1.9370202540742487, 3.2160803235570725, 1.8939544282143053)
    else
        error("Time delay not recognized, the ones available are nothing, 10.0, 20.0, 50.0. Nothing means the default values of the SFR, so no time delay")
    end

    # Physics constant

    clight = uc.clight_kms  #km/s
    Omega0_m = uc.Omega0_m
    Omega0_Lambda = 1 - Omega0_m
    H0 = uc.H0 # km/s 1/Mpc
    H0_yr = 2.25e-18 * 3.154e7 # 1/yr


    Random.seed!(seed)
    theta = acos.(rand(Uniform(-1, 1), n_samples))
    phi = rand(Uniform(0, 2 * pi), n_samples)
    iota = acos.(rand(Uniform(-1, 1), n_samples))
    psi = rand(Uniform(0, pi), n_samples)
    phiCoal = rand(Uniform(0, 2 * pi), n_samples)
    tcoal = rand(Uniform(0, 1), n_samples)  #fraction of a day, same as GWFAST. 
    #No need to consider overlaps now, all the events are independent and we only consider earth rotation around its axis (no needs to consider month and year)

    normalization_redshift = 0.

    if population == "BBH"

        Lambda_1 = zeros(n_samples)
        Lambda_2 = zeros(n_samples)

        chi1 = rand(Beta(1.6, 4.12), n_samples)
        chi2 = rand(Beta(1.6, 4.12), n_samples)

        # from chi to chiz requires a cosine
        # this cosine is extracted from a mixture of two distributions (truncated normal and uniform)
        # defined in arXiv:2111.03634
        xi = 0.66
        sigma_t = 1.5

        dist1 = truncated(Normal(0.0, sigma_t), -1.0, 1.0)
        dist2 = Uniform(-1.0, 1.0)
        dist_theta_spin = MixtureModel([dist1, dist2], [xi, 1 - xi])

        cos_theta_spin_1 = rand(dist_theta_spin,n_samples)
        cos_theta_spin_2 = rand(dist_theta_spin,n_samples)

        chiz1 = chi1 .* cos_theta_spin_1
        chiz2 = chi2 .* cos_theta_spin_2
        

        normalization_BHmass_BBH_first_mass = quadgk(
            M_BH -> BHmass_BBH_first_mass(
                M_BH,
                lambda_peak,
                mass_min,
                mass_max,
                mu_m,
                sigma_m,
                alpha,
                delta_m,
                S,
                smoothing,
            ),
            2.0,
            mass_max,
        )[1]
        m_1 = rejection_sampling_first_mass(
            BHmass_BBH_first_mass,
            lambda_peak,
            mass_min,
            mass_max,
            mu_m,
            sigma_m,
            alpha,
            delta_m,
            S,
            smoothing,
            normalization_BHmass_BBH_first_mass,
            [2.0, mass_max],
            n_samples,
        )
        m_2 = rejection_sampling_second_mass(
            BHmass_BBH_q,
            mass_min,
            delta_m,
            beta_q,
            S,
            smoothing,
            [1e-2, 1.0],
            n_samples,
            m_1,
            normalization_BHmass_BBH_q,
        )

        normalization_redshift = quadgk(
            z -> Redshift(
                z,
                population,
                alpha_z,
                beta_z,
                z_p,
                a_z,
                b_z,
                c_z,
                d_z,
                H0,
                Omega0_m,
                Omega0_Lambda,
                clight,
                get_dV_dz,
                friend,
            ),
            1e-3,
            10.0,
        )[1]
        z = rejection_sampling_redshift(
            Redshift,
            population,
            alpha_z,
            beta_z,
            z_p,
            a_z,
            b_z,
            c_z,
            d_z,
            H0,
            Omega0_m,
            Omega0_Lambda,
            clight,
            get_dV_dz,
            friend,
            normalization_redshift,
            [1e-3, 1e1],
            n_samples,
        )


    elseif population == "BNS"
        Lambda_1 = rand(Uniform(0, 2e3), n_samples)
        Lambda_2 = rand(Uniform(0, 2e3), n_samples)

        chiz1 = rand(Uniform(-0.05, 0.05), n_samples)
        chiz2 = rand(Uniform(-0.05, 0.05), n_samples)

        m_1 = rand(Uniform(1, 2.5), n_samples)
        m_2 = rand(Uniform(1, 2.5), n_samples)
        if m_1 < m_2
            tmp = m_1
            m_1 = m_2
            m_2 = tmp
        end
        normalization_redshift = quadgk(
            z -> Redshift(
                z,
                population,
                alpha_z,
                beta_z,
                z_p,
                a_z,
                b_z,
                c_z,
                d_z,
                H0,
                Omega0_m,
                Omega0_Lambda,
                clight,
                get_dV_dz,
                friend,
            ),
            1e-3,
            10.0,
        )[1]
        z = rejection_sampling_redshift(
            Redshift,
            population,
            alpha_z,
            beta_z,
            z_p,
            a_z,
            b_z,
            c_z,
            d_z,
            H0,
            Omega0_m,
            Omega0_Lambda,
            clight,
            get_dV_dz,
            friend,
            normalization_redshift,
            [1e-3, 1e1],
            n_samples,
        )

    elseif population == "NSBH"
        Lambda_1 = rand(Uniform(0, 2e3), n_samples)
        Lambda_2 = zeros(n_samples)

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
        normalization_redshift = quadgk(
            z -> Redshift(
                z,
                population,
                alpha_z,
                beta_z,
                z_p,
                a_z,
                b_z,
                c_z,
                d_z,
                H0,
                Omega0_m,
                Omega0_Lambda,
                clight,
                get_dV_dz,
                friend,
            ),
            1e-3,
            10.0,
        )[1]
        z = rejection_sampling_redshift(
            Redshift,
            population,
            alpha_z,
            beta_z,
            z_p,
            a_z,
            b_z,
            c_z,
            d_z,
            H0,
            Omega0_m,
            Omega0_Lambda,
            clight,
            get_dV_dz,
            friend,
            normalization_redshift,
            [1e-3, 1e1],
            n_samples,
        )
    end

    if local_rate === nothing
        # source http://arxiv.org/abs/2111.03634
        if population == "BBH"
            local_rate = 17.0
        elseif population == "BNS"
            local_rate = 105.5
        elseif population == "NSBH"
            local_rate = 45.
        end
    end
    ## number of events in an year
    # integrate the rate over redshift
    total_number_sources_yr = normalization_redshift * local_rate / 1e9 # to convert the normalization factor from Mpc^3 to Gpc^3
    # this number is not used in the code, it is just to have an estimate of the total number of sources in a year



    chirp_mass = (m_1 .* m_2) .^ (3 / 5) ./ (m_1 .+ m_2) .^ (1 / 5)
    chirp_mass_detector_frame = chirp_mass .* (1 .+ z)
    eta = (m_1 .* m_2) ./ (m_1 .+ m_2) .^ 2
    dL = get_dL(z, clight, H0, Omega0_m, Omega0_Lambda) ./ 1e3 # Gpc
    if auto_save
        date = Dates.now()
        date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
        println("Name of the catalog: ", name_file)
        h5open(name_file, "w") do file
            attributes(file)["format"] = "GWJulia"
            attributes(file)["number_events"] = nEvents
            attributes(file)["seed"] = seed
            attributes(file)["population"] = population
            attributes(file)["time_delay_in_Myrs"] = time_delay_in_Myr
            attributes(file)["SFR"] = SFR
            attributes(file)["total_number_sources_yr"] = Int(round(total_number_sources_yr, digits=0))
            attributes(file)["local_rate"] = local_rate
            attributes(file)["date"] = date_format

            write(file, "Lambda1", Lambda_1)
            write(file, "Lambda2", Lambda_2)
            write(file, "chi1", chiz1)  # chi1 = chi1z, since chi1x = chi1y = 0
            write(file, "chi2", chiz2)  # chi2 = chi2z, since chi2x = chi2y = 0
            write(file, "mc", chirp_mass_detector_frame)
            write(file, "eta", eta)
            write(file, "dL", dL)
            write(file, "z", z)
            write(file, "theta", theta)
            write(file, "phi", phi)
            write(file, "iota", iota)
            write(file, "psi", psi)
            write(file, "phiCoal", phiCoal)
            write(file, "tcoal", tcoal)
        end
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
    Lambda_1,
    Lambda_2
    

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
        ### see attributes
        println("Attributes: ", keys(attributes(file)))
        attributes_keys = keys(attributes(file))
        for i in 1:length(keys(attributes(file)))
            println(attributes_keys[i], ": ", read(attributes(file)[attributes_keys[i]]))
        end
        println("Parameters: ", keys(file))
        nEvents = read(attributes(file)["number_events"])
        population = read(attributes(file)["population"])
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



    
end
    

