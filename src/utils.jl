module UtilsAndConstants    # this is the name of the module

using ForwardDiff # these are the Julia packages that are used in this module
using LinearAlgebra
using LaTeXStrings
using Base.Threads
using Unitful
using UnitfulAstro 


export GMsun_over_c3, GMsun_over_c2, uGpc, GMsun_over_c2_Gpc, REarth_km, clight_kms, clightGpc, Lamt_delLam_from_Lam12,
         _ra_dec_from_theta_phi_rad, _theta_phi_from_ra_dec_rad, CovMatrix, Errors, SkyArea, _orientationBigCircle, CovMatrix_Lamt_delLam,
         masses_from_chirp_eta, masses_z_to_chirp_eta


##############################################################################
GMsun_over_c3 = 4.925491025543575903411922162094833998e-6 # seconds

GMsun_over_c2 = 1.476625061404649406193430731479084713e3 # meters

uGpc = 3.085677581491367278913937957796471611e25 # meters

GMsun_over_c2_Gpc = GMsun_over_c2 / uGpc # Gpc

REarth_km = 6371.00 # km

clight_kms = 2.99792458e5 # km/s

clightGpc = clight_kms / 3.0856778570831e+22

"""Solar mass"""
MSUN = 1.988409902147041637325262574352366540e30  # kg

"""Geometrized nominal solar mass, m"""
MRSUN = GMsun_over_c2

Omega0_m = 0.3153

Omega0_Lambda = 1 - Omega0_m

H0 = 67.66u"km/s/Mpc" # km/s 1/Mpc


##############################################################################

"""
Compute the masses of the two objects, in the source frame, from the chirp mass, in the detector frame, the symmetric mass ratio and the redshift.
The chirp mass is defined as M_c = (m1*m2)^(3/5)/(m1+m2)^(1/5) * (1+z) and the symmetric mass ratio is defined as eta = m1*m2/(m1+m2)^2.
#### Input arguments:
-  float or vector mc chirp mass
-  float or vector eta symmetric mass ratio
-  float or vector z redshift
#### Outputs:
-  (float, float) or (vector, vector) m1 and m2 Mass of object 1 and object 2
#### Example:
```julia
    mc = 1.0
    eta = 0.25
    z = 0.1
    m1, m2 = masses_from_chirp_eta(mc, eta, z)
    println("m1: ", m1)
    println("m2: ", m2)
``` 
"""
function masses_from_chirp_eta(
    mc::Union{Float64, Vector{Float64}},
    eta::Union{Float64, Vector{Float64}},
    z::Union{Float64, Vector{Float64}}
    )
    # check if mc, eta and z are the same size
    if length(mc) != length(eta) || length(mc) != length(z)
        throw(ArgumentError("mc, eta and z must be the same size"))
    end

    mc_source = @. mc/(1+z)
    M = @. mc_source/eta^(3/5)
    Seta = @. sqrt(1 - 4*eta)
    m1 = @. 0.5 * (1.0 + Seta) * M
    m2 = @. 0.5 * (1.0 - Seta) * M
    return m1, m2
end


"""
Compute the chirp mass and symmetric mass ratio from the masses of the two objects, in the source frame, and the redshift.
The chirp mass is defined as M_c = (m1*m2)^(3/5)/(m1+m2)^(1/5) * (1+z) and the symmetric mass ratio is defined as eta = m1*m2/(m1+m2)^2.
#### Input arguments:
-  float or vector m1 mass of object 1
-  float or vector m2 mass of object 2
-  float or vector z redshift
#### Outputs:
-  (float, float) or (vector, vector) mc and eta chirp mass and symmetric mass ratio
#### Example:
```julia
    m1 = 1.0
    m2 = 0.5
    z = 0.1
    mc, eta = masses_z_to_chirp_eta(m1, m2, z)
    println("mc: ", mc)
    println("eta: ", eta)
```
"""
function masses_z_to_chirp_eta(m1::Union{Float64, Vector{Float64}}, m2::Union{Float64, Vector{Float64}}, z::Union{Float64, Vector{Float64}})
    mc = @. ((m1*m2)^(3/5))/((m1+m2)^(1/5))*(1+z) # chirp mass at detector
    eta = @. m1*m2/(m1 + m2)^2
    return mc, eta
end

r"""
Compute the dimensionless tidal deformability combinations L"\tilde{Lambda}" and L"\delta\tilde{Lambda}", defined in `arXiv:1402.5156 <https://arxiv.org/abs/1402.5156>`_ eq. (5) and (6), as a function of the dimensionless tidal deformabilities of the two objects and the symmetric mass ratio.
Function from GWFAST.

#### Input argument:
-  float Lambda1 Tidal deformability of object 1
-  float Lambda2 Tidal deformability of object 2
-  float eta The symmetric mass ratio of the object.

#### Outputs:
- (float, float) L"\tilde{Lambda}" and L"\delta\tilde{Lambda}", 

"""
function Lamt_delLam_from_Lam12(eta, Lambda1, Lambda2)

    eta2 = eta*eta

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
        
    Lamt = (8. /13.)*((1. + 7. *eta - 31. *eta2)*(Lambda1 + Lambda2) + Seta*(1. + 9. *eta - 11. *eta2)*(Lambda1 - Lambda2))
    
    delLam = 0.5*(Seta*(1. - 13272. /1319. *eta + 8944. /1319. *eta2)*(Lambda1 + Lambda2) + (1. - 15910. /1319. *eta + 32850. /1319. *eta2 + 3380. /1319. *eta2*eta)*(Lambda1 - Lambda2))
    
    return Lamt, delLam

end

r"""
Compute the covariance matrix of L"\tilde{Lambda}" and L"\delta\tilde{Lambda}" from the Fisher matrix of the tidal deformabilities.

#### Input arguments:
-  Fisher Fisher matrix.
-  Lambda1 Tidal deformability of object 1
-  Lambda2 Tidal deformability of object 2
-  eta The symmetric mass ratio of the object.

#### Outputs:
-  Covariance matrix of L"\tilde{Lambda}" and L"\delta\tilde{Lambda}".

"""
function CovMatrix_Lamt_delLam(Cov, eta, Lambda1, Lambda2)

    #Cov_reduced = CovMatrix(Fisher, debug=false)[[2, 12, 13], [2, 12, 13]]
    Cov_reduced = Cov[[2, 12, 13], [2, 12, 13]]
    # jacobian
    J1 = ForwardDiff.gradient(x->Lamt_delLam_from_Lam12(x...)[1], [eta, Lambda1, Lambda2])
    J2 = ForwardDiff.gradient(x->Lamt_delLam_from_Lam12(x...)[2], [eta, Lambda1, Lambda2])
    J = [J1 J2]'
    # covariance matrix of Lambda tilde and delta Lambda
    Cov_Lamt_delLam = J * Cov_reduced * J'
    return Cov_Lamt_delLam
end



"""
Function to obtain RA and DEC from theta and phi in radians. Result is in radians.

#### Input arguments:
-  float theta Polar angle in radians.
-  float phi Azimuthal angle in radians.

#### Outputs:
-  (float, float) RA and DEC in radians.

#### Example:
```julia
ra, dec = _ra_dec_from_theta_phi_rad(0.1, 0.2)
```

"""
function _ra_dec_from_theta_phi_rad(theta::Union{Float64, ForwardDiff.Dual}, phi::Union{Float64, ForwardDiff.Dual})

    ra = phi
    dec = 0.5 * pi - theta
    return ra, dec
end

function _ra_dec_from_theta_phi_rad(theta::AbstractArray, phi::AbstractArray)
    ra = phi
    dec = [0.5 * pi - theta[i] for i in eachindex(theta)]
    return ra, dec
end

function _theta_phi_from_ra_dec_rad(ra::Union{Float64, ForwardDiff.Dual}, dec::Union{Float64, ForwardDiff.Dual})
    theta = 0.5 * pi - dec
    phi = ra
    return theta, phi
end

function _theta_phi_from_ra_dec_rad(ra::AbstractArray, dec::AbstractArray)
    theta = [0.5 * pi - dec[i] for i in eachindex(dec)]
    phi = ra
    return theta, phi
end

"""
Function to compute the covariance matrix from the Fisher matrix. 
The function uses the Cholesky decomposition to invert the Fisher matrix. Some matrices can be difficult to invert, so if the inversion fails, it tries to normalize the Fisher matrix and invert it again. 
If the inversion fails again, it tries to invert the Fisher matrix with 128 bit precision. If also this fails, it returns a zero matrix.
If the Fisher matrix is made of zeros, it returns a zero matrix.
#### Input arguments:
-  Fisher matrix: square matrix of size NxN.

#### Optional arguments:
-  debug: boolean to print debug information. Default is true.
-  threshold: float with the threshold to check if the inversion is correctly done. Default is 5e-2.
-  called_by_3D_function: boolean to check if the function is called by the CovMatrix() function that deals with arrays of covariances. It is used to keep track of the number of failed inversions. Default is false.
-  force_high_precision: boolean to force the inversion with 128 bit precision. Default is false.

#### Outputs:
-  Covariance matrix: square matrix of size NxN.

#### Example:
```julia
    covMatrix = CovMatrix([1 0; 0 1])
```

"""
function CovMatrix(Fisher::Matrix{Float64}; debug = true, threshold = 5e-2, called_by_3D_function = false, force_high_precision = false)
    covMatrix = zeros(size(Fisher))  # Initialize covMatrix as a zero matrix of the same size as Fisher
    
    idx = 0 # index to check if the inversion is correctly done
    # idx = 0 if the inversion is correctly done
    # idx = 1 if the Fisher matrix is not positive definite (try to normalize it)
    # idx = 2 if the Fisher matrix is not positive definite even after normalization
    # idx = 3 if the inversion failed with 128 bit precision

    # we are interested in the case idx = 3, return zero matrix and idx/3 is the Fisher matrix is not positive definite even
    # after normalization
    if Fisher[1,1]==0.
        if called_by_3D_function
            return zeros(size(Fisher)), idx
        else
            return zeros(size(Fisher))
        end
    end

    if force_high_precision
        covMatrix=setprecision(128) do # idea to go 128 bit precision from GWFAST

            try 
                Fisher_BF = BigFloat.(Fisher) # convert Fisher matrix to BigFloat (BF)
                covMatrix = inv(cholesky(Fisher_BF)) 
                return Float64.(covMatrix)
            catch
                idx = 3
                if debug == true
                    println("Inversion failed with 128 bit precision")
                end
                return zeros(size(Fisher))
            end


        end

        if called_by_3D_function
            return covMatrix, idx
        else
            return covMatrix
        end
    end

    # try to invert the Fisher matrix, first with Float64 precision
    try
        covMatrix = inv(cholesky(Fisher))
    catch
        idx = 1
        covMatrix = zeros(size(Fisher))
        if debug == true
            println("Fisher matrix is not positive definite. Try normalize it.")
        end
    end

    inversion_error = maximum(abs.(Fisher * covMatrix - I))

    if idx == 1 || inversion_error > threshold*1e-2
        # normalize the fisher, technique by GWFISH 
        diagonal = diag(Fisher)
        Fisher_ = deepcopy(Fisher) ./ sqrt.(diagonal * diagonal')
        try
            covMatrix_ = inv(cholesky(Fisher_))
            covMatrix_ = covMatrix_ ./ sqrt.(diagonal * diagonal')
            inversion_error_ = maximum(abs.(Fisher * covMatrix_ - I))

            if debug == true
                println("Inversion error before normalization: ", inversion_error)
                println("Inversion error after normalization: ", inversion_error_)
            end

            if inversion_error_ < inversion_error
                
                covMatrix = covMatrix_
                inversion_error = inversion_error_            
            end

        catch
            idx = 2
            covMatrix = zeros(size(Fisher))
            if debug == true
                println("Fisher matrix is not positive definite even after normalization")
            end
        end

    end

    if idx == 2 || inversion_error > threshold*1e-2
        # try to invert the Fisher matrix with 128 bit precision
        covMatrix = setprecision(128) do # idea to go 128 bit precision from GWFAST
            if debug == true
                println("Try to invert with 128 bit precision")
            end

            try 
                Fisher_BF = BigFloat.(Fisher) # convert Fisher matrix to BigFloat (BF)
                covMatrix = inv(cholesky(Fisher_BF))
                return covMatrix
            catch
                idx = 3
                # covMatrix = zeros(size(Fisher))
                if debug == true
                    println("Inversion failed with 128 bit precision")
                end
                return zeros(size(Fisher))
            end


        end

    end

    # check if the inversion is correctly done
    inversion = isapprox(Fisher * covMatrix, I, atol=threshold) 
    if inversion == false # failed inversion, return zero matrix
        idx = 3

        if debug == true
            println("Inversion failed")
            inversion_error = maximum(abs.(Fisher * covMatrix - I))
            println("Inversion error: ", Float64(inversion_error))
        end

        if called_by_3D_function
            return zeros(size(Fisher)), idx
        else
            return zeros(size(Fisher))
        end
    else
        if debug == true
            println("Inversion successful")
            inversion_error = maximum(abs.(Fisher * covMatrix - I))
            println("Inversion error: ", Float64(inversion_error))
        end
        idx = 0 # successful inversion
    end



    # return covMatrix in 64 bit precision
    covMatrix = Float64.(covMatrix)

    if debug == true
        println()
    end

    if called_by_3D_function
        return covMatrix, idx
    else
        return covMatrix
    end

end

"""
Same as CovMatrix(Fisher::Matrix{Float64}) but for a 3D array of Fisher matrices. So array of size MxNxN. where M is the number of Fisher matrices.

"""
function CovMatrix(Fisher::Array{Float64, 3}; threshold = 5e-2, force_high_precision = false, debug = false)
    covMatrix = zeros(size(Fisher))  # Initialize covMatrix as a zero matrix of the same size as Fisher

    n = size(Fisher)[1]
    failed_inv = 0

    @threads for i in 1:n
        idx=0
        covMatrix[i,:,:], idx= CovMatrix(Fisher[i,:,:], debug = debug, threshold = threshold, called_by_3D_function = true, force_high_precision = force_high_precision)
        failed_inv += Int(idx/3)
    end
    println("Failed inversions: ", failed_inv)
    return covMatrix
end
"""
Calculate the errors from the covariance matrix. The errors are the square root of the diagonal elements of the covariance matrix.

#### Input arguments:
-  covMatrix: square matrix of size NxN.

#### Outputs:
-  errors: array of size N.

#### Example:
```julia
    errors = Errors([1 0; 0 4]) # returns [1.0, 2.0]
"""
function Errors(covMatrix::Matrix{Float64})
    return sqrt.(diag(covMatrix))
end

"""
Same as Errors(covMatrix::Matrix{Float64}) but for a 3D array of covariance matrices.
"""
function Errors(covMatrix::Array{Float64, 3})
    
    nCov = size(covMatrix)[1]
    nPar = size(covMatrix)[2]
    errors = zeros(nCov, nPar)

    @threads for ii in 1:nCov
        errors[ii,:] = sqrt.(diag(covMatrix[ii,:,:]))
    end

    return errors
end

"""
Calculate the sky area from the covariance matrix and the catalog of sources. The sky area is the area of the error ellipse in the sky, in square degrees.

#### Input arguments:
-  covMatrix: square matrix of size NxN.
-  thetaCatalog: float with the polar angle of the source in the catalog.

#### Optional arguments:
-  percent_level: float with the confidence level. Default is 90%.

#### Outputs:
-  float with the sky area in square degrees.

#### Example:
```julia
    covMatrix = Matrix{Float64}(I, 11,11)
    sky_area = SkyArea(covMatrix, 0.1) 
"""
function SkyArea(covMatrix::Matrix{Float64}, thetaCatalog; percent_level = 90)
    # take phi and theta from covMatrix and multiply them

    sigmaTheta_sq  = covMatrix[6,6]
    sigmaPhi_sq = covMatrix[7,7]
    sigmaThetaPhi = covMatrix[6,7]
    
    # From Barak, Cutler, PRD 69, 082005 (2004), gr-qc/0310125
    # results in grad
    return - 2*pi*sqrt(sigmaTheta_sq*sigmaPhi_sq - sigmaThetaPhi^2)*abs(sin(thetaCatalog))*log(1-percent_level/100)*(180/pi)^2
end

function SkyArea(covMatrix::Array{Float64,3}, thetaCatalog::Array; percent_level = 90)
    
    nCov = size(covMatrix)[1]
    skyArea = zeros(nCov)

    @threads for ii in 1:nCov
        skyArea[ii] = SkyArea(covMatrix[ii,:,:], thetaCatalog[ii], percent_level=percent_level)
    end

    return skyArea
end

""" 
Calculate the optimal angle for the orientation of the detector to maximize the signal from CBC gravitational wave. 
The function can be used in the definition of 2 detectors network.
It outputs the orientation of the second detectors such that it is at 45 degrees w.r.t. the first detector.
These reasoning works for all orientation conventions and in the GWJulia convention, it is the angle w.r.t. the local East.
If you want to use it for a different orientation, you need to add the angle of the first detector to the output of this function.

#### Input arguments:
- DetectorStructure det1: DetectorStructure object of the first detector.
- DetectorStructure det2: DetectorStructure object of the second detector.

#### Outputs:
-  float orientation_CBC Optimal angle for the optimal orientation of the second detector in radians

#### Example:
```julia
    orientation_CBC = _orientationBigCircle(0.1, 0.2, 0.3, 0.4)
 ```
""" 
function _orientationBigCircle(det1, det2)

    lat1 = det1.latitude_rad
    long1 = det1.longitude_rad
    orientation1 = det1.orientation_rad
    lat2 = det2.latitude_rad
    long2 = det2.longitude_rad


    equatorialRadiusEarth = 6378.137 #from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html 
    polarRadiusEarth = 6356.752 

    eFactor = 1 - (equatorialRadiusEarth^2 - polarRadiusEarth^2) / equatorialRadiusEarth^2  # 1-e^2
    eS = (equatorialRadiusEarth^2 - polarRadiusEarth^2) / equatorialRadiusEarth^2 # e^2

    Lambda1 = eFactor * tan(lat2)/tan(lat1) + eS*sqrt( (1 + eFactor * tan(lat2)^2) / (1 + eFactor * tan(lat1)^2))
    Lambda2 = eFactor * tan(lat1)/tan(lat2) + eS*sqrt( (1 + eFactor * tan(lat1)^2) / (1 + eFactor * tan(lat2)^2))

    alpha1 = atan( sin(long2 - long1) , ((Lambda1 - cos(long2 - long1))*sin(lat1)))
    alpha2 = atan( sin(long1 - long2) , ((Lambda2 - cos(long1 - long2))*sin(lat2))) + pi

    angle_ = orientation1 + alpha1 - alpha2 + pi/4
    return  mod(angle_ + pi, 2pi) - pi

end


end