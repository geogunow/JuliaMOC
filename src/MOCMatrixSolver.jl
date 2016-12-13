#=
    An MOCMatrix is an object that implicitly forms the MOC matrix for use in
    MOC eigenvalue calculations. A multiplication with a vector is equivalent
    to doing a transport sweep with the vector being the neutron sources.
=#
type MOCMatrix

    # The TrackGenerator contains Tracks describing the problem to be solved
    _track_generator::TrackGenerator
    
    # The number of scalar fluxes
    _M::Int64

    # The number of energy groups
    _num_groups::Int64
    
    # The number of polar angles in [0, π/2]
    _num_polar_2::Int64
    
    # An array of all materials in the problem indexed by FSR ID
    _materials::Array{Material}
    
    # The inverse of the sine of each polar angle (1.0 / sin(theta))
    _inv_sin_thetas::Array{Float64}
    
    # The weights for each direction (each azimuthal, polar combination)
    _quad_weights::Array{Float64}

    # The constructor for an MOC matrix
    MOCMatrix(track_generator::TrackGenerator, num_polar_2::Int64) = new(
                                                        track_generator,
                                                        0,
                                                        0,
                                                        num_polar_2,
                                                        Array{Material}(0),
                                                        zeros(0),
                                                        zeros(0))
    
end


#=
    Function to perform a matrix-vector multiplication with an MOCMatrix.
    Note that this is the same as doing an MOC transport sweep. The MOC matrix
    is equivalent to:
        W * inv(T_tilde) * H
    where
        W is the weights matrix
        T_tilde is the transport matrix
        H is the FSR selection matrix
    The vector which the matrix multiplies should contain the neutron sources
=#
function Base.:*(A::MOCMatrix, sources::Array{Float64})
    
    # Allocate a scalar flux array as well as a temporary angular flux array
    scalar_flux = zeros(A._M)
    psi = zeros(A._num_groups * A._num_polar_2)

    # Do the transport sweep. This is the same thing as multiplying by:
    #       W * inv(T_tilde) * H
    # in the SparseSolver
    for a = 1:A._track_generator._num_azim_2
        for t = 1:A._track_generator._num_tracks[a]
            track = A._track_generator._tracks[a][t]
            for d = 1:2
                
                psi *= 0
                    
                ns = length(track._segments)
                for s = 1:ns   

                    seg_idx = (d-1) * (ns - s + 1) + (2-d) * s
                    seg = track._segments[seg_idx]
                    sigma_t = A._materials[seg._id]._sigma_t
                    vol = getFSRVolume(A._track_generator._geometry, seg._id)
                    
                    for g = 1:A._num_groups
                        for p = 1:A._num_polar_2

                            # Compute indexes
                            idx = (g-1) * A._num_polar_2 + p
                            scalar_id = (seg._id-1) * A._num_groups + g
                            weight_index = (a-1) * A._num_polar_2 + p 

                            # Determine the 3D segment length and optical
                            # thickness
                            length_3D = seg._length * A._inv_sin_thetas[p]
                            tau = sigma_t[g] * length_3D
                            
                            # Determine the exponential
                            exponential = exp(-tau)

                            # Calculate matrix components on-the-fly
                            
                            # Calculate the FSR selection matrix values (H)
                            h_value = (1 - exponential) / (4*pi*sigma_t[g])

                            # Calculate weights for incoming and outgoing
                            # angular fluxes (W Matrix Values)
                            weight_out = (A._quad_weights[weight_index] / vol
                                      * length_3D * (1.0/(1.0 - exponential)
                                                     - 1.0/tau))
                            weight_in = (A._quad_weights[weight_index] / vol
                                      * length_3D * (exponential / 
                                                     (1.0 - exponential)
                                                     - 1.0/tau))

                            # Note: an optimization can be made such that
                            # we compute weights as just
                            #       -A._quad_weights[weight_index] / vol *
                            #       length_3D / tau
                            # but add the reduced source - this works because
                            # of cancelling terms

                            # Save the incoming angular flux
                            psi_old = psi[idx]
                           
                            # Determine the outgoing angular flux
                            # Implicitly we are using matrix T whose diagonal
                            # is ones and the off diagonal is the exponential
                            psi[idx] = (h_value * sources[scalar_id] 
                                        + psi[idx] * exponential)

                            # Add the incoming and outgoing angular fluxes
                            # with their weights to the current FSR
                            # This is the same as (W * Ψ)
                            scalar_flux[scalar_id] += weight_out * psi[idx]
                            scalar_flux[scalar_id] -= weight_in * psi_old
                        end
                    end
                end
            end
        end
    end
    return scalar_flux
end


#=
    A class that solves the MOC equations using a full Sparse matrix
    representation
=#
type MOCMatrixSolver

    # The TrackGenerator contains Tracks describing the problem to be solved
    _track_generator::TrackGenerator
    
    # The scalar fission matrix contains all fission terms
    _F_tilde::SparseMatrixCSC{Float64, Int64}
    
    # The scattering matrix contains all scattering source terms for angular
    # fluxes
    _S_tilde::SparseMatrixCSC{Float64, Int64}
    
    # A vector containing the current iteration's scalar fluxes
    _scalar_flux::Array{Float64}

    # A vector containing the previous iteration's scalar fluxes
    _old_scalar_flux::Array{Float64}

    # A matrix which computes pin fission rates from scalar fluxes
    _pin_xs::SparseMatrixCSC{Float64, Int64}

    # A matrix representing the MOC transport sweep
    _transport_sweep::MOCMatrix
    
    # The number of FSRs in the geometry
    _num_FSRs::Int64
    
    # The current iteration's approximation to k-effective
    _k_eff::Float64

    # The number of polar angles in [0, π/2]
    _num_polar_2::Int64
    
    # The threshold on RMS fission rate difference for convergence
    _convergence_thresh::Float64

    # The maximum number of transport sweep iterations allowed
    _max_iters::Int64

    MOCMatrixSolver(track_generator::TrackGenerator, num_polar_2::Int64) = new(
                                                    track_generator,
                                                    spzeros(0,0),
                                                    spzeros(0,0),
                                                    zeros(0),
                                                    zeros(0),
                                                    spzeros(0,0),
                                                    MOCMatrix(track_generator,
                                                              num_polar_2),
                                                    0,
                                                    1.0,
                                                    num_polar_2,
                                                    1e-5,
                                                    10000)
end


#=
    A function that computes the sparse matrices used in the eigenvalue solve
=#
function constructMatrices(self::MOCMatrixSolver)

    println("Constructing matrices...")

    # Get the number of energy groups
    geometry = self._track_generator._geometry
    num_groups = getNumEnergyGroups(geometry)
    self._transport_sweep._num_groups = num_groups
    
    # Get the number of scalar fluxes
    num_FSRs = getNumCells(geometry)
    self._num_FSRs = num_FSRs
    M = num_FSRs * num_groups
    self._transport_sweep._M = M
    
    # Assign materials to FSR IDs
    self._transport_sweep._materials = Array(Material, self._num_FSRs)
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            self._transport_sweep._materials[id] = geometry._mesh[i, j]
        end
    end
    
    # Get the polar weights for the current number of polar angles
    sin_thetas, polar_weights = getPolarQuadrature(self._num_polar_2)
    
    # Allocate inverse sin thetas
    inv_sin_thetas = Array(Float64, self._num_polar_2)
    for p = 1:self._num_polar_2
        inv_sin_thetas[p] = 1.0 / sin_thetas[p]
    end
    self._transport_sweep._inv_sin_thetas = inv_sin_thetas

    # Compute the total weights
    num_azim_2 = self._track_generator._num_azim_2
    weights = Array(Float64, num_azim_2 * self._num_polar_2)
    for a = 1:num_azim_2
        for p = 1:self._num_polar_2
            weight = (4.0 * pi * self._track_generator._azim_weights[a] *
                      self._track_generator._azim_spacings[a] * polar_weights[p])
            weight *= sin_thetas[p]
            weights[(a-1) * self._num_polar_2 + p] = weight
        end
    end 
    self._transport_sweep._quad_weights = weights

    # Initialize sparse matrices and arrays
    self._F_tilde = spzeros(M, M)
    self._pin_xs = spzeros(num_FSRs, M)
    self._scalar_flux = zeros(M)
    self._old_scalar_flux = zeros(M)

    # Create temporary matrices
    S_tilde = spzeros(M, M)

    # Fill the reaction rate matrices (F_tilde, S_tilde, and pin xs)
    for x = 1:geometry._num_x
        for y = 1:geometry._num_y
            r = getFSRId(geometry, x, y)
            mat = geometry._mesh[x,y]
            for g = 1:num_groups
                i = (r-1) * num_groups + g
                for gp = 1:num_groups
                    j = (r-1) * num_groups + gp
                    self._F_tilde[i, j] = mat._chi[g] * mat._nu_sigma_f[gp]
                    S_tilde[i, j] = mat._sigma_s[gp, g]
                end
                self._pin_xs[r, i] = mat._nu_sigma_f[g]
            end
        end
    end

    # Pre-allocate polar arrays
    # NOTE: this looping doesn't make sense when we don't tally to the local FSR
    length_3D = Array{Float64}(self._num_polar_2)
    optical_thickness = Array{Float64}(self._num_polar_2)
    exponentials = Array{Float64}(self._num_polar_2)
    
    # Calculate initial scalar flux
    self._scalar_flux = ones(M)

    # Calculate MOC matrices
    self._S_tilde = S_tilde
end


#=
    Normalizes angular and scalar fluxes by the total fission source
=#    
function normalizeFluxes(self::MOCMatrixSolver)

    # Calculate the total fission rate
    M, N = size(self._F_tilde)
    num_groups = convert(Int64, M/self._num_FSRs)
    volumes = spzeros(M, M)
    for r = 1:self._num_FSRs
        for g = 1:num_groups
            idx = (r-1)*num_groups+g
            volumes[idx,idx] =
                getFSRVolume(self._track_generator._geometry, r) 
        end
    end
    fission_rate = sum(volumes * self._F_tilde * self._scalar_flux)
    
    # Normalize fluxes so that the sum of the fission rates is 1.0
    self._scalar_flux /= fission_rate
    self._old_scalar_flux /= fission_rate
    
end


#=
    Stores the previous iteration's scalar flux
=#
function storeFSRFluxes(self::MOCMatrixSolver)
    temp = self._old_scalar_flux
    self._old_scalar_flux = self._scalar_flux
    self._scalar_flux = self._old_scalar_flux
end


#=
    Calculates the eigenvalue from fission rates using scalar fluxes
=#
function computeKeff(self::MOCMatrixSolver)
    
    # Calculate the total fission rate
    M, N = size(self._F_tilde)
    num_groups = convert(Int64, M/self._num_FSRs)
    volumes = spzeros(M, M)
    for r = 1:self._num_FSRs
        for g = 1:num_groups
            idx = (r-1)*num_groups+g
            volumes[idx,idx] =
                getFSRVolume(self._track_generator._geometry, r) 
        end
    end
    fission_rate = sum(volumes * self._F_tilde * self._scalar_flux)

    # Compute new k-effective
    self._k_eff *= fission_rate

end


#=
    Calculates the residual using the RMS difference in fission rates
=#
function computeResidual(self::MOCMatrixSolver)

    # Calcualte fission rates from scalar fluxes
    pin_fission_rate = self._pin_xs * self._scalar_flux
    old_pin_fission_rate = self._pin_xs * self._old_scalar_flux
    for i = 1:length(pin_fission_rate)
        if pin_fission_rate[i] == 0.0
            pin_fission_rate[i] = 1e-12
        end
        if old_pin_fission_rate[i] == 0.0
            old_pin_fission_rate[i] = 1e-12
        end
    end
    
    # Calculate the residual
    res = sqrt(sum((pin_fission_rate ./ old_pin_fission_rate - 1).^2))
    res /= sqrt(self._num_FSRs)
    
    return res
end


#=
    Calculates the eigenvlaue of the problem represented by the TrackGenerator
=#
function computeEigenvalue(self::MOCMatrixSolver)

    # Construct the MOC matrices and ensure every FSR is traversed
    constructMatrices(self)
    checkFSRCrossings(self)
    
    println("[RESULT]\t\tIteration\tk-eff\t\tResidual")

    # Use "power method"
    for iter = 1:self._max_iters
 
        # Normalize the fluxes for the next iteration
        normalizeFluxes(self)
        
        # Calculate the source vector
        sources = self._F_tilde * self._scalar_flux
        sources /= self._k_eff
        sources += self._S_tilde * self._scalar_flux

        # Do a "transport sweep"
        self._scalar_flux = self._transport_sweep * sources

        # Calculate the new scalar flux and eigenvalue
        computeKeff(self)
        
        # Print the current eigenvalue
        res = computeResidual(self)
        storeFSRFluxes(self)
        println("[RESULT]", "\t\t", iter, "\t\t", round(self._k_eff,5),
                "\t\t", res)
        if abs(self._k_eff) > 100.0
            println("[ERROR]\t\tFailed to converge eigenvalue")
            break
        end
        
        # Check for convergence
        if (iter > 2 && res < self._convergence_thresh)
            break
        end
    end
end
