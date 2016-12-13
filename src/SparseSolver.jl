#=
    A class that solves the MOC equations using a full Sparse matrix
    representation
=#
type SparseSolver

    # The TrackGenerator contains Tracks describing the problem to be solved
    _track_generator::TrackGenerator
    
    # The transport matrix containing all non-fission terms
    _T::SparseMatrixCSC{Float64, Int64}
    
    # The scalar fission matrix contains all fission terms
    _F_tilde::SparseMatrixCSC{Float64, Int64}
    
    # The fission matrix contains all fission source terms for angular fluxes
    _F::SparseMatrixCSC{Float64, Int64}
    
    # The weighting matrix converts angular fluxes to scalar fluxes
    _W::SparseMatrixCSC{Float64, Int64}
    
    # The scattering matrix contains all scattering source terms for angular
    # fluxes
    _S::SparseMatrixCSC{Float64, Int64}
    
    # A vector containing all outgoing angular fluxes for each segment
    _Ψ::Array{Float64}

    # A vector containing the current iteration's scalar fluxes
    _scalar_flux::Array{Float64}

    # A vector containing the previous iteration's scalar fluxes
    _old_scalar_flux::Array{Float64}
    
    # A matrix which computes pin fission rates from scalar fluxes
    _pin_xs::SparseMatrixCSC{Float64, Int64}
    
    # The number of angular fluxes
    _N::Int64
    
    # The number of FSRs in the geometry
    _num_FSRs::Int64
    
    # An array of all materials in the problem indexed by FSR ID
    _materials::Array{Material}
    
    # The current iteration's approximation to k-effective
    _k_eff::Float64

    # The number of polar angles in [0, π/2]
    _num_polar_2::Int64
    
    # The threshold on RMS fission rate difference for convergence
    _convergence_thresh::Float64

    # The maximum number of transport sweep iterations allowed
    _max_iters::Int64

    # Whether to lag the scattering source during transport solves
    _lag_scattering::Bool
    
    # Whether to use the "eigs" solver in Julia to compute the eigenvalues
    _use_eigs::Bool
    
    SparseSolver(track_generator::TrackGenerator, num_polar_2::Int64) = new(
                                                        track_generator,
                                                        spzeros(0,0),
                                                        spzeros(0,0),
                                                        spzeros(0,0),
                                                        spzeros(0,0),
                                                        spzeros(0,0),
                                                        zeros(0),
                                                        zeros(0),
                                                        zeros(0),
                                                        spzeros(0,0),
                                                        0,
                                                        0,
                                                        Array{Material}(0),
                                                        1.0,
                                                        num_polar_2,
                                                        1e-5,
                                                        10000,
                                                        true,
                                                        false)
end


#=
    A function that computes the sparse matrices used in the eigenvalue solve
=#
function constructMatrices(self::SparseSolver)

    println("Constructing matrices...")

    # Get the number of energy groups
    geometry = self._track_generator._geometry
    num_groups = getNumEnergyGroups(geometry)
    
    # Determine the number of angular fluxes
    num_tracks = 0
    num_segments = 0
    for a = 1:self._track_generator._num_azim_2
        num_tracks += self._track_generator._num_tracks[a]
        for i = 1:self._track_generator._num_tracks[a]
            track = self._track_generator._tracks[a][i]
            num_segments += length(track._segments)
        end
    end
    self._N = 2 * num_segments * num_groups * self._num_polar_2

    # Get the number of scalar fluxes
    num_FSRs = getNumCells(geometry)
    self._num_FSRs = num_FSRs
    M = num_FSRs * num_groups
    
    # Assign materials to FSR IDs
    self._materials = Array(Material, self._num_FSRs)
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            self._materials[id] = geometry._mesh[i, j]
        end
    end

    # Get the polar weights for the current number of polar angles
    sin_thetas, polar_weights = getPolarQuadrature(self._num_polar_2)
    
    # Allocate inverse sin thetas
    inv_sin_thetas = Array(Float64, self._num_polar_2)
    for p = 1:self._num_polar_2
        inv_sin_thetas[p] = 1.0 / sin_thetas[p]
    end

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

    # Initialize sparse matrices and arrays
    self._T = spzeros(self._N, self._N)
    self._F = spzeros(self._N, self._N)
    self._F_tilde = spzeros(M, M)
    self._W = spzeros(M, self._N)
    self._Ψ = zeros(self._N)
    self._pin_xs = spzeros(num_FSRs, M)
    self._scalar_flux = zeros(M)
    self._old_scalar_flux = zeros(M)

    # Create temporary matrices
    T_tilde = spzeros(self._N, self._N)
    S_tilde = spzeros(M, M)
    H = spzeros(self._N, M)

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

    # Fill the FSR selection matrix (H)
    i = 0
    for a = 1:num_azim_2
        for t = 1:self._track_generator._num_tracks[a]
            track = self._track_generator._tracks[a][t]
            for d = 1:2
                ns = length(track._segments)
                for s = 1:length(track._segments)
                    seg_idx = (d-1) * (ns - s + 1) + (2-d) * s
                    seg = track._segments[seg_idx]
                    sigma_t = self._materials[seg._id]._sigma_t
                    for g = 1:num_groups
                        for p = 1:self._num_polar_2
                            
                            i += 1
                            j = (seg._id-1) * num_groups + g

                            length_3D = seg._length * inv_sin_thetas[p]
                            tau = sigma_t[g] * length_3D
                            val = (1 - exp(-tau)) / (4*pi*sigma_t[g])
                            
                            H[i, j] = val
                        end
                    end
                end
            end
        end
    end

    # Pre-allocate polar arrays
    # NOTE: this looping doesn't make sense when we don't tally to the local FSR
    length_3D = Array{Float64}(self._num_polar_2)
    optical_thickness = Array{Float64}(self._num_polar_2)
    exponentials = Array{Float64}(self._num_polar_2)
    
    # Fill the weight matrix (W)
    j = 0
    for a = 1:num_azim_2
        for t = 1:self._track_generator._num_tracks[a]
            track = self._track_generator._tracks[a][t]
            for d = 1:2
                ns = length(track._segments)
                for s = 1:length(track._segments)
                    seg_idx = (d-1) * (ns - s + 1) + (2-d) * s
                    seg = track._segments[seg_idx]
                    sigma_t = self._materials[seg._id]._sigma_t
                    for g = 1:num_groups

                        # Pre-calculate exponentials
                        for p = 1:self._num_polar_2
                            length_3D[p] = seg._length * inv_sin_thetas[p]
                            optical_thickness[p] = length_3D[p] * sigma_t[g]
                            exponentials[p] = exp(-optical_thickness[p]) 
                        end
                        
                        # Compte weights
                        for p = 1:self._num_polar_2

                            # Calculate indexes and extract information
                            j += 1
                            i = (seg._id-1) * num_groups + g
                            weight = weights[(a-1) * self._num_polar_2 + p]
                            tau = optical_thickness[p]
                            exp_val = exponentials[p]
                            vol = getFSRVolume(geometry, seg._id)

                            # Set weights
                            #FIXME (when using source iteration)
                            self._W[i, j] += (weight * length_3D[p] / vol * 
                                        (1.0/(1.0 - exp_val) - 1.0/tau))
                            
                            if s != 1
                                j0 = j - self._num_polar_2 * num_groups
                                self._W[i, j0] += (weight * length_3D[p] /vol *
                                             (1.0/tau - exp_val / 
                                              (1.0 - exp_val)))
                            end
                        end
                    end
                end
            end
        end
    end

    # Fill the transport matrix (T_tilde)
    i = 0
    for a = 1:num_azim_2
        for t = 1:self._track_generator._num_tracks[a]
            track = self._track_generator._tracks[a][t]
            for d = 1:2
                ns = length(track._segments)
                for s = 1:length(track._segments)
                    seg_idx = (d-1) * (ns - s + 1) + (2-d) * s
                    seg = track._segments[seg_idx]
                    sigma_t = self._materials[seg._id]._sigma_t
                    for g = 1:num_groups
                        for p = 1:self._num_polar_2

                            i += 1
                            T_tilde[i,i] = 1.0

                            if s != 1
                                j = i - self._num_polar_2 * num_groups
                                length_3D = seg._length * inv_sin_thetas[p]
                                tau = sigma_t[g] * length_3D
                                T_tilde[i, j] = -exp(-tau)
                            end
                        end
                    end
                end
            end
        end
    end

    # Calculate initial scalar flux
    self._scalar_flux = ones(M)
    self._Ψ = self._W \ self._scalar_flux

    # Calculate MOC matrices
    self._S = H * S_tilde * self._W
    self._F = H * self._F_tilde * self._W

    # Make the transport matrix according to how scattering is treated
    if (self._lag_scattering)
        self._T = T_tilde
    else
        self._T = T_tilde - self._S
    end
end


#=
    Normalizes angular and scalar fluxes by the total fission source
=#    
function normalizeFluxes(self::SparseSolver)

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
    self._Ψ /= fission_rate
    self._scalar_flux /= fission_rate
    self._old_scalar_flux /= fission_rate
    
end


#=
    Stores the previous iteration's scalar flux
=#
function storeFSRFluxes(self::SparseSolver)
    temp = self._old_scalar_flux
    self._old_scalar_flux = self._scalar_flux
    self._scalar_flux = self._old_scalar_flux
end


#=
    Calculates the eigenvalue from fission rates using scalar fluxes
=#
function computeKeff(self::SparseSolver)
    
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
function computeResidual(self::SparseSolver)

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
    Checks to ensure that all FSRs have at least one segment traversing them
=#
function checkFSRCrossings(self)
    crossed = falses(self._num_FSRs)
    for a = 1:self._track_generator._num_azim_2
        for t = 1:self._track_generator._num_tracks[a]
            track = self._track_generator._tracks[a][t]
            for d = 1:2
                ns = length(track._segments)
                for s = 1:length(track._segments)
                    seg = track._segments[s]
                    crossed[seg._id] = true
                end
            end
        end
    end
    for i = 1:self._num_FSRs
        if !crossed[i]
            println("ERROR: No Tracks cross FSR ", i)
            exit()
        end
    end
end


#=
    Calculates the eigenvlaue of the problem represented by the TrackGenerator
=#
function computeEigenvalue(self::SparseSolver)

    # Construct the MOC matrices and ensure every FSR is traversed
    constructMatrices(self)
    checkFSRCrossings(self)

    println("[RESULT]\t\tIteration\tk-eff\t\tResidual")

    # Determine which eigenvalue method to use
    if self._use_eigs

        # Check if we can use Julia's eigs
        if self._lag_scattering
            println("ERROR: Can only use eigs if scattering is not lagged")
            exit()
        end

        # Compute eigenvalues with eigs
        sol = eigs(self._F, self._T, which=:LM, nev=1)
        self._k_eff = convert(Float64,real(sol[1][1]))
        comp_fluxes = sol[2]
        iter = sol[4]
        nmult = sol[5]
        resid = sol[6]
        println("[RESULT]", "\t\t", iter, "\t\t", round(self._k_eff,5),
                "\t\t", sum(resid))

    else

        # Use "power method"
        for iter = 1:self._max_iters

            # Normalize the fluxes for the next iteration
            normalizeFluxes(self)

            # Perform a "transport sweep"
            if (self._lag_scattering)
                total_source = ((self._F * self._Ψ) / self._k_eff
                                + self._S * self._Ψ)
                self._Ψ = self._T \ total_source
            else
                fission_source = (self._F * self._Ψ) / self._k_eff
                self._Ψ = self._T \ fission_source
            end

            # Calculate the new scalar flux and eigenvalue
            self._scalar_flux = self._W * self._Ψ
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
end
