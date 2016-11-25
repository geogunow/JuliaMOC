#=
    A class that solves the MOC equations in the traditional way
=#
type StandardSolver

    # The TrackGenerator contains Tracks describing the problem to be solved
    _track_generator::TrackGenerator

    # An array of neutron angular fluxes for each Track
    _angular_fluxes::Array{Float64}

    # An array of the scalar neutron fluxes for each FSR
    _scalar_fluxes::Array{Float64}

    # An array of the scalar neutron form the previous iteration
    _old_scalar_fluxes::Array{Float64}

    # An array of the total neutron sources for each FSR divided by the total
    # cross-section and volume
    _reduced_sources::Array{Float64}

    # The current iteration's approximation to k-effective
    _k_eff::Float64

    # An array of all materials in the problem indexed by FSR ID
    _materials::Array{Material}

    # The weights for each direction (each azimuthal, polar combination)
    _weights::Array{Float64}

    # The inverse of the sine of each polar angle (1.0 / sin(theta))
    _inv_sin_thetas::Array{Float64}

    # Indicates whether arrays have been allocated for the current solver
    _memory_allocated::Bool

    # The number of FSRs in the geometry
    _num_FSRs::Int64

    # The number of Tracks in the geometry
    _num_tracks::Int64

    # The number of polar angles in [0, π/2]
    _num_polar_2::Int64

    # The number of energy groups
    _num_groups::Int64

    # The threshold on RMS fission rate difference for convergence
    _convergence_thresh::Float64

    # The maximum number of transport sweep iterations allowed
    _max_iters::Int64

    # The constructor for the standard solver
    StandardSolver(track_generator::TrackGenerator, num_polar_2::Int64) =
                                                new(track_generator,
                                                Array{Float64}(0),
                                                Array{Float64}(0),
                                                Array{Float64}(0),
                                                Array{Float64}(0),
                                                1.0,
                                                Array{Material}(0),
                                                Array{Float64}(0),
                                                Array{Float64}(0),
                                                false,
                                                0,
                                                0,
                                                num_polar_2,
                                                0,
                                                1e-5,
                                                10000)

end


#=
    Allocates memory for angular and scalar fluxes
=#
function allocateMemory(self::StandardSolver)

    # Determine the number of FSRs
    geometry = self._track_generator._geometry
    self._num_FSRs = geometry._num_x * geometry._num_y
    self._num_groups = getNumEnergyGroups(geometry)

    # Allocate angular fluxes
    geometry = self._track_generator._geometry
    self._num_tracks = getNumTotalTracks(self._track_generator)
    self._angular_fluxes = Array(Float64, 2 * self._num_tracks *
                                 self._num_groups * self._num_polar_2)

    # Allocate scalar fluxes
    self._num_FSRs = getNumCells(geometry)
    self._scalar_fluxes = Array(Float64, self._num_FSRs * self._num_groups)
    self._old_scalar_fluxes = Array(Float64, self._num_FSRs * self._num_groups)

    # Allocate sources
    self._reduced_sources = Array(Float64, self._num_FSRs * self._num_groups)

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

    # Allocate weights
    num_azim_2 = self._track_generator._num_azim_2
    self._weights = Array(Float64, num_azim_2 * self._num_polar_2)
    for a = 1:num_azim_2
        for p = 1:self._num_polar_2
            weight = (4.0 * pi * self._track_generator._azim_weights[a] *
                      self._track_generator._azim_spacings[a] * polar_weights[p])
            weight *= sin_thetas[p]
            self._weights[(a-1) * self._num_polar_2 + p] = weight
        end
    end
    
    # Allocate inverse sin thetas
    self._inv_sin_thetas = Array(Float64, self._num_polar_2)
    for p = 1:self._num_polar_2
        self._inv_sin_thetas[p] = 1.0 / sin_thetas[p]
    end
end


#=
    A simple function that sets all scalar fluxes to a set value
=#
function flattenFSRFluxes(self::StandardSolver, value::Float64)
    for i = 1:self._num_FSRs * self._num_groups
        self._scalar_fluxes[i] = value
    end
end


#=
    A simple function that sets all angular fluxes to zero
=#
function zeroAngularFluxes(self::StandardSolver)
    num_fluxes = 2 * self._num_tracks * self._num_groups * self._num_polar_2
    for i = 1:num_fluxes
        self._angular_fluxes[i] = 0
    end
end


#=
    Normalizes scalar fluxes by the total fission source
=#
function normalizeFluxes(self::StandardSolver)

    # Compute the total fission source
    geometry = self._track_generator._geometry
    fission_source = 0.0
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            mat = geometry._mesh[i, j]
            vol = getFSRVolume(geometry, id)
            local_fission_source = 0.0
            for g = 1:self._num_groups
                index = (id-1) * self._num_groups + g
                local_fission_source += (self._scalar_fluxes[index]
                                         * mat._nu_sigma_f[g])
            end
            fission_source += local_fission_source * vol
        end
    end

    # Normalize fluxes by the total fission source
    for i = 1:self._num_FSRs * self._num_groups
        self._scalar_fluxes[i] /= fission_source
        self._old_scalar_fluxes[i] /= fission_source
    end
end


#=
    Calculates the FSR sources
=#
function computeFSRSources(self::StandardSolver)

    # Zero FSR sources
    for i = 1:self._num_FSRs * self._num_groups
        self._reduced_sources[i] = 0.0
    end

    # Compute each FSR source
    geometry = self._track_generator._geometry
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y

            # Get FSR material
            id = getFSRId(geometry, i, j)
            mat = geometry._mesh[i, j]

            # Compute fission source
            fission_source = 0.0
            for g = 1:self._num_groups
                index = (id-1) * self._num_groups + g
                fission_source += (self._scalar_fluxes[index]
                                   * mat._nu_sigma_f[g])
            end
            fission_source /= self._k_eff

            # Compute the scattering and total source
            for g = 1:self._num_groups

                index = (id-1) * self._num_groups + g

                in_scatter_source = 0.0
                for gp = 1:self._num_groups
                    index_p = (id-1) * self._num_groups + gp
                    in_scatter_source += (mat._sigma_s[gp, g] *
                                          self._scalar_fluxes[index_p])
                end

                self._reduced_sources[index] = (in_scatter_source +
                                                mat._chi[g] * fission_source)
                self._reduced_sources[index] /= (4 * pi * mat._sigma_t[g])
            end
        end
    end
end


#=
    Adds the source term contribution to the FSR scalar flux
=#
function addSourceToScalarFlux(self::StandardSolver)

    geometry = self._track_generator._geometry
    
    # Add source term and normalize flux for each FSR
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            vol = getFSRVolume(geometry, id)
            mat = geometry._mesh[i, j]
            for g = 1:self._num_groups
                index = (id-1) * self._num_groups + g
                self._scalar_fluxes[index] /= (vol * mat._sigma_t[g])
                self._scalar_fluxes[index] += (4 * pi * 
                                               self._reduced_sources[index])
            end
        end
    end
end


#=
    A function that computes the angular flux attenuation and its contribution
    to the local scalar flux for a given segment
=#
function tallyScalarFlux(self::StandardSolver, seg::Segment, azim::Int64, 
                         index::Int64)

    # Extract segment and material information
    sigma_t = self._materials[seg._id]._sigma_t
    ang_index = index
    sc_index = (seg._id - 1) * self._num_groups

    # Loop over groups
    for g = 1:self._num_groups
        
        tau = sigma_t[g] * seg._length
        sc_index += 1

        # Loop over polar angles
        for p = 1:self._num_polar_2
            ang_index += 1
            exponential = 1.0 - exp(-tau * self._inv_sin_thetas[p])
            delta_psi = exponential * (self._angular_fluxes[ang_index] -
                                       self._reduced_sources[sc_index])
            weight = self._weights[(azim-1) * self._num_polar_2 + p] 
            self._scalar_fluxes[sc_index] += delta_psi * weight
            self._angular_fluxes[ang_index] -= delta_psi
        end
    end
end


#=
    Loops over all tracks and segments, applying the MOC equations
=#
function transportSweep(self::StandardSolver)

    # Zero scalar fluxes
    flattenFSRFluxes(self, 0.0)

    # Loop over tracks
    index = 0
    track_count = 0
    for a = 1:self._track_generator._num_azim_2
        for i = 1:self._track_generator._num_tracks[a]

            # Extract track information
            track = self._track_generator._tracks[a][i]

            # Loop over segments forward
            track_count += 1
            for s = 1:length(track._segments)
                tallyScalarFlux(self, track._segments[s], a, index)
            end
            index += self._num_groups * self._num_polar_2

            # Loop over segments backward
            track_count += 1
            for s = length(track._segments):-1:1
                tallyScalarFlux(self, track._segments[s], a, index)
            end
            index += self._num_groups * self._num_polar_2
        end
    end
end


#=
    Computes the eigenvalue from fission rates
=#
function computeKeff(self::StandardSolver)

    # Compute the total fission source
    geometry = self._track_generator._geometry
    fission_source = 0.0
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            mat = geometry._mesh[i, j]
            vol = getFSRVolume(geometry, id)
            tmp_fission_source = 0.0
            for g = 1:self._num_groups
                index = (id-1) * self._num_groups + g
                tmp_fission_source += (self._scalar_fluxes[index]
                                       * mat._nu_sigma_f[g])
            end
            fission_source += tmp_fission_source * vol
        end
    end

    self._k_eff *= fission_source

end


#=
    Calculates the residual using the RMS difference in fission rates
=#
function computeResidual(self::StandardSolver)

    # Compute the RMS fission rate difference
    geometry = self._track_generator._geometry
    res = 0.0
    for i = 1:geometry._num_x
        for j = 1:geometry._num_y
            id = getFSRId(geometry, i, j)
            mat = geometry._mesh[i, j]
            old_fission_source = 0.0
            new_fission_source = 0.0
            for g = 1:self._num_groups
                index = (id-1) * self._num_groups + g
                old_fission_source += (self._old_scalar_fluxes[index]
                                       * mat._nu_sigma_f[g])
                new_fission_source += (self._scalar_fluxes[index]
                                       * mat._nu_sigma_f[g])
            end
            diff = new_fission_source / old_fission_source - 1.0
            res += diff * diff
        end
    end
    return sqrt(res / self._num_FSRs)
end


#=
    Saves the scalar flux to the old flux array
=#
function storeFSRFluxes(self::StandardSolver)
    for i = 1:self._num_FSRs * self._num_groups
        self._old_scalar_fluxes[i] = self._scalar_fluxes[i]
    end
end


#=
    Calculates the eigenvlaue of the problem represented by the TrackGenerator
=#
function computeEigenvalue(self::StandardSolver)

    # Allocate memory for angular and scalar fluxes
    if (!self._memory_allocated)
        allocateMemory(self)
    end

    geometry = self._track_generator._geometry
    flattenFSRFluxes(self, 1.0)
    storeFSRFluxes(self)

    # Source iteration loop
    println("[RESULT]\t\tIteration\tk-eff\t\tResidual")
    for iter = 1:self._max_iters

        # Setup sources and zero angular fluxes (assumes vacuum boundaries)
        normalizeFluxes(self)
        computeFSRSources(self)
        zeroAngularFluxes(self)

        # Compute new angular and scalar fluxes
        transportSweep(self)
        addSourceToScalarFlux(self)
        computeKeff(self)

        # Print the current eigenvalue
        res = computeResidual(self)
        storeFSRFluxes(self)
        println("[RESULT]", "\t\t", iter, "\t\t", round(self._k_eff,5),
                "\t\t", res)

        # Check for convergence
        if (iter > 2 && res < self._convergence_thresh)
            break
        end
    end
end
