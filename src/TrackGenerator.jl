using PyPlot
using Colors

#=
    A TrackGenerator contains information about the Tracks laid across the
    given Geometry for use in an MOC solver
=#
type TrackGenerator

    # The number of azimuthal angles in [0, π]
    _num_azim_2::Int64

    # The number of Tracks intersecting the x-axis for each azimuthal angle
    _num_x::Array{Int64}
    
    # The number of Tracks intersecting the y-axis for each azimuthal angle
    _num_y::Array{Int64}

    # The total number of Tracks laid down for each azimuthal angle
    _num_tracks::Array{Int64}

    # The Tracks laid down across the Geometry
    _tracks::Array{Array{Track}}

    # The azimuthal angles for which Tracks are laid down across the Geometry
    _phi::Array{Float64}

    # The perpendicular spacing between Tracks for each azimuthal angle
    _azim_spacings::Array{Float64}

    # The weight of each azimuthal angle
    _azim_weights::Array{Float64}

    # The Geometry across which Tracks are laid down
    _geometry::Geometry

    # The constructor for the TrackGenerator
    TrackGenerator(geometry::Geometry) = new(0, 
                                             Array(Int64,0), 
                                             Array(Int64,0),
                                             Array(Int64,0),
                                             Array(Array{Track}, 0),
                                             Array(Float64,0),
                                             Array(Float64,0),
                                             Array(Float64,0),
                                             geometry)

end


#= 
    Lays Tracks across the associated Geometry with a defined azimuthal ray
    spacing and number of azimuthal angles
=#
function generateTracks(self::TrackGenerator, azim_spacing::Float64,
                        num_azim::Int64)

    # Initialize arrays
    self._num_azim_2 = floor(num_azim / 2)
    self._num_x = Array{Int64}(self._num_azim_2)
    self._num_y = Array{Int64}(self._num_azim_2)
    self._num_tracks = Array{Int64}(self._num_azim_2)
    self._tracks = Array{Array{Track}}(self._num_azim_2)
    self._phi = Array{Float64}(self._num_azim_2)
    self._azim_spacings = Array{Float64}(self._num_azim_2)
    self._azim_weights = Array{Float64}(self._num_azim_2)

    # Extract geometry infromation
    width_x = getWidthX(self._geometry)
    width_y = getWidthY(self._geometry)
    min_x = self._geometry._min_x
    min_y = self._geometry._min_y

    # Create tracks
    for a = 1:convert(Int64, self._num_azim_2/2)

        # Compute the desired azimuthal angle
        phi = pi / self._num_azim_2 * (a - 0.5)
        
        # Compute the number of intersections with x and y axes
        self._num_x[a] = floor(width_x / azim_spacing * sin(phi)) + 1
        self._num_y[a] = floor(width_x / azim_spacing * cos(phi)) + 1

        # Record the total number of tracks for this azimuthal angle
        self._num_tracks[a] = self._num_x[a] + self._num_y[a]

        # Compute and record the corrected azimuthal angle
        phi = atan((width_y * self._num_x[a]) / (width_x * self._num_y[a]))
        self._phi[a] = phi

        # Record complimentary angle information
        self._num_x[self._num_azim_2-a+1] = self._num_x[a]
        self._num_y[self._num_azim_2-a+1] = self._num_y[a]
        self._num_tracks[self._num_azim_2-a+1] = self._num_tracks[a]
        self._phi[self._num_azim_2-a+1] = pi - phi

        # Compute ray spacing
        dx_eff = width_x / self._num_x[a]
        dy_eff = width_y / self._num_y[a]

        # Compute perpendicular ray spacing
        self._azim_spacings[a] = dx_eff * sin(phi)
        self._azim_spacings[self._num_azim_2-a+1] = dx_eff * sin(phi)

        # Allocate coordinates
        self._tracks[a] = Array{Track}(self._num_tracks[a])
        self._tracks[self._num_azim_2-a+1] = Array{Track}(self._num_tracks[a])

        # Compute starting points for tracks starting on the x-axis
        for i = 1:self._num_x[a]
            self._tracks[a][i] = Track()
            self._tracks[a][i]._coords = (min_x + dx_eff *
                                          (self._num_x[a] - i + 0.5), min_y)
            self._tracks[self._num_azim_2-a+1][i] = Track()
            self._tracks[self._num_azim_2-a+1][i]._coords = (min_x + dx_eff *
                                                             (i - 0.5), min_y)
        end
        
        # Compute starting points for tracks starting on the y-axis
        for i = 1:self._num_y[a]
            ind = self._num_x[a] + i
            self._tracks[a][ind] = Track()
            self._tracks[a][ind]._coords = (min_x, min_y + dy_eff * (i - 0.5))
            self._tracks[self._num_azim_2-a+1][ind] = Track()
            self._tracks[self._num_azim_2-a+1][ind]._coords = (min_x + width_x,
                                                               min_y + dy_eff *
                                                               (i - 0.5))
        end
    end

    # Compute azimuthal weights
    for a = 1:convert(Int64, self._num_azim_2/2)
        
        phi_pos = 0
        phi_neg = 0

        if a == 1
            phi_neg = self._phi[a]
        else
            phi_neg = self._phi[a] - self._phi[a-1]
        end

        if a == self._num_azim_2/2
            phi_pos = pi / 2 - self._phi[a]
        else
            phi_pos = self._phi[a+1] - self._phi[a]
        end

        weight = (phi_pos + phi_neg) / pi
        self._azim_weights[a] = weight
        self._azim_weights[self._num_azim_2-a+1] = weight

    end
end


#=
    Traces the Tracks to determine and allocate their associated Segments
=#
function traceTracks(self::TrackGenerator)
    
    # Determine the width of each cell in both the x and y directions
    total_width_x = getWidthX(self._geometry)
    total_width_y = getWidthY(self._geometry)
    width_x = total_width_x / self._geometry._num_x
    width_y = total_width_y / self._geometry._num_y

    # Loop over all azimuthal angles
    for a = 1:self._num_azim_2
      
        # Determine the direction of the Tracks
        unit_x = cos(self._phi[a])
        unit_y = sin(self._phi[a])

        # Determine the local boundaries of the cell in the x-direction
        bound_x = Inf
        reset_x = Inf
        dir_x = 0
        if unit_x > 0
            bound_x = width_x
            reset_x = 0.0
            dir_x = 1.0
        elseif unit_x < 0
            bound_x = 0.0
            reset_x = width_x
            dir_x = -1.0
        end
        
        # Determine the local boundaries of the cell in the y-direction
        bound_y = Inf
        reset_y = Inf
        dir_y = 0
        if unit_y > 0
            bound_y = width_y
            reset_y = 0.0
            dir_y = 1.0
        elseif unit_y < 0
            bound_y = 0.0
            reset_y = width_y
            dir_y = -1.0
        end

        # Loop over all Tracks for the current azimuthal angle
        for i = 1:self._num_tracks[a]
        
            # Determine the starting x and y coordinates
            x_rel = self._tracks[a][i]._coords[1] - self._geometry._min_x
            y_rel = self._tracks[a][i]._coords[2] - self._geometry._min_y

            # Determine the 
            cell_x = floor(x_rel / width_x) + 1
            cell_y = floor(y_rel / width_y) + 1

            # Determine if the Track is on a max x-boundary
            dist_to_edge = abs(x_rel - total_width_x)
            if (cell_x == self._geometry._num_x+1 && dist_to_edge < 1e-8)
                cell_x = self._geometry._num_x
            end
            
            # Determine if the Track is on a max y-boundary
            dist_to_edge = abs(y_rel - total_width_y)
            if (cell_y == self._geometry._num_y+1 && dist_to_edge < 1e-8)
                cell_y = self._geometry._num_y
            end

            # Check to ensure the x and y cells are legitamate
            if cell_x < 0 || cell_x > self._geometry._num_x
                println("Found Cell ", cell_x, " at ", self._coords[a][i])
            end 
            if cell_y < 0 || cell_y > self._geometry._num_y
                println("Found Cell ", cell_y, " at ", self._coords[a][i])
            end
            
            # Determine the local x and y coordinates in the current cell
            x_rel -= (cell_x - 1) * width_x
            y_rel -= (cell_y - 1) * width_y

            # Trace the Track while inside the Geometry
            while cell_x > 0 && cell_x <= self._geometry._num_x &&
                  cell_y > 0 && cell_y <= self._geometry._num_y

                # Find the current FSR ID
                fsr_id = getFSRId(self._geometry, cell_x, cell_y)

                # Calculate the distances to x and y intersections
                dist = 0.0
                dist_x = (bound_x - x_rel) / unit_x
                dist_y = (bound_y - y_rel) / unit_y
                
                # Determine which distance is closer
                if (dist_x < dist_y)
                    dist = dist_x
                    x_rel = reset_x
                    y_rel += unit_y * dist
                    cell_x += dir_x
                else
                    dist = dist_y
                    x_rel += unit_x * dist
                    y_rel = reset_y
                    cell_y += dir_y
                end

                # Create the Segment
                seg = Segment(fsr_id, dist)
                push!(self._tracks[a][i]._segments, seg)
            end
        end
    end
end


#=
    Plots all Tracks laid across the Geometry
=#
function plotTracks(self::TrackGenerator)
    
    surf_x = [self._geometry._min_x, self._geometry._max_x]
    surf_y = [self._geometry._min_y, self._geometry._max_y]

    # Loop over all azimuthal angles
    for azim = 1:self._num_azim_2
        track_ends = Array{Tuple{Float64, Float64}}(self._num_tracks[azim])
        
        # Calculate the direction of the Tracks
        unit_x = cos(self._phi[azim])
        unit_y = sin(self._phi[azim])

        # Loop over all Tracks with this azimuthal angle
        for i = 1:self._num_tracks[azim]

            # Calculate the distance to an x-boundary crossing
            distance_x = Inf
            if unit_x < 0
                distance_x = (self._geometry._min_x - 
                              self._tracks[azim][i]._coords[1])
                distance_x /= unit_x
            elseif unit_x > 0
                distance_x = (self._geometry._max_x - 
                              self._tracks[azim][i]._coords[1])
                distance_x /= unit_x
            end
            
            # Calculate the distance to an y-boundary crossing
            distance_y = Inf
            if unit_y < 0
                distance_y = (self._geometry._min_y - 
                              self._tracks[azim][i]._coords[2])
                distance_y /= unit_y
            elseif unit_y > 0
                distance_y = (self._geometry._max_y - 
                              self._tracks[azim][i]._coords[2])
                distance_y /= unit_y
            end

            dist = min(distance_x, distance_y)

            track_ends[i] = (dist*unit_x + self._tracks[azim][i]._coords[1], 
                             dist*unit_y + self._tracks[azim][i]._coords[2])
        end

        # Plot the Tracks across the geometry
        for i = 1:self._num_tracks[azim]
            start = self._tracks[azim][i]._coords
            term = track_ends[i]
            plot([start[1], term[1]], [start[2], term[2]], "k")
            plot([start[1]], [start[2]], "bo")
            plot([term[1]], [term[2]], "ro")
        end
    end
    axis("tight")
    show()
end


#=
    Plots all Segments laid across the Geometry
=#
function plotSegments(self::TrackGenerator)
    
    # Plot the boundaries of the cells
    width_x = getWidthX(self._geometry) / self._geometry._num_x
    width_y = getWidthY(self._geometry) / self._geometry._num_y
    num_x = self._geometry._num_x
    num_y = self._geometry._num_y
    for i = 1:num_x
        for j = 1:num_y

            surf_x = [self._geometry._min_x + (i-1)*width_x, 
                      self._geometry._min_x + i*width_x]
            surf_y = [self._geometry._min_y + (j-1)*width_y,
                      self._geometry._min_y + j*width_y]
    
            plot([surf_x[1], surf_x[1]], [surf_y[1], surf_y[2]], "k")
            plot([surf_x[2], surf_x[2]], [surf_y[1], surf_y[2]], "k")
            plot([surf_x[1], surf_x[2]], [surf_y[1], surf_y[1]], "k")
            plot([surf_x[1], surf_x[2]], [surf_y[2], surf_y[2]], "k")
        end
    end

    # Create a map to color by cell
    color_map = distinguishable_colors(num_x*num_y)
    
    # Get the Geometry boundaries
    surf_x = [self._geometry._min_x, self._geometry._max_x]
    surf_y = [self._geometry._min_y, self._geometry._max_y]

    # Plot the segments
    for azim = 1:self._num_azim_2
        track_ends = Array{Tuple{Float64, Float64}}(self._num_tracks[azim])
        
        unit_x = cos(self._phi[azim])
        unit_y = sin(self._phi[azim])

        for i = 1:self._num_tracks[azim]

            x_curr = self._tracks[azim][i]._coords[1]
            y_curr = self._tracks[azim][i]._coords[2]
            for s = 1:length(self._tracks[azim][i]._segments)
                seg = self._tracks[azim][i]._segments[s]
                id = seg._id
                dist = seg._length
                rgb = color_map[id]
                c = [rgb.r, rgb.g, rgb.b]
                x_end = x_curr + unit_x * dist
                y_end = y_curr + unit_y * dist
                plot([x_curr, x_end], [y_curr, y_end], color=c)
                plot([x_curr, x_end], [y_curr, y_end], "ko")
                x_curr = x_end
                y_curr = y_end
            end
        end
    end
    axis("tight")
    show()
end

#=
    Returns the total number of Tracks laid across the geometry
=#
function getNumTotalTracks(self::TrackGenerator)
    num_tracks = 0
    for a = 1:convert(Int64, self._num_azim_2/2)
        num_tracks += 2*self._num_tracks[a]
    end
    return num_tracks
end
