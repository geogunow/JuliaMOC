#=
    A class describing a uniform grid geometry
=#
type Geometry
    
    # The minimum x-coordinate in the geometry
    _min_x::Float64
    
    # The maximum x-coordinate in the geometry
    _max_x::Float64

    # The minimum y-coordinate in the geometry
    _min_y::Float64
    
    # The maximum y-coordinate in the geometry
    _max_y::Float64

    # The number of mesh cells in the x-direction
    _num_x::Int64
    
    # The number of mesh cells in the y-direction
    _num_y::Int64

    # The mesh containing material information in the geometry
    _mesh::Array{Material,2}

    # The constructor for the Geometry
    Geometry(min_x, max_x, min_y, max_y) = new(min_x, max_x, min_y, max_y, 1,
                                               1, Array{Material,2}(0,0))
end


#=
    Returns the width of the entire geometry in the x-direction
=#
function getWidthX(self::Geometry)
    return self._max_x - self._min_x
end


#=
Returns the width of the entire geometry in the y-direction
=#
function getWidthY(self::Geometry)
    return self._max_y - self._min_y
end


#=
    Sets the mesh and associated materials in the geometry
=#
function setMesh(self::Geometry, mesh::Array{Material, 2})
    self._num_x, self._num_y = size(mesh)
    self._mesh = mesh
end


#=
    Returns the number of energy groups in the problem and checks that
    all materials in the mesh have the same number of energy groups
=#
function getNumEnergyGroups(self::Geometry)
    num_energy_groups = self._mesh[1, 1]._num_groups
    for i = 1:self._num_x
        for j = 1:self._num_y
            if self._mesh[i,j]._num_groups != num_energy_groups
                println("ERROR: energy group mismatch at ", i, ", ", j)
            end
        end
    end
    return num_energy_groups
end


#=
    Returns the FSR ID from x and y indexes of the FSR
=#
function getFSRId(self::Geometry, cell_x, cell_y)
    return convert(Int64, (cell_y-1) * self._num_x + cell_x)
end


#=
    Returns the FSR volume associated with a given FSR ID
=#
function getFSRVolume(self::Geometry, id::Int64)
    return getWidthX(self) * getWidthY(self) / getNumCells(self)
end


#=
    Returns the total number of cells (FSRs) within the geometry
=#
function getNumCells(self::Geometry)
    return self._num_x * self._num_y
end
