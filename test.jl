# Import the required classes
include("src/Material.jl")
include("src/Track.jl")
include("src/Geometry.jl")
include("src/TrackGenerator.jl")
include("src/PolarQuad.jl")
include("src/StandardSolver.jl")
include("src/SparseSolver.jl")

# Define material properties
fuel = Material(1)
fuel._sigma_t = [0.452648699]
fuel._sigma_s = 0.383259177 * eye(1)
fuel._sigma_f = [0.0414198575]
fuel._nu_sigma_f = [0.0994076580]
fuel._chi = [1.0]

# Create the mesh
nx = 3
ny = 3
mesh = Array{Material}(nx*ny)
for i = 1:nx
    for j = 1:ny
        mesh[nx*(i-1) + j] = fuel
    end
end
mesh = reshape(mesh, (nx,ny))

# Create the geometry
g = Geometry(-100,100,-100,100)
setMesh(g, mesh)

# Create the TrackGenerator
T = TrackGenerator(g)
generateTracks(T, 10.0, 4)
traceTracks(T)

# Create the Solver and compute the eigenvalue
#S = StandardSolver(T, 3)
S = SparseSolver(T, 3)
computeEigenvalue(S)
