# Import the required classes
include("src/Material.jl")
include("src/Track.jl")
include("src/Geometry.jl")
include("src/TrackGenerator.jl")
include("src/PolarQuad.jl")
include("src/StandardSolver.jl")

# Define material properties
fuel = Material(1)
fuel._sigma_t = [0.452648699]
fuel._sigma_s = 0.383259177 * eye(1)
fuel._sigma_f = [0.0414198575]
fuel._nu_sigma_f = [0.0994076580]
fuel._chi = [1.0]

# Create the mesh (1x1)
mesh = [fuel]
mesh = reshape(mesh, (1,1))

# Create the geometry
g = Geometry(-100,100,-100,100)
setMesh(g, mesh)

# Create the TrackGenerator
T = TrackGenerator(g)
generateTracks(T, 0.1, 4)
traceTracks(T)

# Create the Solver and compute the eigenvalue
S = StandardSolver(T, 3)
computeEigenvalue(S)
