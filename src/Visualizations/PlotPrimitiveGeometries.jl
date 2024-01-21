
using Makie

# Sample data
black_points = rand(Point3f0, 100) # Random set of points
red_points = [Point3f0(x,y,z) for x in 0:1, y in 0:1, z in 0:1] # Points of a cube

# Create a scene
scene = Scene()

# Plotting the first set of points (black and small)
scatter!(scene, black_points, color = :black, markersize = 0.05)

# Plotting the second set of points (red and large)
scatter!(scene, red_points, color = :red, markersize = 0.2)

# Creating a wireframe model of the cube
lines!(scene, red_points[1:4], color = :red)
lines!(scene, red_points[5:8], color = :red)
for i in 1:4
    lines!(scene, [red_points[i], red_points[i+4]], color = :red)
end

# Display the scene
scene
