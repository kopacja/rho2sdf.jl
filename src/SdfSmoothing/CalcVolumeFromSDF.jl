"""
    calculate_volume_from_sdf(
        fine_sdf::Array{Float32,3}, 
        fine_grid::Array{Vector{Float32},3}; 
        iso_threshold::Float32=0.0f0, 
        detailed_quad_order::Int=9
    ) -> Float32

Calculates volume of geometry defined by SDF iso-surface.

# Arguments
- `fine_sdf`: SDF values at regular grid nodes
- `fine_grid`: Coordinates of grid nodes, each element is [x,y,z]
- `iso_threshold`: Iso-surface value (default = 0.0)
- `detailed_quad_order`: Gauss-Legendre quadrature order

# Returns
- Volume of the geometry

# Description
For each grid element:
1. Skip if outside iso-surface (all SDF values < threshold)
2. Add full volume if inside iso-surface (all SDF values ≥ threshold)
3. Use quadrature for partial elements intersecting the iso-surface
"""
function calculate_volume_from_sdf(
    fine_sdf::Array{Float32,3},
    fine_grid::Array{Vector{Float32},3};
    iso_threshold::Float32 = 0.0f0,
    detailed_quad_order::Int = 9
)
    # Determine grid dimensions
    nx, ny, nz = size(fine_sdf)
    @assert size(fine_grid) == (nx, ny, nz) "Dimensions of fine_sdf and fine_grid must match"
    
    # Calculate element size (assuming all elements are cubes with the same edge length)
    edge_vector = fine_grid[2, 1, 1] .- fine_grid[1, 1, 1]
    element_edge_length = norm(edge_vector)
    element_volume = element_edge_length^3
    
    # Set up Gauss-Legendre quadrature
    gp_double, w_double = FastGaussQuadrature.gausslegendre(detailed_quad_order)
    gp = Float32.(gp_double)
    w = Float32.(w_double)
    
    # Initialize total volume using atomic variable for parallel processing
    total_volume = Atomic{Float32}(0.0f0)
    
    # Precompute Jacobian for cubic element
    # For mapping from [-1,1]^3 to a cube with edge a, the Jacobian determinant is (a/2)^3 = a^3/8
    jacobian_det = element_volume / 8.0f0
    
    # Parallel iteration over all elements
    @threads for k in 1:nz-1
        for j in 1:ny-1
            for i in 1:nx-1
                # Get SDF values at element nodes
                c000 = fine_sdf[i, j, k]
                c100 = fine_sdf[i+1, j, k]
                c010 = fine_sdf[i, j+1, k]
                c110 = fine_sdf[i+1, j+1, k]
                c001 = fine_sdf[i, j, k+1]
                c101 = fine_sdf[i+1, j, k+1]
                c011 = fine_sdf[i, j+1, k+1]
                c111 = fine_sdf[i+1, j+1, k+1]
                
                # Check if element can contain the iso-surface
                min_value = min(c000, c100, c010, c110, c001, c101, c011, c111)
                max_value = max(c000, c100, c010, c110, c001, c101, c011, c111)
                
                # Skip element if completely outside the iso-surface
                if max_value < iso_threshold
                    continue
                end
                
                # If element is completely inside the iso-surface, add its full volume
                if min_value >= iso_threshold
                    atomic_add!(total_volume, element_volume)
                    continue
                end
                
                # Element intersects the iso-surface, perform numerical integration
                element_volume_partial = 0.0f0
                
                # Perform numerical integration over all quadrature points
                for kq in 1:detailed_quad_order
                    zeta_m = gp[kq]  # Gauss point in [-1,1]
                    zeta = (zeta_m + 1) / 2  # Convert to [0,1] for interpolation
                    
                    for jq in 1:detailed_quad_order
                        eta_m = gp[jq]
                        eta = (eta_m + 1) / 2
                        
                        for iq in 1:detailed_quad_order
                            xi_m = gp[iq]
                            xi = (xi_m + 1) / 2
                            
                            # Trilinear interpolation of SDF value at Gauss point
                            c00 = c000 * (1.0f0 - xi) + c100 * xi
                            c01 = c001 * (1.0f0 - xi) + c101 * xi
                            c10 = c010 * (1.0f0 - xi) + c110 * xi
                            c11 = c011 * (1.0f0 - xi) + c111 * xi
                            
                            c0 = c00 * (1.0f0 - eta) + c10 * eta
                            c1 = c01 * (1.0f0 - eta) + c11 * eta
                            
                            point_sdf = c0 * (1.0f0 - zeta) + c1 * zeta
                            
                            # Add contribution only if point is inside (or on) the iso-surface
                            if point_sdf >= iso_threshold
                                weight = w[iq] * w[jq] * w[kq]
                                element_volume_partial += weight * jacobian_det
                            end
                        end
                    end
                end
                
                # Add to total volume
                atomic_add!(total_volume, element_volume_partial)
            end
        end
    end
    
    return total_volume[]
end

# calculate_volume_from_sdf(fine_sdf, fine_grid)

function calculate_volume_simple(
    fine_sdf::Array{Float32,3},
    fine_grid::Array{Vector{Float32},3};
    iso_threshold::Float32 = 0.0f0
)
    # Determine grid dimensions
    nx, ny, nz = size(fine_sdf)
    @assert size(fine_grid) == (nx, ny, nz) "Dimensions of fine_sdf and fine_grid must match"
    
    # Calculate element size
    edge_vector = fine_grid[2, 1, 1] .- fine_grid[1, 1, 1]
    element_edge_length = norm(edge_vector)
    element_volume = Float32(element_edge_length^3)
    
    # Initialize total volume
    total_volume = 0.0f0
    
    # Counters for statistics
    full_elements = 0
    partial_elements = 0
    empty_elements = 0
    
    # Iterate over all elements
    for k in 1:nz-1
        for j in 1:ny-1
            for i in 1:nx-1
                # Get SDF values at element nodes
                values = [
                    fine_sdf[i, j, k],      # c000
                    fine_sdf[i+1, j, k],    # c100
                    fine_sdf[i, j+1, k],    # c010
                    fine_sdf[i+1, j+1, k],  # c110
                    fine_sdf[i, j, k+1],    # c001
                    fine_sdf[i+1, j, k+1],  # c101
                    fine_sdf[i, j+1, k+1],  # c011
                    fine_sdf[i+1, j+1, k+1] # c111
                ]
                
                # Count nodes inside material (SDF >= iso_threshold)
                nodes_inside = count(v -> v >= iso_threshold, values)
                
                if nodes_inside == 8
                    # All nodes inside material - add full volume
                    total_volume += element_volume
                    full_elements += 1
                elseif nodes_inside > 0
                    # Some nodes inside - add proportional volume
                    total_volume += (nodes_inside / 8.0f0) * element_volume
                    partial_elements += 1
                else
                    # No nodes inside
                    empty_elements += 1
                end
            end
        end
    end
    
    # Print statistics
    println("Volume calculation statistics:")
    println("  Fully inside elements: $full_elements")
    println("  Partially inside elements: $partial_elements")
    println("  Completely outside elements: $empty_elements")
    println("  Total elements: $(full_elements + partial_elements + empty_elements)")
    println("  Calculated volume: $total_volume")
    
    return total_volume
end

# volume = calculate_volume_simple(fine_sdf, fine_grid)
