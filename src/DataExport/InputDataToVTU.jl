function InputDataToVTU(mesh::Mesh, filename::String)

  # Extract node coordinates
  points = mesh.X

  # Check the dimensionality of the coordinates
  nsd_actual = size(points, 1)

  # If mesh is 2D but we need 3D coordinates for VTK
  if nsd_actual == 2
      points = vcat(points, zeros(Float64, 1, size(points, 2)))
  end

  # Determine the VTK cell type based on the number of element nodes
  vtk_cell_type = if mesh.nen == 8
      VTKCellTypes.VTK_HEXAHEDRON  # 8-node hexahedron
  elseif mesh.nen == 4 && nsd_actual == 3
      VTKCellTypes.VTK_TETRA  # 4-node tetrahedron
  elseif mesh.nen == 4 && nsd_actual == 2
      VTKCellTypes.VTK_QUAD  # 4-node quadrilateral
  elseif mesh.nen == 3
      VTKCellTypes.VTK_TRIANGLE  # 3-node triangle
  else
      error("Unsupported element type: nen=$(mesh.nen), nsd=$(nsd_actual)")
  end

  # Create cell array for VTK (convert from 1-based to 0-based indexing)
  cells = [MeshCell(vtk_cell_type, Vector(mesh.IEN[:, i])) for i in 1:mesh.nel]

  # Create the VTK unstructured grid
  vtkfile = vtk_grid(filename, points, cells)

  # Add element density as a cell field
  vtk_cell_data(vtkfile, mesh.rho, "density")

  # Write the file
  vtk_save(vtkfile)

  println("VTU file written successfully: $(filename).vtu")

  return vtkfile
end
