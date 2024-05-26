import open3d as o3d
import numpy as np
import pandas as pd
import time

# Načtení souřadnic bodů sítě z CSV souboru
points = pd.read_csv('sdf_grid_points.csv', header=None).values

# Výstup tvaru matice před transpozicí
print(f"Shape of points before transpose: {points.shape}")

# Transpozice dat, aby tvar byl (216, 3)
points = points.T

# Výstup tvaru matice po transpozici
print(f"Shape of points after transpose: {points.shape}")

# Nahrajte STL soubor
mesh = o3d.io.read_triangle_mesh('chapadloSTL.stl')
mesh.compute_vertex_normals()

# Převod na Tensor geometrie Open3D
mesh_t = o3d.t.geometry.TriangleMesh.from_legacy(mesh)

# Vytvoření RaycastingScene a přidání mesh
scene = o3d.t.geometry.RaycastingScene()
_ = scene.add_triangles(mesh_t)

# Ujistěte se, že query_points mají správný tvar (num_points, 3)
if points.shape[1] != 3:
    raise ValueError(f"Expected points to have shape (num_points, 3) but got {points.shape}")

query_points = o3d.core.Tensor(points, dtype=o3d.core.Dtype.Float32)

# Změřte čas potřebný pro raycasting
start_time = time.time()
# Vypočítejte occupancy, aby se zjistilo, zda jsou body uvnitř mesh
occupancy = scene.compute_occupancy(query_points)
end_time = time.time()

# Přiřazení +1 bodům uvnitř a -1 bodům vně tělesa
occupancy_results = np.where(occupancy.numpy() == 1, 1, -1)

# Výstup výsledků
np.savetxt('occupancy_results.csv', occupancy_results, delimiter=',')
print(f"Raycasting took {end_time - start_time:.4f} seconds")
