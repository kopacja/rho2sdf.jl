module Rho2sdf

using LinearAlgebra
using Statistics
using DelimitedFiles
using Einsum
using BenchmarkTools
using ProgressMeter

# Predefined shape functions and its derivatives:
include("ShapeFunctions/ShapeFunctions.jl")
using .ShapeFunctions

include("MeshGrid/MeshGrid.jl")
using .MeshGrid

include("SignedDistances/SignedDistances.jl")
using .SignedDistances

include("PrimitiveGeometries/PrimitiveGeometries.jl")
using .PrimitiveGeometries

include("MarchingCubes/MarchingCubes.jl")
using .MarchingCubes


include("DataExport/DataExport.jl")
using .DataExport

function exportToVTU(
    fileName::String,
    X::Vector{Vector{Float64}},
    IEN::Vector{Vector{Int64}},
    VTK_CODE::Int64,
    rho::Vector{Float64}=nothing,
)

    nnp = length(X)
    nsd = length(X[1])

    nel = length(IEN)
    nen = length(IEN[1])

    io = open(fileName, "w")

    println(
        io,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">",
    )
    println(io, "  <UnstructuredGrid>")
    println(
        io,
        "    <Piece NumberOfPoints=\"",
        nnp,
        "\" NumberOfCells=\"",
        nel,
        "\">",
    )

    println(io, "	  <Points>")
    println(
        io,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">",
    )
    for A = 1:nnp
        print(io, "          ")
        for i = 1:nsd
            if (abs(X[A][i]) < 1.0e-20)
                X[A][i] = 0.0
            end
            print(io, " ", X[A][i])
        end
        print(io, "\n")
    end
    println(io, "        </DataArray>")
    println(io, "	  </Points>")

    println(io, "      <Cells>")
    println(
        io,
        "		  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">",
    )

    for el = 1:nel
        print(io, "         ")
        for a = 1:nen
            print(io, " ", IEN[el][a] - 1)
        end
        print(io, "\n")
    end

    println(io, "        </DataArray>")
    println(
        io,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">",
    )
    for i = nen:nen:nen*nel
        println(io, "          ", i)
    end
    println(io, "        </DataArray>")
    println(
        io,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">",
    )
    for el = 1:nel
        println(io, "          ", VTK_CODE)
    end
    println(io, "        </DataArray>")
    println(io, "      </Cells>")

    if (rho !== nothing)
    println(io, "      <PointData Scalars=\"scalars\">")
    println(io, "           <DataArray type=\"Float32\" Name=\"density\" Format=\"ascii\">")
    for A = 1:nnp
        println(io, "             ", rho[A])
    end
    println(io, "           </DataArray>")
    println(io, "      </PointData>")
    end


    println(io, "    </Piece>")
    println(io, "  </UnstructuredGrid>")
    println(io, "</VTKFile>")

    close(io)
end

function exportStructuredPointsToVTK(
    fileName::String,
    grid::Grid,
    vals::Vector{Float64},
    valLabel::String,
)

    dim = grid.N .+ 1
    org = grid.AABB_min
    spacing = grid.cell_size

    io = open(fileName, "w")

    write(io, "# vtk DataFile Version 1.0\n")
    write(
        io,
        "Texture map for thresholding data (use boolean textures for 2D map)\n",
    )

    write(io, "ASCII\n\n")
    write(io, "DATASET STRUCTURED_POINTS\n")
    dim_x = dim[1]
    dim_y = dim[2]
    dim_z = dim[3]
    write(io, "DIMENSIONS $dim_x $dim_y $dim_z\n") # dimenze pravidelné sítě
    
    #     spacing_x = spacing[1]
    #     spacing_y = spacing[2]
    #     spacing_z = spacing[3]

    #     write(io, "SPACING $spacing_x $spacing_y $spacing_z\n") # krok sítě ve 3 směrech jiný

    write(io, "SPACING $spacing $spacing $spacing\n") # krok sítě ve 3 směrech stejný

    org_x = org[1]
    org_y = org[2]
    org_z = org[3]
    write(io, "ORIGIN $org_x $org_y $org_z\n\n") # souřadnice počátku

    n = prod(dim)
    write(io, "POINT_DATA $n\n") # počet uzlů pravidelné sítě
    write(io, "SCALARS $valLabel float 1\n") # druh hodnoty v uzlech (vzdálenost)
    write(io, "LOOKUP_TABLE default\n") # ??

    for val ∈ vals
        write(io, "$val\n") # hodnota vzdálenostní funkce v jednotlivých bodech
    end
    close(io)
end



end # module Rho2sdf
