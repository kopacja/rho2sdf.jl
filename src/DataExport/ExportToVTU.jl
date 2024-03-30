
function exportToVTU(
    fileName::String,
    X::Vector{Vector{Float64}},
    IEN::Vector{Vector{Int64}},
    VTK_CODE::Int64,
    rho::Union{Vector{Float64}, Nothing} = nothing
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


