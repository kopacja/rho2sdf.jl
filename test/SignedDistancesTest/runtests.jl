module SignedDistancesTest

using Test

@testset "PseudoNormalsTest" begin include("PseudoNormalsTest.jl") end

@testset "sdfOnDensityFieldTest" begin include("sdfOnDensityFieldTest.jl") end

@testset "sdfOnTriangularMeshTest" begin include("sdfOnTriangularMeshTest.jl") end

end
