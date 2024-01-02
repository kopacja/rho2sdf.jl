module MeshGridTest

using Test

@testset "InputDataTests" begin include("InputDataTests.jl") end
@testset "GridTests" begin include("GridTests.jl") end

end
