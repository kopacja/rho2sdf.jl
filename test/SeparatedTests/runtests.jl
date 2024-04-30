module SeperatedTests

using Test

# @testset "TestsOnSingleElement" begin include("TestsOnSingleElement.jl") end
@testset "TestOnLegGripper" begin include("TestOnLegGripper.jl") end
@testset "TestOnLegGripperSTL" begin include("TestOnLegGripperSTL.jl") end

end
