using ArnoldiVandermonde
using Test
using Aqua

@testset "ArnoldiVandermonde.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ArnoldiVandermonde)
    end
    # Write your tests here.
end
