using ArnoldiVandermonde, LinearAlgebra
using Test
using Aqua

@testset "ArnoldiVandermonde.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ArnoldiVandermonde, ambiguities=false)
    end

    @testset "ArnoldiBasis construction" begin
        for z in (rand(10), rand(11) .+ im*rand(11))
            n = length(z)
            m = 5
            B = ArnoldiBasis(z, m)
            for _ = 1:3
                Q = vectors(B)
                @test nodes(B) == z
                @test size(Q) == (n, m+1)
                @test norm(Q' * Q - n * I) < n * 1e-10
                @test B.H isa UpperHessenberg
                increment!(B)
                m += 1
            end
        end
    end

    @testset "ArnoldiPolynomial construction" begin
        z = rand(10)
        B = ArnoldiBasis(z, 5)
        coeff = rand(6)
        p = ArnoldiPolynomial(coeff, B)
        @test p.coeff == coeff
        @test length(p) == length(coeff)
        @test degree(p) == 5
    end

    @testset "ArnoldiPolynomial evaluation" begin
        z = rand(10)
        B = ArnoldiBasis(z, 5)
        coeff = rand(6)
        p = ArnoldiPolynomial(coeff, B)
        for i in 1:length(z)
            @test p(z[i]) ≈ sum(coeff[j] * B.Q[i, j] for j in 1:6)
        end
        y = p(z)
        @test y ≈ B.Q[:, 1:6] * coeff
    end

    @testset "Approximation" begin
        f(x) = exp(x)
        x = range(0, 1, 40)
        B = ArnoldiBasis(x, 8)
        p = B \ f
        @test p.(x) ≈ f.(x) atol=1e-9
        p = project(f, 0, 1; tol=1e-13)
        x = nodes(p)
        @test p.(x) ≈ f.(x) atol=1e-13
        p = project(f, BigFloat(0), 1; tol=1e-50)
        x = nodes(p)
        @test p(x) ≈ f.(x) atol=1e-48
        # force a node refinement, off-grid evaluation
        f(x) = sin(16x)
        p = project(f, 0, 1; tol=1e-13)
        @test all(abs(p(x) - f(x)) < 1e-13 for x in range(0, 1, 1000))
    end

    @testset "Type promotion" begin
        z = rand(10)
        B = ArnoldiBasis(z, 5)
        coeff = rand(6)
        p = ArnoldiPolynomial(coeff, B)
        p32 = ArnoldiPolynomial{Float32}(p)
        @test eltype(p32) == Float32
        @test eltype(p32.basis) == Float32
    end

    @testset "Error handling" begin
        @test_throws ArgumentError ArnoldiBasis(rand(5), 5)
        @test_throws ArgumentError ArnoldiPolynomial(rand(5), ArnoldiBasis(rand(10), 5))
    end
end
