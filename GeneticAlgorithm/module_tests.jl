# Tests for NDGeneticAlgorithm module

using NDGeneticAlgorithm,
      Test

@testset "General module tests" begin

    @testset "Tests for `Individual` struct" begin
        ### Test `Individual` instantiation
        a = Individual(5)
        b = Individual(BitArray(zeros(5)), -Inf)

        # Equality of variables
        @test a.genome == b.genome
        @test a.fitness == b.fitness

        # Equality of `Individual`s
        @test a == b
    end

    @testset "Tests for `crossover`" begin
        c = Individual(BitArray([0,0,0,1,1,1]))
        d = Individual(BitArray([1,1,1,0,0,0]))
        off1, off2 = crossover(c, d, x_point=4)

        @test off1.genome == BitArray(zeros(6))
        @test off2.genome == BitArray(ones(6))
    end

    @testset "Tests for mutation!" begin
        a = Individual(5)
        NDGeneticAlgorithm.mutate!(a, BitArray([1,0,0,0,0]))

        @test a.genome == BitArray([1,0,0,0,0])

        b = Individual(5)
        mutation!(b, 0.7)

        @test typeof(b) == Individual
        @test typeof(b.genome) == BitArray{1}
    end

    @testset "Tests for fitness" begin

    end
end
