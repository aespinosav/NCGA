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

    @testset "Tests for selection" begin
        # Fitness proportional selection

        # Tournament selection (ordinal selection)
        a = Individual(3)
        a.fitness = 10
        b = Individual(3)
        b.fitness = 1
        victor = tournament_selection([a, b], 2)

        @test victor == a
    end

    @testset "Tests for fitness" begin

    end

    @testset "Tests for elitism" begin
        pop1 = [Individual(3) for i in 1:10]
        pop2 = [Individual(3) for i in 1:10]

        for i in 1:10
            pop1[i].fitness = 20 + i
            pop2[i].fitness = i
        end

        fitness_sort!(pop1)
        fitness_sort!(pop2)

        new_pop = elite_conservation(pop1, pop2, 0.1)

        @test new_pop[1].fitness == 30
        @test new_pop[2].fitness == 10
    end
end
