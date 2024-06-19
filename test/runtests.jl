using PriceIndexes
using Test

p0 = [2.0, 1.0, 1.0, 0.5]
p1 = [1.75, 0.5, 0.95, 0.55]
q0 = [0.384615385, 1.538461538, 1.538461538, 12.30769231]
q1 = [0.58466256, 7.162116362, 1.983965751, 11.83820886]

@testset "PriceIndexes.jl" begin
    @test PriceIndexes.BilateralIndexFormulas.carli(p1, p0) ≈ 0.85625 atol = 1e-5
    @test PriceIndexes.jevons(p1, p0) ≈ 0.822287308 atol = 1e-5
    @test PriceIndexes.laspeyres(p1, p0, q0) ≈ 0.967307692 atol = 1e-5

end

