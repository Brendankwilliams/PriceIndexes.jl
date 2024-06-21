using PriceIndexes
using Test

p0 = [2.0, 1.0, 1.0, 0.5]
p1 = [1.75, 0.5, 0.95, 0.55]
q0 = [0.384615385, 1.538461538, 1.538461538, 12.30769231]
q1 = [0.58466256, 7.162116362, 1.983965751, 11.83820886]
σ = 0.7

@testset "Bilateral Formulas" begin
    @test PriceIndexes.bmw(p1, p0) ≈ 0.82147799421931 atol=1e-5
    @test PriceIndexes.carli(p1, p0) ≈ 0.85625 atol=1e-5
    @test PriceIndexes.cswd(p1, p0) ≈ 0.819125217607815 atol=1e-5
    @test PriceIndexes.dutot(p1, p0) ≈ 0.833333333333333 atol=1e-5
    @test PriceIndexes.jevons(p1, p0) ≈ 0.822287307949555 atol=1e-5
    @test PriceIndexes.harmonic(p1, p0) ≈ 0.783610069630423 atol=1e-5
    @test PriceIndexes.banerjee(p1, p0, q1, q0) ≈ 0.874822206457738 atol=1e-5
  # @test PriceIndexes.bialek() ≈ 0.880989388595227 atol=1e-5
    @test PriceIndexes.davies(p1, p0, q1, q0) ≈ 0.887379354876721 atol=1e-5
    @test PriceIndexes.drobisch(p1, p0, q1, q0) ≈ 0.884035461190647 atol=1e-5
    @test PriceIndexes.fisher(p1, p0, q1, q0) ≈ 0.880104784765626 atol=1e-5
    @test PriceIndexes.geary_khamis(p1, p0, q1, q0) ≈ 0.92290239682701 atol=1e-5
    @test PriceIndexes.laspeyres(p1, p0, q0) ≈ 0.967307692307692 atol=1e-5
    @test PriceIndexes.lehr(p1, p0, q1, q0) ≈ 0.925890032785708 atol=1e-5
    @test PriceIndexes.marshall_edgeworth(p1, p0, q1, q0) ≈ 0.864246196464236 atol=1e-5
    @test PriceIndexes.paasche(p1, p0, q1) ≈ 0.800763230073602 atol=1e-5
    @test PriceIndexes.palgrave(p1, p0, q1) ≈ 0.895264545283374 atol=1e-5
    @test PriceIndexes.sato_vartia(p1, p0, q1, q0) ≈ 0.895264545283374 atol=1e-5
    @test PriceIndexes.stuvel(p1, p0, q1, q0) ≈ 0.858364550544791 atol=1e-5
    @test PriceIndexes.tornqvist(p1, p0, q1, q0) ≈ 0.892571497899412 atol=1e-5
    @test PriceIndexes.montgomery_vartia(p1, p0, q1, q0) ≈ 0.894193937668097 atol=1e-5
    @test PriceIndexes.walsh(p1, p0, q1, q0) ≈ 0.895264545283374 atol=1e-5
    @test PriceIndexes.lloyd_moulton(p1, p0, q0, σ) ≈ 0.9463576 atol=1e-5

end

