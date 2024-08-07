using PriceIndexes
using Test
using Dates
using DataFrames

p0 = [2.0, 1.0, 1.0, 0.5]
p1 = [1.75, 0.5, 0.95, 0.55]
p_end = [1.5, 1.1, 0.85, 0.65] #end period where p1 and q1 are used as base
q0 = [0.384615385, 1.538461538, 1.538461538, 12.30769231]
q1 = [0.58466256, 7.162116362, 1.983965751, 11.83820886]
q_end = [0.9149417, 1.7013378, 2.8492993, 9.7449408] #end period where p1 and q1 are used as base
σ = 0.7 #elasticity for Lloyd-Moulton 

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
    @test PriceIndexes.geary_khamis_bi(p1, p0, q1, q0) ≈ 0.92290239682701 atol=1e-5
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
    @test PriceIndexes.lloyd_moulton(p1, p0, q1, q0; σ = 0.7) ≈ 0.9463576 atol=1e-5
    @test PriceIndexes.ag_mean(p1, p0, q0, 0.7) ≈ 0.9453866 atol=1e-5
    @test PriceIndexes.geolaspeyres(p1, p0, q0) ≈ 0.9359918 atol=1e-5
    @test PriceIndexes.geopaasche(p1, p0, q1) ≈ 0.8511654 atol=1e-5
    # base quantity formula tests
    @test PriceIndexes.lowe(p_end, p0, q1) ≈ 1.117159 atol=1e-5
    #@test PriceIndexes.geolowe(p_end, p0, q0) ≈ 1.101997 atol=1e-5
    #@test PriceIndexes.geoyoung(p_end, p0, p0, q0) ≈ 1.117903 atol=1e-5
    @test PriceIndexes.young(p_end, p0; p_b = p1, q_b = q1) ≈ 1.136377 atol=1e-5


end

price_matrix = [
    2.0   1.0   1.0   0.5;
    1.75  0.5   0.95  0.55;
    1.6   1.05  0.9   0.6;
    1.5   1.1   0.85  0.65;
    1.45  1.12  0.4   0.7;
    1.4   1.15  0.8   0.75;
    1.35  1.18  0.75  0.7;
    1.3   0.6   0.72  0.65;
    1.25  1.2   0.7   0.7;
    1.2   1.25  0.4   0.75;
    1.15  1.28  0.7   0.75;
    1.1   1.3   0.65  0.8
]

quantity_matrix = [
    0.384615  1.53846   1.53846  12.3077;
    0.584663  7.16212   1.98397  11.8382;
    0.71355   1.65686   2.25517  10.1483;
    0.914942  1.70134   2.8493    9.74494;
    1.02806   1.72313  13.5093    8.82241;
    1.20582   1.78708   3.69283   8.40325;
    1.32933   1.73995   4.30702   9.88858;
    1.45749   6.84211   4.75146  11.6599;
    1.62188   1.75986   5.17182  10.3436;
    1.83824   1.69412  16.5441    9.41176;
    2.1055    1.69954   5.68269   9.90052;
    2.4576    1.75959   7.03834   9.29281
]

@testset "MultilateralIndexFormulas Formulas" begin
    @test PriceIndexes.geary_khamis(price_matrix, quantity_matrix)[12] ≈ 1.2324531452384582 atol=1e-5
    @test PriceIndexes.CCDI(price_matrix, quantity_matrix)[12] ≈ 1.1328205 atol=1e-5
    @test PriceIndexes.GEKS(price_matrix, quantity_matrix)[12] ≈ 1.1152731 atol=1e-5
    @test PriceIndexes.GEKS_general(price_matrix, quantity_matrix, PriceIndexes.walsh)[12] ≈ 1.1370031 atol=1e-5
    @test PriceIndexes.similarity_linking(price_matrix, quantity_matrix, 
        PriceIndexes.predicted_share_diff, PriceIndexes.laspeyres)[12] ≈  1.18350857 atol=1e-5
    @test PriceIndexes.similarity_linking(price_matrix, quantity_matrix, 
        PriceIndexes.paasche_laspeyres_spread, PriceIndexes.tornqvist)[12] ≈  1.1366808 atol=1e-5
    @test PriceIndexes.spq(price_matrix, quantity_matrix, PriceIndexes.fisher)[12] ≈  1.135090 atol=1e-5
end # @testset "MultilateralIndexFormulas Formulas"

@testset "Date Calculations" begin
    @test PriceIndexes.date_to_period(Date(2024, 7, 30), Date(1993, 10, 11), "days") == 11250
    @test PriceIndexes.period_to_date(11250, Date(1993, 10, 11), "days") == Date(2024, 7, 30)
end

# Define the columns as vectors
time = ["1/1/2015", "2/1/2015", "3/1/2015", "4/1/2015", "5/1/2015", "6/1/2015", "7/1/2015", "8/1/2015", "9/1/2015", "10/1/2015", "11/1/2015", "12/1/2015",
        "1/1/2015", "2/1/2015", "3/1/2015", "4/1/2015", "5/1/2015", "6/1/2015", "7/1/2015", "8/1/2015", "9/1/2015", "10/1/2015", "11/1/2015", "12/1/2015",
        "1/1/2015", "2/1/2015", "3/1/2015", "4/1/2015", "5/1/2015", "6/1/2015", "7/1/2015", "8/1/2015", "9/1/2015", "10/1/2015", "11/1/2015", "12/1/2015",
        "1/1/2015", "2/1/2015", "3/1/2015", "4/1/2015", "5/1/2015", "6/1/2015", "7/1/2015", "8/1/2015", "9/1/2015", "10/1/2015", "11/1/2015", "12/1/2015"]

prices = [2, 1.75, 1.6, 1.5, 1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1,
          1, 0.5, 1.05, 1.1, 1.12, 1.15, 1.18, 0.6, 1.2, 1.25, 1.28, 1.3,
          1, 0.95, 0.9, 0.85, 0.4, 0.8, 0.75, 0.72, 0.7, 0.4, 0.7, 0.65,
          0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.7, 0.65, 0.7, 0.75, 0.75, 0.8]

quantities = [0.384615385, 0.58466256, 0.713550247, 0.914941664, 1.0280574, 1.205823442, 1.329328335, 1.457489879, 1.621882542, 1.838235294, 2.105497116, 2.457603559,
              1.538461538, 7.162116362, 1.656860438, 1.701337804, 1.723127139, 1.787080488, 1.739946058, 6.842105263, 1.759855189, 1.694117647, 1.699536094, 1.75958598,
              1.538461538, 1.983965751, 2.255171151, 2.849299298, 13.50931677, 3.69283429, 4.307023806, 4.751461988, 5.171819331, 16.54411765, 5.682693746, 7.038343919,
              12.30769231, 11.83820886, 10.14827018, 9.744940795, 8.822410952, 8.403249584, 9.888575064, 11.65991903, 10.34363866, 9.411764706, 9.900515327, 9.292813456]

prodID = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
          3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

# Combine into a DataFrame
cesDF = DataFrame(time = time, prices = prices, quantities = quantities, prodID = prodID)

@testset "PriceRelativeFrame Operations" begin
    prf = PriceIndexes.PriceFrame(cesDF, Date(2000,1,1),"months",1, Symbol[], [:prodID], :prices, :quantities)
    @test PriceIndexes.direct_index(prf, PriceIndexes.walsh).price_df.index_relative[12] ≈ 113.70031 atol=1e-5
    @test PriceIndexes.direct_index(prf, PriceIndexes.jevons).price_df.index_relative[12] ≈ 92.86132 atol=1e-5
end # @testset "PriceRelativeFrame Operations"