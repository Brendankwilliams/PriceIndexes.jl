module MultilateralIndexFormulas
using LinearAlgebra, DataFrames

export geary_khamis

function geary_khamis(df::DataFrame)
    # Unstack the DataFrame to create a matrix for the values of x
    df_wide = unstack(df, :time, :prodID, :prices)
    df_wideQ = unstack(df, :time, :prodID, :quantities)

    # Extract the matrix from the unstacked DataFrame
    # Matrices are time period x product
    price_matrix = Matrix(df_wide[:, Not(:time)])
    quantity_matrix = Matrix(df_wideQ[:, Not(:time)])
    exp_matrix = price_matrix .* quantity_matrix
    exp_share = exp_matrix ./ sum(exp_matrix, dims=2) #divide expenditure by row-level (time period total)

    prod_sums = sum(quantity_matrix, dims=1) #Sum by column to get total product counts across time

    D = inv(diagm(vec(prod_sums))) #Create a Diagonal Matrix with the product level sums and take the inverse

    prod_num = length(prod_sums)  #Get number of products
    time_num = size(exp_share)[1] #Get number of time periods
    sq = zeros(prod_num, prod_num)
    for i in 1:time_num
        sq = sq + quantity_matrix[i, :] * transpose(exp_share[i, :])
    end

    C = D * transpose(sq) #Note indexNumr does not use tranpose, but using it here matches it

    R = zeros(prod_num, prod_num)
    R[1, :] .= 1
    b = inv(I-C+R)[:, 1]

    p_index = Vector{Float64}(undef, time_num)
    for i in 1:time_num
        tempDF = filter(row -> row.time == i, df)
        p_index[i] = sum(tempDF[:,:prices] .* tempDF[:, :quantities])/ sum(b .* tempDF[:, :quantities])
    end

    return p_index/p_index[1] #Normalize all periods by period 1, first period to 1
end #function geary_khamis

function geary_khamis(price_matrix::Matrix{Float64}, quantity_matrix::Matrix{Float64})
    # Matrices are time period x product
    exp_matrix = price_matrix .* quantity_matrix
    exp_share = exp_matrix ./ sum(exp_matrix, dims=2) #divide expenditure by row-level (time period total)

    prod_sums = sum(quantity_matrix, dims=1) #Sum by column to get total product counts across time

    D = inv(diagm(vec(prod_sums))) #Create a Diagonal Matrix with the product level sums and take the inverse

    prod_num = length(prod_sums)  #Get number of products
    time_num = size(exp_share)[1] #Get number of time periods
    sq = zeros(prod_num, prod_num)
    for i in 1:time_num
        sq = sq + quantity_matrix[i, :] * transpose(exp_share[i, :])
    end

    C = D * transpose(sq) #Note indexNumr does not use tranpose, but using it here matches it

    R = zeros(prod_num, prod_num)
    R[1, :] .= 1
    b = inv(I-C+R)[:, 1]

    p_index = Vector{Float64}(undef, time_num)
    for i in 1:time_num
        p_index[i] = sum(price_matrix[i, :] .* quantity_matrix[i, :])/ sum(b .* quantity_matrix[i, :])
    end

    return p_index/p_index[1] #Normalize all periods by period 1, first period to 1
end #function geary_khamis

end #module MultilateralIndexFormulas