module MultilateralIndexFormulas
using ..BilateralIndexFormulas # Import the functions from the BilateralIndexFormulas module
using DataFrames
using LinearAlgebra

export geary_khamis, GEKS, GEKS_general, CCDI

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
    time_num, prod_num = size(price_matrix) # Matrices are time period x product
    exp_matrix = price_matrix .* quantity_matrix
    exp_share = exp_matrix ./ sum(exp_matrix, dims=2) #divide expenditure by row-level (time period total)

    prod_sums = sum(quantity_matrix, dims=1) #Sum by column to get total product counts across time
    D = inv(diagm(vec(prod_sums))) #Create a Diagonal Matrix with the product level sums and take the inverse

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

function CCDI(price_matrix::Matrix{Float64}, quantity_matrix::Matrix{Float64})
    time_num, prod_num = size(price_matrix) 

    rel_mat = ones(time_num, time_num)
    
    for t in 1:time_num
        for a in (t+1):time_num #Only fill lower triangle
            rel_mat[a, t] = tornqvist(price_matrix[t, :], price_matrix[a, :], 
            quantity_matrix[t, :], quantity_matrix[a, :])
        end
    end
    #Tranpose and invert elements to fill upper triangle
    bilat_matrix = transpose((rel_mat).^(-1)) .* rel_mat 

    #Alternative geometric mean calculation exp.(cumsum(log.(bilat_matrix), dims =1)./12)
    results = cumprod(bilat_matrix, dims=1)[12, :].^(1/12) #Get geometric average
    return results/results[1] #Normalize to first period and return results
end

function GEKS(price_matrix::Matrix{Float64}, quantity_matrix::Matrix{Float64})
    time_num, prod_num = size(price_matrix) 

    rel_mat = ones(time_num, time_num)
    
    for t in 1:time_num
        for a in (t+1):time_num #Only fill lower triangle
            rel_mat[a, t] = fisher(price_matrix[t, :], price_matrix[a, :], 
            quantity_matrix[t, :], quantity_matrix[a, :])
        end
    end
    #Tranpose and invert elements to fill upper triangle
    bilat_matrix = transpose((rel_mat).^(-1)) .* rel_mat 

    #Alternative geometric mean calculation exp.(cumsum(log.(bilat_matrix), dims =1)./12)
    results = cumprod(bilat_matrix, dims=1)[12, :].^(1/12) #Get geometric average
    return results/results[1] #Normalize to first period and return results
end

"""
    GEKS_general(price_matrix::Matrix{Float64}, quantity_matrix::Matrix{Float64}, index_formula::Function)

    GEKS type multilateral calculation that takes a bilateral index_formula as an argument. When Tornqvist is
    used this is equivalent to CCDI. When Fisher is used this is equivalent to GEKS. Only bilateral Formulas
    with four arguments (price1, price0, quantity1, quantity0) are supported.
"""
function GEKS_general(price_matrix::Matrix{Float64}, quantity_matrix::Matrix{Float64}, index_formula::Function)
    time_num, prod_num = size(price_matrix) 

    rel_mat = ones(time_num, time_num)
    
    for t in 1:time_num
        for a in (t+1):time_num #Only fill lower triangle
            rel_mat[a, t] = index_formula(price_matrix[t, :], price_matrix[a, :], quantity_matrix[t, :], 
            quantity_matrix[a, :])
        end
    end
    #Tranpose and invert elements to fill upper triangle
    bilat_matrix = transpose((rel_mat).^(-1)) .* rel_mat 

    #Alternative geometric mean calculation exp.(cumsum(log.(bilat_matrix), dims =1)./12)
    results = cumprod(bilat_matrix, dims=1)[12, :].^(1/12) #Get geometric average
    return results/results[1] #Normalize to first period and return results
end

end #module MultilateralIndexFormulas