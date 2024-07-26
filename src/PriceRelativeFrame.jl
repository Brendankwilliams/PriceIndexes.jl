module PriceRelativeFrame
using DataFrames, Dates

#Module for creating a PriceRelativeFrame struct and related functions and structures


#Use in pricerelativeFrame constuctor
#function sort_dataframe(df::DataFrame, sort_by::Union{Symbol, Vector{Symbol}}=:ColName)
mutable struct PriceRelativeFrame 
    price_df::DataFrame
    epoch_date::Date
    time_unit::String
    time_frequency::Integer
    output_variables::Vector{Symbol}
    product_definition::Vector{Symbol}
    price_variable::Symbol
    quantity_variable::Symbol

    function PriceRelativeFrame(df::DataFrame, epoch_date::Date, time_unit::String, time_frequency::Int64, 
        product_definition::Vector{Symbol}, output_variables::Vector{Symbol}, 
        price_variable::Symbol, quantity_variable::Symbol)
        
        df.time = Date.(df.time, "m/d/yyyy")
        #Convert to time index
        df.time_count = calculate_periods.(df.time, epoch_date, time_unit)
        
        #calculate price relative and expenditure
        new(df, epoch_date, time_unit, time_frequency, 
            output_variables, product_definition, price_variable, quantity_variable)
    end
end

"""
    merge_periods(df1::DataFrame, time_frequency::Integer,
    productDefCols::Vector{Symbol}, outputLevelCols::Vector{Symbol})

    Merge dataframe with itself after incrementing the time index by frequency. The resulting dataframe
    has price and weight information a single row to set up price index calclations. This version works across
    all time periods in the dataframe and excludes unmatched prices. See merge_one_period for a version 
    that works incrementaly to accomodate imputation.
    Used in DataFrame portion of a PriceRelativeFrame.
"""
function merge_periods(df1::DataFrame, time_frequency::Integer,
    productDefCols::Vector{Symbol}, outputLevelCols::Vector{Symbol})
    df1.expend = df1.prices .* df1.quantities
    df2 = deepcopy(df1) #Use a deep copy of the data and lag it in order to merge
    df2.time_count  = df2.time_count .+ time_frequency

    # merge data frames with dynamic byvals
    # [;] combines two vectors of symbols
    joinCols = [productDefCols; outputLevelCols]

    comboDF = innerjoin(df1, df2, on=joinCols, makeunique=true, renamecols="_0" => "_1")

    return comboDF
end

"""
    merge_one_period(df1::DataFrame, time_period::Integer, time_frequency::Integer,
    productDefCols::Vector{Symbol}, outputLevelCols::Vector{Symbol})

    Incremental version of merge_periods that is designed to be called with a loop. Facilitates price index
    calclations that use imputation.
"""
function merge_one_period(df1::DataFrame, time_period::Integer, time_frequency::Integer,
    productDefCols::Vector{Symbol}, outputLevelCols::Vector{Symbol})
    per1 =  filter(row -> row[:time_count] == time_period, df1) 
    per2 =  filter(row -> row[:time_count] == (time_period-time_frequency), df1)
    per2.time_count  = per2.time_count .+ time_frequency

    # [;] combines two vectors of symbols
    joinCols = [productDefCols; outputLevelCols]

    comboDF = innerjoin(per1, per2, on=joinCols, makeunique=true, renamecols="_0" => "_1")
    #Antijoins capture missing values, stored as separate dfs in struct. Can be ignored if no imputation
    #or appended and then used iteratively if imputation is called
    curr_miss_df = antijoin(per2, per1, on=joinCols, makeunique=true)
    return comboDF, curr_miss_df
end

"""
    calculate_periods(date::Date, epoch::Date, time_unit::String)

    Take a date and a time unit and calculate the number of periods since epoch and create a time
    period count that can faciltate price index calculations without having to manipulate dates.
"""
function calculate_periods(date::Date, epoch::Date, time_unit::String)
    if time_unit == "days" 
        return Dates.value(date - epoch)
    elseif time_unit == "weeks" #Defintion of week is determined by the day of the week of epoch
        return div(Dates.value(date - epoch), 7)
    elseif time_unit == "months"
        return (year(date) - year(epoch)) * 12 + (month(date) - month(epoch))
    elseif time_unit == "quarters" #Defintion of quarter is determined by the month of epoch
        return (year(date) - year(epoch)) * 4 + div(month(date)-1, 3) - div(month(epoch)-1, 3)
    elseif time_unit == "years"
        return div((year(date) - year(epoch)) * 12 + (month(date) - month(epoch)), 12)
    else
        throw(ArgumentError("Invalid period specified. Choose from days, weeks, months, quarters, or years"))
    end
end 


function price_index(prf::PriceRelativeFrame, index_formula::Function; imputation=nothing, kwargs...)
    #If no imputes, just apply formula time as by variable
    if imputation == nothing
        #Merge pricing periods and ignore prices that are unmatched to make relatives
        matchDF = merge_periods(prf.price_df, prf.time_frequency,
            prf.product_definition, [prf.output_variables; :time_count])
        #Outlier removal TO BE ADDED
        
        #Group by for relative output level
        gdf = groupby(matchDF, vcat(prf.output_variables, :time_count))
        ans = combine(gdf,   [Symbol(prf.price_variable,"_1"), Symbol(prf.price_variable,"_0"),
                            Symbol(prf.quantity_variable,"_1"), Symbol(prf.quantity_variable,"_0") ] => 
                            ( (p1, p0, q1, q0) -> index_formula(p1, p0, q1, q0; kwargs...) ) => :price_index)
    end
    
    #HANDLING OF IMPUTATION TO BE ADDED
    return ans
                
    
    #return a price relative frame where prices are index levels, expenditure calculate
end

function ksigma_outliers(prf::PriceRelativeFrame, groupby_col::Symbol = nothing)
    if isnothing(groupby_col)
        mean_val = mean(df[!, column])
        std_val = std(df[!, column])
        return df[abs.(df[!, column] .- mean_val) .<= x * std_val, :]
    else
        return combine(groupby(df, groupby_col)) do subdf
            mean_val = mean(subdf[!, column])
            std_val = std(subdf[!, column])
            return subdf[abs.(subdf[!, column] .- mean_val) .<= x * std_val, :]
        end
    end
end



end #end module PriceRelativeFrame