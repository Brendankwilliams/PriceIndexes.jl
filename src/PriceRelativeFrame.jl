module PriceRelativeFrame
using DataFrames, Dates

#Module for creating a PriceRelativeFrame struct and related functions and structures

#merge data frames with dynamic byvals
#Use in pricerelativeFrame constuctor
#function sort_dataframe(df::DataFrame, sort_by::Union{Symbol, Vector{Symbol}}=:ColName)
struct PriceRelativeFrame 
    rel_df::DataFrame
    curr_miss_df::DataFrame
    prior_miss_df::DataFrame
    epoch_date::Date
    time_unit::String
    time_frequency::Int32
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
        
        #Merge dataframe by frequency, outputlevel, and product definition, and then group by output level
        rel_df, curr_miss_df, prior_miss_df = mergePricePeriods(df, time_frequency, 
            product_definition, output_variables)
        #calculate price relative and expenditure
        new(rel_df, curr_miss_df, prior_miss_df, epoch_date, time_unit, time_frequency, 
            output_variables, product_definition, price_variable, quantity_variable)
    end
end

function mergePricePeriods(df1::DataFrame, time_frequency::Int64,
    productDefCols::Union{Symbol,Vector{Symbol}}=:prodID, outputLevelCols::Union{Symbol,Vector{Symbol}}=:time)
    df1.expend = df1.prices .* df1.quantities
    df2 = deepcopy(df1) #Use a deep copy of the data and lag it in order to merge
    df2.time_count  = df2.time_count .+ time_frequency

    # [;] combines two vectors of symbols
    joinCols = [productDefCols; outputLevelCols]

    comboDF = innerjoin(df1, df2, on=joinCols, makeunique=true, renamecols="_0" => "_1")
    #Antijoins capture missing values, stored as separate dfs in struct. Can be ignored if no imputation
    #or appended and then used iteratively if imputation is called
    curr_miss_df = antijoin(df1, df2, on=joinCols, makeunique=true)
    prior_miss_df = antijoin(df2, df1, on=joinCols, makeunique=true)

    comboDF.prcRelative = comboDF.prices_1 ./ comboDF.prices_0
    return comboDF, curr_miss_df, prior_miss_df
end

#  function to calculate periods
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



end #end module PriceRelativeFrame