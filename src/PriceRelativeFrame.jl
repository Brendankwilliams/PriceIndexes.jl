module PriceRelativeFrame
using DataFrames, Dates

export date_to_period, period_to_date, merge_periods, merge_one_period, merge_direct_periods, IndexFrame, direct_index,
    PriceRelativeFrame

#Module for creating a PriceRelativeFrame struct and related functions and structures
abstract type AbstractPriceFrame end

#merge data frames with dynamic byvals
#Use in pricerelativeFrame constuctor
#function sort_dataframe(df::DataFrame, sort_by::Union{Symbol, Vector{Symbol}}=:ColName)
mutable struct PriceRelativeFrame <: AbstractPriceFrame
    price_df::DataFrame
    epoch_date::Date
    time_unit::String
    time_frequency::Integer
    output_variables::Vector{Symbol}
    product_definition::Vector{Symbol}
    price_variable::Symbol
    quantity_variable::Symbol


    function PriceRelativeFrame(df::DataFrame, epoch_date::Date, time_unit::String, time_frequency::Integer, 
        product_definition::Vector{Symbol}, output_variables::Vector{Symbol}, 
        price_variable::Symbol, quantity_variable::Symbol)
        
        df.time = Date.(df.time, "m/d/yyyy")
        #Convert to time index
        df.time_count = calculate_periods.(df.time, epoch_date, time_unit)
        
        #Merge dataframe by frequency, outputlevel, and product definition, and then group by output level

        #calculate price relative and expenditure
        new(df, epoch_date, time_unit, time_frequency, 
            output_variables, product_definition, price_variable, quantity_variable)
    end
end

function PriceRelativeFrame(df::DataFrame; epoch_date::Date=Date(1978, 1, 1), time_unit::String="Month",
    time_frequency::Integer=1, product_definition::Vector{Symbol}=:product_id, output_variables::Vector{Symbol}=Symbol[], 
    price_variable::Symbol=:price, quantity_variable::Symbol=:quantity))
    return PriceRelativeFrame(df, epoch_date, time_unit, time_frequency, product_definition, 
        output_variables, price_variable, quantity_variable)
end

#Create IndexFrame as a type of PriceRelativeFrame. Should inherit functions from PRF but allow new functions
#e.g. to convert relatives 
mutable struct IndexFrame <: AbstractPriceFrame 
    price_df::DataFrame #Need same name as PriceRelativeFZrame
    epoch_date::Date
    time_unit::String
    time_frequency::Integer
    output_variables::Vector{Symbol}
    price_variable::Symbol
    quantity_variable::Symbol
    relative_variable::Symbol #Added variable compared to PriceRelativeFrame
    index_type::String #chained, direct, multilateral?
    base_date::Date #Base date for index to be set at 100
    
    function IndexFrame(index_df::DataFrame, epoch_date::Date, time_unit::String, time_frequency::Integer, output_variables::Vector{Symbol},
            price_variable::Symbol, quantity_variable::Symbol, relative_variable::Symbol, index_type::String, base_date::Date)
        if index_type ∉ ["chained", "direct", "multilateral"]
            error("Invalid identifier. Must be one of \"chained\", \"direct\", or \"multilateral\".")
        end
        
        #Convert time_count back to date format
        index_df[!, :Date] = period_to_date.(index_df[!, :time_count], epoch_date, time_unit)
        
        new(index_df, epoch_date, time_unit, time_frequency, output_variables, price_variable, 
            quantity_variable, relative_variable, index_type, base_date)
    end
end

#Constructor that inherits most values from a parent price relative frame
#Use a January 1, 1978 epoch that starts on a Sunday to set a week definition with weeks starting on Sunday
#NEED TO FIX EXPENDITURE
function IndexFrame(index_df::DataFrame, parent_df::AbstractPriceFrame; index_type::String, base_date::Date)
    return IndexFrame(index_df, parent_df.epoch_date, parent_df.time_unit, parent_df.time_frequency, parent_df.output_variables,
        :index_level, :expenditure, :index_relative, index_type, base_date)
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
    merge_direct_periods(df1::DataFrame, productDefCols::Vector{Symbol}; 
    outputLevelCols::Vector{Symbol}=Symbol[])

    Merge the first period in a dataframe with all other rows to prepare for the calcualtion of a direct price index.
"""
function merge_direct_periods(df1::DataFrame, productDefCols::Vector{Symbol}; 
    outputLevelCols::Vector{Symbol}=Symbol[])

    base_df = df1[df1.time .== minimum(df1.time), :]

    # [;] combines two vectors of symbols
    joinCols = [productDefCols; outputLevelCols]

    comboDF = innerjoin(base_df, df1, on=joinCols, makeunique=true, renamecols="_0" => "_1")
    rename!(comboDF, :time_count_1 => :time_count)
    return comboDF
end

"""
    date_to_period(date::Date, epoch::Date, time_unit::String)

    Take a date and a time unit and calculate the number of periods since epoch and create a time
    period count that can faciltate price index calculations without having to manipulate dates. For periods
    greater than "day," the first day of the period should represent the whole period. The epoch date 
    establishes beginning of the time unit. For example, an epoch date on a Monday establishes a week
    definition that starts on Monday, and an epoch date in September establishes a year that starts in September
"""
function date_to_period(date::Date, epoch::Date, time_unit::String)
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
        throw(ArgumentError("Invalid period specified. Choose from days, weeks, months, quarters, or years."))
    end
end

"""
    period_to_date(period_count::Integer, epoch::Date, time_unit::String)

    Convert a period count back to a date given an epoch date and time unit. 
    After price index calculations have been completed, this function will convert calculation 
    friendly period counts to user friendly dates.
"""
function period_to_date(period_count::Integer, epoch::Date, time_unit::String)
    if time_unit == "days" 
        return epoch + Day(period_count)
    elseif time_unit == "weeks"
        return epoch + Week(period_count)
    elseif time_unit == "months"
        return epoch + Month(period_count)
    elseif time_unit == "quarters"
        return epoch + Month(3*period_count)
    elseif time_unit == "years"
        return epoch + Year(period_count)
    else
        throw(ArgumentError("Invalid period specified. Choose from days, weeks, months, quarters, or years."))
    end
end

function direct_index(prf::AbstractPriceFrame, index_formula::Function; base_value::Real = 100, kwargs...)
    merge_df = merge_direct_periods(prf.price_df, prf.product_definition)
    #Group by for relative output level
    gdf = groupby(merge_df, vcat(prf.output_variables, :time_count))
    ans = combine(gdf,   [Symbol(prf.price_variable,"_1"), Symbol(prf.price_variable,"_0"),
                            Symbol(prf.quantity_variable,"_1"), Symbol(prf.quantity_variable,"_0") ] => 
                            ( (p1, p0, q1, q0) -> index_formula(p1, p0, q1, q0; kwargs...) ) => :index_relative)
    ans[!, :index_relative] *= base_value
    base_date = period_to_date(minimum(prf.price_df.time_count), prf.epoch_date, prf.time_unit) 
    IFdf = IndexFrame(ans, prf; index_type = "direct", base_date = base_date)
    return IFdf
end


end #end module PriceRelativeFrame