module PriceRelativeFrame
using DataFrames, Dates

export date_to_period, period_to_date, merge_periods, merge_one_period

#Module for creating a PriceRelativeFrame struct and related functions and structures




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
    date_to_period(date::Date, epoch::Date, time_unit::String)

    Take a date and a time unit and calculate the number of periods since epoch and create a time
    period count that can faciltate price index calculations without having to manipulate dates.
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
        throw(ArgumentError("Invalid period specified. Choose from days, weeks, months, quarters, or years"))
    end
end

function period_to_date(period_count::Integer, epoch::Date, time_unit::String)
    if time_unit == "days" 
        return epoch + Day(period_count)
    elseif time_unit == "weeks"
        return epoch + Week(period_count)
    elseif time_unit == "months"
        return epoch + Month(period_count)
    elseif time_unit == "years"
        return epoch + Year(period_count)
    else
        error("Unsupported unit. Use day, week, month, or year.")
    end
end





end #end module PriceRelativeFrame