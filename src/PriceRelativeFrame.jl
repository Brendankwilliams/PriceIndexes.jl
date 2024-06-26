module PriceRelativeFrame
using DataFrames, Dates

#Module for creating a PriceRelativeFrame struct and related functions and structures

#merge data frames with dynamic byvals
#Use in pricerelativeFrame constuctor
#function sort_dataframe(df::DataFrame, sort_by::Union{Symbol, Vector{Symbol}}=:ColName)
struct PriceRelativeFrame 
    gdf::GroupedDataFrame
    time_unit::String
    time_frequency::Int32
    output_variables::Vector{Symbol}
    product_definition::Vector{Symbol}
    price_variable::Symbol
    quantity_variable::Symbol

    function PriceRelativeFrame(df::DataFrame, time_unit::String, time_frequency::Int32, 
        product_definition, output_variables::Vector{Symbol}, price_variable, quantity_variable)
        #Convert to time index
        #Merge dataframe by frequency, outputlevel, and product definition, and then group by output level
        gdf = groupby(df, output_variables)
        #calculate price relative and expenditure
        new(gdf, time_unit, time_frequency, output_variables, product_definition, price_variable, quantity_variable)
    end
end

function mergePricePeriods(df1::DataFrame, productDefCols::Union{Symbol,Vector{Symbol}}=:prodID, outputLevelCols::Union{Symbol,Vector{Symbol}}=:time)
    df2 = deepcopy(df1) #Use a deep copy of the data and lag it in order to merge
    df2.time = df2.time .+ 1

    # [;] combines two vectors of symbols
    joinCols = [productDefCols; outputLevelCols]

    comboDF = innerjoin(df1, df2, on=joinCols, makeunique=true, renamecols="_0" => "_1")

    comboDF.prcRelative = comboDF.prices_1 ./ comboDF.prices_0
    comboDF.expend1 = comboDF.prices_1 .* comboDF.quantities_1
    comboDF.expend0 = comboDF.prices_0 .* comboDF.quantities_0
    return comboDF
end


# Function to calculate periods since epoch
function time_periods_since_epoch(df::DataFrame, date_column::Symbol, epoch::Date, period::Symbol)

    # Calculate the number of periods since the epoch
    periods_since_epoch = [calculate_periods(row[date_column], epoch, period) for row in eachrow(df)]

    # Return a new DataFrame with the calculated periods
    return hcat(df, DataFrame(PeriodsSinceEpoch = periods_since_epoch))
end

# Helper function to calculate periods
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



end