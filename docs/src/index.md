# PriceIndexes.jl

Documentation for PriceIndexes.jl

## PriceFrames
The PriceFrames module contains the basic structs and functions for users to set up data and generate price indexes.

```@autodocs
    Modules = [PriceIndexes.PriceFrames]
```
## Bilateral Index Formulas
Functions to apply price index formulas to price and quantity data.

Most functions take four arguments: current and previous prices and quantities. Not all price index formulas use all four arguments, but each has four argument call to enable standard calls from general price index functions. Multiple dispatch is used to distinguish, which function call to use.

```@autodocs
    Modules = [PriceIndexes.BilateralIndexFormulas]
```
## Multilateral Index Formulas
Multilateral index formulas make comparisons across multiple periods. The methods are most commonly used to address chain drift.

```@autodocs
    Modules = [PriceIndexes.MultilateralIndexFormulas]
```

