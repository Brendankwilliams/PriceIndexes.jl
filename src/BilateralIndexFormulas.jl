module BilateralIndexFormulas
using DataFrames, Statistics, StatsBase

export laspeyres, paasche, fisher, carli, dutot, cswd, bmw, drobish, tornqvist, jevons

"""
    laspeyres(p1::Vector{float64}, p0::Vector{float64}, q0::Vector{Float64} )
Laspeyres function with only vectors as inputs
"""
function laspeyres(p1::Vector{Float64}, p0::Vector{Float64}, q0::Vector{Float64})
    return sum(p1 .* q0) / (sum(p0 .* q0))
end

"""
    paasche(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64} )
    paasche function with only vectors as inputs
"""
function paasche(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64})
    return sum(p1 .* q1) / (sum(p0 .* q1))
end

"""
    fisher(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64} )
    Fisher price index
"""
function fisher(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    return sqrt(laspeyres(p1, p0, q0) * paasche(p1, p0, q1))
end


@doc raw"""
    carli(p1::Vector{Float64}, p0::Vector{Float64})
Carli price index, the arithmetic mean of price relatives
``P_{Carli}(p^{1},p^{0}) =
\frac{1}{N}\sum_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)``
"""
function carli(p1::Vector{Float64}, p0::Vector{Float64})
    return mean((p1 ./ p0))
end

"""
    dutot(p1::Vector{Float64}, p0::Vector{Float64})
    ``P_{Dutot}(p^{1},p^{0}) = ``
The Dutot index, the ratio of the arithmetic means of prices in each period
"""
function dutot(p1::Vector{Float64}, p0::Vector{Float64})
    return mean(p1) / mean(p0)
end
"""
    jevons(p1::Vector{Float64}, p0::Vector{Float64})


"""
function jevons(p1::Vector{Float64}, p0::Vector{Float64})
    return geomean((p1 ./ p0))
end
"""
    harmonic(p1::Vector{Float64}, p0::Vector{Float64})

Also known as Coggeshall, the Harmonic mean of price relatives using the StatsBase function

"""
function harmonic(p1::Vector{Float64}, p0::Vector{Float64})
    return harmmean((p1 ./ p0))
end
"""
    cswd(p1::Vector{Float64}, p0::Vector{Float64})

The Carruthers, Sellwood, Ward, Dalén index.
     A geometric average of the Carli and the harmonic mean of price relatives index.
"""
function cswd(p1::Vector{Float64}, p0::Vector{Float64})
    return sqrt(harmonic(p1, p0) * carli(p1, p0))
end

"""
    bmw(p1::Vector{Float64}, p0::Vector{Float64})
    Balk-Mehrhoff-Walsh (BMW) index

"""
function bmw(p1::Vector{Float64}, p0::Vector{Float64})
    numer = sum((p1 ./ p0) .* sqrt.(p0 ./ p1))
    denom = sum(sqrt.(p0 ./ p1))
    return numer / denom
end
"""
    drobish(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64} )


"""
function drobish(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    return (laspeyres(p1, p0, q0) + paasche(p1, p0, q1)) / 2.0
end

function tornqvist(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    total_expend0 = sum(expend0)
    total_expend1 = sum(expend1)

    avg_share = 0.5 * ((expend1 ./total_expend1) + (expend0 ./ total_expend0))

    return exp(sum(log.( p1 ./ p0 ) .* avg_share))
end

function walsh(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    avg_q = sqrt.(q1 .* q0)
    return sum(avg_q .* p1) / sum(avg_q .* p0)
end

"""
    geary_khamis(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
"""
function geary_khamis(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    #Calculate the harmonic mean of q for each element of q1 and the corresponding q0
    harm_q = 1 ./ (1 ./ q1 + 1 ./ q0)
    return sum(harm_q .* p1) / sum(harm_q .* p0)
end
"""
    palgrave(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64})

"""
function palgrave(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64})
    expend1 = p1 .* q1
    total_expend1 = sum(expend1)

    expend_share = (expend1 ./ total_expend1)

    return sum(expend_share .* p1 ./ p0) 
end

"""
    marshall_edgeworth(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
"""
function marshall_edgeworth(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    return sum((q0+q1) .* p1) / sum((q0+q1) .* p0)
end

"""
    sato_vartia(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})

"""
function sato_vartia(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    total_expend0 = sum(expend0)
    total_expend1 = sum(expend1)

    weight_share1 = expend1 ./ total_expend1
    weight_share0 = expend0 ./ total_expend0
    #sv_weight_share = (weight_share1-weight_share0) ./ (log.(weight_share1)-log.(weight_share0))
    sv_weight_share = vartia_L(weight_share1, weight_share0)
    #Normalize weight share
    sv_weight_share_norm = sv_weight_share/sum(sv_weight_share)

    return exp(sum(log.( p1 ./ p0 ) .* sv_weight_share_norm))
end
"""
    vartia_L(a::Vector{Float64}, b::Vector{Float64})

    Helper function related to Sato-Vartia and montgomery-vartia
"""
function vartia_L(a::Float64, b::Float64)
    if b == a 
        return b
    else
        return (a-b) ./ (log.(a)-log.(b))
    end
end

"""
    montgomery_vartia(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})

    Montgomery-Vartia or sometimes just Vartia or Vartia-I, but preceeded by Montgomery (1937).
    Implemented based on https://doi.org/10.1111/rssa.12633 von Auer and Wengenroth (2020).
"""
function montgomery_vartia(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    denom = vartia_L(sum(expend1), sum(expend0))
    
    return exp(sum(log.(p1 ./ p0).*vartia_L.(expend1, expend0)/denom))
end

"""
    lloyd_moulton(p1::Vector{Float64}, p0::Vector{Float64}, q0::Vector{Float64}, σ::Float64)

"""
function lloyd_moulton(p1::Vector{Float64}, p0::Vector{Float64}, q0::Vector{Float64}, σ::Float64)
    expend0 = p0 .* q0
    total_expend0 = sum(expend0)
    weight_share0 = expend0 ./ total_expend0
    return sum(weight_share0  .* ((p1 ./ p0) .^ (1 - σ))) ^ (1 /(1-σ))
end

"""
    stuvel(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})

"""
function stuvel(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
   # lasp_Q = sum(p0 .* q1) / sum(p0 .* q0)
   # lasp_P = sum(p1 .* q0) / sum(p0 .* q0) 
    #The term A is the difference between the laspeyres price and quuantity indexes, divided by 2.0
    A = (sum(p1 .* q0) / sum(p0 .* q0) - sum(p0 .* q1) / sum(p0 .* q0) )/2.0
    v1 = sum(p1 .* q1)
    v0 = sum(p0 .* q0)
    return A + sqrt(A^2.0 + v1/v0)
end

"""
    lehr(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})

TBW
"""
function lehr(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    numer = sum( (expend0 + expend1) .* q0 ./ (q0+q1) )
    denom = sum( (expend0 + expend1) .* q1 ./ (q0+q1) )

    return ( sum(expend1)/denom * numer/sum(expend0) )
end
"""
    banerjee(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})

TBW
"""
function banerjee(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    return (sum(expend1)/sum(expend0)) * sum(expend0 .* (1.0 .+ p1 ./ p0))/sum(expend1 .* (1.0 .+ p0 ./ p1))
end

function davies(p1::Vector{Float64}, p0::Vector{Float64}, q1::Vector{Float64}, q0::Vector{Float64})
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    return (sum(expend1)/sum(expend0)) * sum(expend0 .* sqrt.(p1 ./ p0))/sum(expend1 .* sqrt.(p0 ./ p1))
end

function raito_of_harmonic_means()
    
end

function białek()
end

end