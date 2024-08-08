module BilateralIndexFormulas
using DataFrames, Statistics, StatsBase

export laspeyres, paasche, fisher, carli, dutot, cswd, bmw, drobisch, tornqvist, jevons,
    sato_vartia, marshall_edgeworth, geary_khamis_bi, palgrave, harmonic, banerjee, bialek, davies, lehr, stuvel,
    walsh, lloyd_moulton, montgomery_vartia, lowe, geolowe, young, geoyoung, geolaspeyres, geopaasche, ag_mean


@doc raw""" 
    laspeyres(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray )

Laspeyres function with only vectors as inputs

```math
    P_{Laspeyres}(p^{0},p^{1},q^{0}) =
    \frac{\sum_{n=1}^{N}p^{1}_{n}q^{0}_{n}}{\sum_{n=1}^{N}p^{0}_{n}q^{0}_{n}}
```
"""
function laspeyres(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray)
    return sum(p1 .* q0) / (sum(p0 .* q0))
end
function laspeyres(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, q0::AbstractArray)
    return laspeyres(p1, p0, q0)
end

@doc raw"""
    paasche(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray )
Paasche function with only vectors as inputs

```math
    P_{Paasche}(p^{0},p^{1},q^{1}) =
    \frac{\sum_{n=1}^{N}p^{1}_{n}q^{1}_{n}}{\sum_{n=1}^{N}p^{0}_{n}q^{1}_{n}}

```

#### Arguments

- `p1::AbstractArray` Array with later period prices
- `p0::AbstractArray` Array with earlier period prices
- `q1::AbstractArray` Array with later period quantities

"""
function paasche(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray)
    return sum(p1 .* q1) / (sum(p0 .* q1))
end
function paasche(p1::AbstractArray, p0::AbstractArray,  q1::AbstractArray, args...)
    return paasche(p1, p0, q1)
end

@doc raw"""
    fisher(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray )
The Fisher price index, the geometric average of Laspeyres and Paasche indexes. The Fisher price index
was proposed by Irving Fisher as the "Ideal" price index based on its axiomatic properties. This index
is also a "superlative" index.
```math
     P_{Fisher}(p^{0},p^{1},q^{0},q^{1}) = \sqrt{P_{Paasche}\cdot P_{Laspeyres}}
```
"""
function fisher(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    return sqrt(laspeyres(p1, p0, q0) * paasche(p1, p0, q1))
end


@doc raw"""
    carli(p1::AbstractArray, p0::AbstractArray)
    
The Carli price index, the arithmetic mean of price relatives.
```math
    P_{Carli}(p^{1},p^{0}) =
    \frac{1}{N}\sum_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)
```
"""
function carli(p1::AbstractArray, p0::AbstractArray)
    return mean((p1 ./ p0))
end
function carli(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, ::AbstractArray)
    return carli(p1, p0)
end

@doc raw"""
    dutot(p1::AbstractArray, p0::AbstractArray)

The Dutot index, the ratio of the arithmetic means of prices in each period.

```math
  P_{Dutot}(p^{0},p^{1}) = \frac{\sum_{n=1}^{N}p^{1}_{n}}{\sum_{n=1}^{N}p^{0}_{n}}
```
"""
function dutot(p1::AbstractArray, p0::AbstractArray)
    return mean(p1) / mean(p0)
end
function dutot(p1::AbstractArray, p0::AbstractArray, args...)
    return dutot(p1, p0)
end

@doc raw"""
    jevons(p1::AbstractArray, p0::AbstractArray)
```math
  P_{Jevons}(p^{0},p^{1}) =
    \prod_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{(1/N)}
```
"""
function jevons(p1::AbstractArray, p0::AbstractArray)
    return geomean((p1 ./ p0))
end
function jevons(p1::AbstractArray, p0::AbstractArray, args...)
    return jevons(p1, p0)
end
@doc raw"""
    harmonic(p1::AbstractArray, p0::AbstractArray)

Also known as Coggeshall, the Harmonic mean of price relatives using the StatsBase function
```math
P_{Harmonic}(p^{0},p^{1}) = \bigl\{\frac{1}{N}\sum_{n=1}^{N} \left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{-1}\bigr\}^{-1}
```
"""
function harmonic(p1::AbstractArray, p0::AbstractArray)
    return harmmean((p1 ./ p0))
end
function harmonic(p1::AbstractArray, p0::AbstractArray, args...)
    return harmonic(p1, p0)
end

@doc raw"""
    cswd(p1::AbstractArray, p0::AbstractArray)

The Carruthers, Sellwood, Ward, Dalén index. A geometric average of the Carli and the 
harmonic mean of price relatives index.

```math
P_{CSWD}(p^{0},p^{1}) = \sqrt{P_{Carli}\cdot P_{Harmonic}}
```
"""
function cswd(p1::AbstractArray, p0::AbstractArray)
    return sqrt(harmonic(p1, p0) * carli(p1, p0))
end
function cswd(p1::AbstractArray, p0::AbstractArray, args...)
    return cswd(p1, p0)
end

@doc raw"""
    bmw(p1::AbstractArray, p0::AbstractArray)

Balk-Mehrhoff-Walsh (BMW) index, an index of price relatives weighted by
the square root of the inverse price relatives. The BMW is the linear approximation to the Jevons index. The BMW index was
developed in Mehrhoff (2007).  Balk (2008) also developed the index as the unweighted Walsh.

```math
P_{BMW}(p^{0},p^{1})=\sum_{n=1} \frac{p_{i}^{1}}{p_{i}^{0}} 
\frac{\sqrt{\frac{p_{i}^{0}}{p_{i}^{1}}}}{\sum_{n=1} \sqrt{\frac{p_{i}^{0}}{p_{i}^{1}}}}
```
"""
function bmw(p1::AbstractArray, p0::AbstractArray)
    numer = sum((p1 ./ p0) .* sqrt.(p0 ./ p1))
    denom = sum(sqrt.(p0 ./ p1))
    return numer / denom
end
function bmw(p1::AbstractArray, p0::AbstractArray, args...)
    return bmw(p1, p0)
end

@doc raw"""
    drobisch(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray )

```math
    P_{Drobisch}(p^{0},p^{1},q^{0},q^{1}) = (P_{Laspeyres}+P_{Paasche})/2
```

"""
function drobisch(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    return (laspeyres(p1, p0, q0) + paasche(p1, p0, q1)) / 2.0
end

@doc raw"""
    tornqvist(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

Törnqvist index, the geometric average of price relatives weighted by the arithmetic average of weight shares
between two periods. The Törnqvist is "superlative" and one of the most commonly used indexes when timely
quantity or expenditure information is available.

```math
     P_{Törnqvist}(p^{0},p^{1},q^{0},q^{1}) =
    \prod_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{\left(s^{0}_{n}+s^{1}_{n}\right)/2}
```

"""
function tornqvist(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    total_expend0 = sum(expend0)
    total_expend1 = sum(expend1)

    avg_share = 0.5 * ((expend1 ./total_expend1) + (expend0 ./ total_expend0))

    return exp(sum(log.( p1 ./ p0 ) .* avg_share))
end

@doc raw"""
    walsh(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

```math
    P_{Walsh}(p^{0},p^{1},q^{0},q^{1}) =
    \frac{\sum_{n=1}^{N}\sqrt{q^{0}_{n}q^{1}_{n}}\cdot
    p^{1}_{n}}{\sum_{n=1}^{N}\sqrt{q^{0}_{n}q^{1}_{n}}\cdot p^{0}_{n}}
```
"""
function walsh(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    avg_q = sqrt.(q1 .* q0)
    return sum(avg_q .* p1) / sum(avg_q .* p0)
end

@doc raw"""
    geary_khamis_bi(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

```math
    P_{Geary-Khamis}(p^{0},p^{1},q^{0},q^{1}) = \frac{\sum_{n=1}^{N}h(q^{0}_{n},
    q^{1}_{n})p^{1}_{n}}{\sum_{n=1}^{N}h(q^{0}_{n}, q^{1}_{n})p^{0}_{n}}
```
"""
function geary_khamis_bi(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    #Calculate the harmonic mean of q for each element of q1 and the corresponding q0
    harm_q = 1 ./ (1 ./ q1 + 1 ./ q0)
    return sum(harm_q .* p1) / sum(harm_q .* p0)
end

@doc raw"""
    palgrave(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray)

```math
    P_{Palgrave}(p^{0},p^{1},q^{0},q^{1}) =
    \sum_{n=1}^{N}s^{1}_{n}\frac{p^{1}_{n}}{p^{0}_{n}}
```

"""
function palgrave(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray)
    expend1 = p1 .* q1
    total_expend1 = sum(expend1)

    expend_share = (expend1 ./ total_expend1)

    return sum(expend_share .* p1 ./ p0) 
end
function palgrave(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, args...)
    return palgrave(p1, p0, q1)
end

@doc raw"""
    marshall_edgeworth(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

```math
    P_{Marshall-Edgeworth}(p^{0},p^{1},q^{0},q^{1}) =
    \frac{\sum_{n=1}^{N}p_{n}^{1}(q_{n}^{0}+q_{n}^{1})}{\sum_{n=1}^{N}p_{n}^{0}(q_{n}^{0}+q_{n}^{1})}
```
"""
function marshall_edgeworth(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    return sum((q0+q1) .* p1) / sum((q0+q1) .* p0)
end

@doc raw"""
    sato_vartia(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
```math
  P_{Sato-Vartia}(p^{0},p^{1},q^{0},q^{1}) =
\prod_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{
\frac{s^{1}_{n}-s^{0}_{n}}{\log (s^{1}_{n}) - \log
(s^{0}_{n})} \cdot \left ( \sum_{n=1}^{N}\frac{s^{1}_{n}-s^{0}_{n}}{\log (s^{1}_{n}) - \log
(s^{0}_{n})} \right )^{-1}
}
```
"""
function sato_vartia(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    total_expend0 = sum(expend0)
    total_expend1 = sum(expend1)

    weight_share1 = expend1 ./ total_expend1
    weight_share0 = expend0 ./ total_expend0

    sv_weight_share = vartia_L(weight_share1, weight_share0)
    #Normalize weight share
    sv_weight_share_norm = sv_weight_share/sum(sv_weight_share)

    return exp(sum(log.( p1 ./ p0 ) .* sv_weight_share_norm))
end
"""
    vartia_L(a::AbstractArray, b::AbstractArray)

    Helper function related to Sato-Vartia and montgomery-vartia
"""
function vartia_L(a::AbstractArray, b::AbstractArray)
    if b == a 
        return b
    else
        return (a-b) ./ (log.(a)-log.(b))
    end
end

"""
    vartia_L(a::Float64}, b::Float64)

montgomery-vartia uses the method with float, while sato-vartia uses the method with float vectors
"""
function vartia_L(a::Float64, b::Float64)
    if b == a 
        return b
    else
        return (a-b) ./ (log.(a)-log.(b))
    end
end
@doc raw"""
    montgomery_vartia(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

    Montgomery-Vartia or sometimes just Vartia or Vartia-I, but preceeded by Montgomery (1937).
    Implemented based on https://doi.org/10.1111/rssa.12633 von Auer and Wengenroth (2020).
"""
function montgomery_vartia(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1
    denom = vartia_L(sum(expend1), sum(expend0))
    
    return exp(sum(log.(p1 ./ p0).*vartia_L.(expend1, expend0)/denom))
end

@doc raw"""
    lloyd_moulton(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray, σ::Float64)

```math
     P_{Lloyd-Moulton}(p^{0},p^{1},q^{0}) =
\left[\sum_{n=1}^{N}s_{n}^{0}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{(1-\sigma)}\right]^{\left(\frac{1}{1-\sigma}\right)},
```
"""
function lloyd_moulton(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray, σ::Float64)
    expend0 = p0 .* q0
    total_expend0 = sum(expend0)
    weight_share0 = expend0 ./ total_expend0
    return sum(weight_share0  .* ((p1 ./ p0) .^ (1 - σ))) ^ (1 /(1-σ))
end
#Version of lloyd_moulton that takes a general call with kwargs and extracts the σ
function lloyd_moulton(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, q0::AbstractArray; kwargs...)
    σ = get(kwargs, :σ, 0.7)  #Set 0.7 as default
    return lloyd_moulton(p1, p0, q0, σ)
end

@doc raw"""
    stuvel(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

The Stuvel index unlike many superlative indexes, the Stuvel is consistent-in-aggregation. However,
it is not linearly homogeneous. Proposed by G. Stuvel in the 1957 *Econometrica* article
“A New Index Number Formula.” 
```math
    P_{Stuvel}(p^{{0}},p^{1},q^{{0}},q^{1}) = A + \sqrt{A^2 + V^{1}/V^{0}} \\
```
where
```math
    A = (P_{Laspeyres}-Q_{Laspeyres})/2 \\
    Q_{Laspeyres} = \frac{\sum_{n=1}^{N}q^{1}_{n}p^{0}_{n}}{\sum_{n=1}^{N}q^{0}_{n}p^{0}_{n}}\\
    V^{t} = p^{t} \cdot q^{t}
```
"""
function stuvel(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    # A is the difference between the laspeyres price and quuantity indexes, divided by 2.0
    A = (sum(p1 .* q0) / sum(p0 .* q0) - sum(p0 .* q1) / sum(p0 .* q0) )/2.0
    v1 = sum(p1 .* q1)
    v0 = sum(p0 .* q0)
    return A + sqrt(A^2.0 + v1/v0)
end

@doc raw"""
    lehr(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

```math
P_{Lehr}(p^{0},p^{1},q^{0},q^{1}) = \frac{\sum v_i^1}{\sum v_i^0} \cdot \frac{\sum \left( v_i^0 \left( v_i^0 + v_i^1 \right) \left( v_i^0 + \frac{v_i^1}{r} \right) \right)}
{\sum \left( v_i^1 \left( v_i^0 + v_i^1 \right) \left( v_i^0 r_i + v_i^1 \right) \right)}
```
where
```math
v_i^{t} = p_i^{t} \cdot q_i^{t} \\
r_i = \frac{p_i^1}{p_i^0}
```

"""
function lehr(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    numer = sum( (expend0 + expend1) .* q0 ./ (q0+q1) )
    denom = sum( (expend0 + expend1) .* q1 ./ (q0+q1) )

    return ( sum(expend1)/denom * numer/sum(expend0) )
end

@doc raw"""
    banerjee(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

On the factorial approach providing the true index of cost of living: an interpretation of a special pair of equations
K.S. Banerjee (1977).

```math
P_{Banerjee} = \frac{\sum v_i^1}{\sum v_i^0} \cdot 
\frac{\sum \left( v_i^0 \left( 1 + r_i \right) \right)}
{\sum \left( v_i^1 \left( 1 + \frac{1}{r_i} \right) \right)}
```
where
```math
v_i^{t} = p_i^{t} \cdot q_i^{t} \\
r_i = \frac{p_i^1}{p_i^0}
```
"""
function banerjee(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    return (sum(expend1)/sum(expend0)) * sum(expend0 .* (1.0 .+ p1 ./ p0))/sum(expend1 .* (1.0 .+ p0 ./ p1))
end

@doc raw"""
    davies(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)

```math
P_{Davies} = \frac{\sum v_i^1}{\sum v_i^0} \cdot 
\frac{\sum v_i^0 \sqrt{r_i}} 
{\sum \left( v_i^1 \sqrt{\frac{1}{r_i}} \right) } 
```
where
```math
v_i^{t} = p_i^{t} \cdot q_i^{t} \\
r_i = \frac{p_i^1}{p_i^0}
```
"""
function davies(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    expend1 = p1 .* q1

    return (sum(expend1)/sum(expend0)) * sum(expend0 .* sqrt.(p1 ./ p0))/sum(expend1 .* sqrt.(p0 ./ p1))
end

@doc raw"""
    lowe(p1::AbstractArray, p0::AbstractArray, q_b::AbstractArray,)

```math
      P_{Lowe}(p^{0},p^{1},q^{b}) =
    \frac{\sum_{n=1}^{N}p^{1}_{n}q^{b}_{n}}{\sum_{n=1}^{N}p^{0}_{n}q^{b}_{n}}
```

TBW
"""
function lowe(p1::AbstractArray, p0::AbstractArray, q_b::AbstractArray)
    return sum(p1 .* q_b) / (sum(p0 .* q_b))
end
function lowe(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, ::AbstractArray; kwargs...)
    q_b = get(kwargs, :q_b, 0)  
    return lowe(p1, p0, q_b)
end

@doc raw"""
    young(p1::AbstractArray, p0::AbstractArray; p_b::AbstractArray, q_b::AbstractArray)

```math
    P_{Young}(p^{0},p^{1},p^{b},q^{b}) =
    \sum_{n=1}^{N}s^{b}_{n}\frac{p^{1}_{n}}{p^{0}_{n}}
```
"""
function young(p1::AbstractArray, p0::AbstractArray; p_b::AbstractArray, q_b::AbstractArray)
    expendB = p_b .* q_b
    total_expendB = sum(expendB)
    weight_shareB = expendB ./ total_expendB

    return sum(weight_shareB .* (p1 ./  p0))
end

@doc raw"""
function young(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, ::AbstractArray; kwargs...)
    p_b = get(kwargs, :p_b, 0)  
    q_b = get(kwargs, :q_b, 0)  
    return young(p1, p0, p_b, q_b)
end
"""

@doc raw"""
    geolaspeyres(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray)

TBW
"""
function geolaspeyres(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray)
    expend0 = p0 .* q0
    total_expend0 = sum(expend0)
    weight_share0 = expend0 ./ total_expend0

    return exp(sum(weight_share0 .* log.(p1 ./ p0)))
end
function geolaspeyres(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, q0::AbstractArray)
    return geolaspeyres(p1, p0, q0)
end

@doc raw"""
    geopaasche(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray)

TBW
"""
function geopaasche(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray)
    expend1 = p1 .* q1
    total_expend1 = sum(expend1)
    weight_share1 = expend1 ./ total_expend1

    return exp(sum(weight_share1 .* log.(p1 ./ p0)))
end
function geopaasche(p1::AbstractArray, p0::AbstractArray, q1::AbstractArray, ::AbstractArray)
    return geopaasche(p1, p0, q1)
end

@doc raw"""
    ag_mean(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray, η::Float64)

The AG Mean is the weighted average of the Laspeyres and geometric Laspeyres indexes with η indicating the
relative weight on each. Lent and Dorfman (2009) introduced this formula and it is discussed in Armknect and
Silver (2012).

```math
P_{AGmean}(p^{0},p^{1},q^{0}) =
\eta^{\tau}
\prod_{n=1}^{N}\left(\frac{p^{1}_{n}}{p^{0}_{n}}\right)^{s^{0}_{n}}
+
\left (1 - \eta^{\tau} \right )\sum_{n=1}^{N}s^{0}_{n}\frac{p^{1}_{n}}{p^{0}_{n}}
```

"""
function ag_mean(p1::AbstractArray, p0::AbstractArray, q0::AbstractArray, η::Float64)
    return η * geolaspeyres(p1, p0, q0) + (1-η) * laspeyres(p1, p0, q0)
end
function ag_mean(p1::AbstractArray, p0::AbstractArray, ::AbstractArray, q0::AbstractArray; kwargs...)
    η = get(kwargs, :η, 0) 
    return ag_mean(p1, p0, q0, η)
end

function raito_of_harmonic_means()
    
end

function białek()
end

end # module BilateralIndexFormulas