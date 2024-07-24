using Statistics

function im_median(collection)
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = median(real(i))
        im_part = median(imag(i))
        return real_part + im_part * im
    end

    # If all NaN
    return NaN
end

"""
## im_median_magnitude
    Take the median of a collection of complex numbers using the absolute magnitude.
    This great function was created by the Julia Community, specifically:
    @PeterSimmon & @mbauman
"""
function im_median_magnitude(collection::AbstractArray{Complex{T}}) where {T}
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN
    n = length(i)
    v = partialsort!(i, div(n+1, 2, RoundDown):div(n+1, 2, RoundUp); by=abs2)
    return sum(v)/length(v)
end

nan_matrix = [
    0.5488135 + 0.71518937im  NaN 0.42365480 + 0.64589411im;
    0.43758721 + 0.89177300im  0.96366276 + 0.38344152im  NaN;
    0.56804456 + 0.92559664im  NaN 0.02021840 + 0.83261985im
]

test_m = [
    1.0 + 0.0im	1.0 + 0.0im	1.0 + 0.0im;
    1.0 + 0.0im	NaN + 0.0im	0.0 + 0.0im;
    0.0 + 0.0im	0.0 + 0.0im	0.0 + 0.0im
]


im_median_magnitude(test_m[:])
