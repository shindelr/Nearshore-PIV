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

    # sort!(i, by = x -> (abs2(x), angle(x)))
    # no2 = n ÷ 2
    # isodd(n) && return i[no2+1]
    # return (i[no2] + i[no2+1]) / 2

    # v = partialsort!(i, div(n+1, 2, RoundDown):div(n+1, 2, RoundUp); by=abs2)
    v = partialsort!(i, div(n+1, 2, RoundDown):div(n+1, 2, RoundUp); by=x -> (abs2(x), angle(x)))

    return sum(v)/length(v)
end

# test_m = [
#     0.0 + 0.0im	0.0 + 0.0im	0.0 + 0.0im;
#     2.0 + 0.0im	NaN + 0.0im	0.0 + 0.0im;
#     2.0 + 1.0im	2.0 + 1.0im	0.0 + 0.0im
# ]

# test_m = [
#     NaN + NaN*im	NaN + NaN*im	NaN + NaN*im;
#     1.0 + 0.0im	NaN + 0.0im	2.0 + 0.0im;
#     1.0 + 0.0im	1.0 + 1.0im	2.0 + 1.0im
# ]

# jj:6 ii:9
# test_m = [
#     0.0 + 0.0im	0.0 + 1.0im	1.0 + 1.0im;
#     0.0 + 0.0im	NaN + 0.0im	1.0 + 0.0im;
#     0.0 + 0.0im	1.0 + 0.0im	1.0 + 0.0im
# ]

# Iter: jj:10, ii:3
# test_m = [
#     1.0 + 1.0im	1.0 + 1.0im	1.0 + 0.0im;
#     0.0 + 1.0im	NaN + 0.0im	0.0 + 0.0im;
#     0.0 + 1.0im	0.0 + 0.0im	0.0 + 0.0im
# ]
# Histo[10, 3] = 0.0 + 1.0im

# Iter: jj:21, ii:3
# test_m = [
#     1.0 + 0.0im	2.0 + 0.0im	1.0 + 0.0im;
#     1.0 + 1.0im	NaN + 0.0im	3.0 + 0.0im;
#     1.0 + 1.0im	2.0 + 0.0im	3.0 + 0.0im
# ]
# # Histo[21, 3] = 2.5 + 0.0im

# Iter: jj:39, ii:2
test_m = [
    NaN + NaN*im	4.0 + 2.0im	5.0 + 0.0im;
    NaN + NaN*im	NaN + 0.0im	4.0 + 3.0im;
    NaN + NaN*im	4.0 + 3.0im	5.0 + 3.0im
]
# Histo[39, 2] = 5.0 + 0.0im


im_median_magnitude(test_m[:])
