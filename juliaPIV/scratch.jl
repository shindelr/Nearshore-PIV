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

function im_median_magnitude(collection)
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN

    sort!(i; by=abs)
    n = length(i)
    no2 = n ÷ 2
    isodd(n) && return i[no2 + 1]
    return (i[no2] + i[no2+1]) / 2
end


og_matrix = [
    0.5488135 + 0.71518937im  0.60276338 + 0.54488318im  0.42365480 + 0.64589411im;
    0.43758721 + 0.89177300im  0.96366276 + 0.38344152im  0.79172504 + 0.52889492im;
    0.56804456 + 0.92559664im  0.07103606 + 0.08712930im  0.02021840 + 0.83261985im
]

nan_matrix = [
    0.5488135 + 0.71518937im  NaN 0.42365480 + 0.64589411im;
    0.43758721 + 0.89177300im  0.96366276 + 0.38344152im  NaN;
    0.56804456 + 0.92559664im  NaN 0.02021840 + 0.83261985im
]

# nan_m = [
#     0.5488 + 0.7152im      NaN + 0.0000im   0.4237 + 0.6459im;
#     0.4376 + 0.8918im   0.9637 + 0.3834im     NaN + 0.0000im;
#     0.5680 + 0.9256im      NaN + 0.0000im   0.0202 + 0.8326im
# ]

im_part = [0.7152, 0.8918, 0.9256, 0.3834, 0.6459, 0.8326]
real_part = [0.5488, 0.4376, 0.5680, 0.9637, 0.4237, 0.0202]
# median(real_part)

im_median_magnitude(nan_matrix[:])
