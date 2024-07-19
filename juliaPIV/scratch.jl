using Statistics

function im_median(collection)
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = median(real(i))
        im_part = median(imag(i))
        return real_part + im_part * im

        # return median(i)

    end

    # If all NaN
    return NaN
end

og_matrix = [
    0.5488135 + 0.71518937im  0.60276338 + 0.54488318im  0.42365480 + 0.64589411im;
    0.43758721 + 0.89177300im  0.96366276 + 0.38344152im  0.79172504 + 0.52889492im;
    0.56804456 + 0.92559664im  0.07103606 + 0.08712930im  0.02021840 + 0.83261985im
]

nan_matrix = [
    0.5488135 + 0.71518937im  NaN + NaN*im              0.42365480 + 0.64589411im;
    0.43758721 + 0.89177300im  0.96366276 + 0.38344152im  NaN + NaN*im;
    0.56804456 + 0.92559664im  NaN + NaN*im              0.02021840 + 0.83261985im
]

usum = im_median(nan_matrix[:])
