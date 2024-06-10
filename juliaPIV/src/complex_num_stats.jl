"""
Author: Robin Shindelman
Date: 2024-06-10
Some basic statistics designed specifically for working with the wave PIV
program. I had to retool the basic statistics packages to be able to work
with the complex numbers in the PIV matrices.
"""
module ComplexNumStats
import Statistics: median, mean, std

"""
### im_median
    Find the median of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_median(collection)
    if length(collection) < 1
        return NaN
    end
    real_part = median(real.(collection))
    im_part = median(imag.(collection))
    return real_part + im_part * im
end

"""
### im_mean
    Find the mean of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_mean(collection)
    if length(collection) < 1
        return NaN
    end
    real_part = mean(real.(collection))
    im_part = mean(imag.(collection))
    return real_part + im_part * im
end

"""
### im_median
    Find the std dev of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_std(collection)
    if length(collection) < 1
        return NaN
    end
    real_part = std(real(collection), corrected=false)
    im_part = std(imag(collection), corrected=false)
    return real_part + im_part * im
end

# End module
end