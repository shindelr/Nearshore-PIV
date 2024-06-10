"""
All the filtering functions for the wave PIV program are grouped together here
in this module.
"""
module PIVFilters
include("./complex_num_stats.jl")
using .ComplexNumStats: im_median, im_mean, im_std

"""
### localfilt
    Filter out vectors that deviate from the median or the mean of their \
    surrounding neighbors by the factor `threshold` times the standard \
    deviation of the neighbors.\n
    Parameters:
    -----------
    - x, y, u, v : `Matrices`
    - threshold : `Int`
        Specifies the point at which a vector has deviated too 
        far from the specified statistical mean or median.
    - median_bool: `Bool`
        If true, specifies that the median should be the turning
        point for the data to be filtered out on. Defaults
        to true. If specified as false, the mean value will be 
        used instead.
        method.
    - m : `Int`
        Defines the number of vectors contributing to the median 
        or mean value of each vector. Defaults to 3, though the
        original implementation mentions that 5 is a good number
        too. Also known as "kernelsize"
    - mask : `Matrix`
        Use to mask out areas of the given matrices to improve
        computation times. Default to an empty matrix.
    Returns:
    --------
    - hu, hv : `Matrices`
            Successfully filtered matrices. New versions of u
            and v.
"""
function localfilt(x, y, u, v, threshold, median_bool=true, m=3, mask=[])
    method =  median_bool ? "median" : "mean"
    IN = zeros(eltype(u), size(u))
    # !!!! Should handle mask being a file here !!!! #

    nu_dim1 = round(Int, size(u, 1) + 2 * floor(m/2))
    nu_dim2 = round(Int, size(u, 2) + 2 * floor(m/2))
    nu = zeros(eltype(u), (nu_dim1, nu_dim2)) * NaN
    nv = zeros(eltype(u), (nu_dim1, nu_dim2)) * NaN
    
    # Transfer over data
    from_cols = round(Int, floor(m/2) + 1)
    minus_rows = round(Int, floor(m/2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v

    # Testing: same five errors detected from before with datax, datay
    # writedlm("tests/juliaOut/JtestNV.csv", nv, ',')
    
    INx = zeros(eltype(nu), size(nu))
    INx[from_cols: end - minus_rows, from_cols: end - minus_rows] = IN
    # Testing: Success!
    
    # Could be a little problem area here. Not sure any of these vars are used.
    prev = isnan.(nu)
    previndex = findall(prev)
    teller = true

    U2 = nu .+ im .* nv
    # Testing: U2 Looks okay, but might not be, it's hard to tell with im's. 
    ma, na = size(U2)
    histostd = zeros(ComplexF64, size(nu)) 
    histo = zeros(ComplexF64, size(nu))
    hista = zeros(eltype(nu), size(nu)) 
    histastd = zeros(eltype(nu), size(nu)) 

    iter = ProgressBar(m - 1:na - m + 2)
    for p in iter  # Looks gnar, but just a bar!
        for ii in m - 1:1:na - m + 2
            for jj in m - 1:1:ma - m + 2

                if INx[jj, ii] != 1
                    m_floor_two = floor(Int, m / 2)
                    tmp = U2[round(Int, jj - m_floor_two): round(Int, jj + m_floor_two),
                            round(Int, ii - m_floor_two): round(Int, ii + m_floor_two)] 
                    tmp[ceil(Int, m / 2), ceil(Int, m / 2)] = NaN;

                    # Create a collection of all elements without NaN values
                    usum_prep = collect(Skipper.skip(x -> isnan(x), tmp[:]))
                    
                    # !!!!!!!! NEEDS TESTING STILL !!!!!!!! #
                    # Run the appropriate stat depending on method arg.
                    usum = median_bool ? im_median(usum_prep) : mean(usum_prep)
                    histostd[jj, ii] = im_std(usum_prep)

                else
                    usum = tmp = histostd[jj, ii] = NaN
                end
                histo[jj, ii] = usum
            end
        end 
        set_description(iter, "Local $method filter running: ")
    end

    # TESTING: Success!
    # writedlm("tests/juliaOut/JtestHISTOSTD.csv", histostd, ',')
    
    # Locate gridpoints w/higher value than the threshold
    coords = findall(
        (real(U2) .> real(histo) .+ threshold .* real(histostd)) .|
        (imag(U2) .> imag(histo) .+ threshold .* imag(histostd)) .|
        (real(U2) .< real(histo) .- threshold .* real(histostd)) .|
        (imag(U2) .< imag(histo) .- threshold .* imag(histostd)))
    
    # Then "filter" those points out by changing them to NaN!
    for jj in eachindex(coords)
        nu[coords[jj]] = NaN
        nv[coords[jj]] = NaN
    end

    # TESTING: Showing 13 differences, nu & nv are directly related to datax
    #  and datay. So maybe this problem stems from firstpass()

    # writedlm("tests/juliaOut/JtestNUFilt.csv", nu, ',')
    # writedlm("tests/juliaOut/JtestNVFilt.csv", nv, ',')

    # Skipped print statement about how many vectors were filtered.
    # Skpped checking for 'interp' arg, because the actual program wasn't using
    # it. We call naninterp explicitly right after this function.

    m_ceil_two = ceil(Int, m/2)
    m_floor_two = floor(Int, m/2)
    hu = nu[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]
    hv = nv[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]
    return hu, hv
end

"""
### naninterp
    Interpolates NaN's in a vectorfield. Sorts all spurious 
    vectors based on the number of spurious neighbors to a 
    point. The function replaces NaN values in the input 
    vector fields U and V using linear interpolation.\n
    Interpolation starts with the ones that have the least 
    number of outliers in their neighborhood and loops until no 
    NaN's are present in the field.\n
    NOTE: This function is a completely gutted version of the
    original naninterp combined with naninterp2. The majority
    of naninterp was varargin code figuring out which method
    to execute naninterp2 with. Since arguments to naninterp
    were hardcoded in, the desired form of naninterp was just
    turned into a specific function. If necessary, other 
    methods can be written to emulate the functionality of
    the original functions. Moreover, though a mask was being
    passed in, the masking branch was not being executed in
    either functions. Therefore, in this implementation, masking
    was left out.\n

    Parameters:
    -----------
    - u, v: `Matrices`\n

    Original author:
    ----------------
    J. Kristian Sveen (jks@math.uio.no)
    Department of Mathematics, Mechanics Division, University of Oslo, Norway
    Copyright 1999 - 2001
    For use with MatPIV 1.6, Copyright
    Distributed under the terms of the GNU - GPL license
"""
function linear_naninterp(u, v)
    
end

# End module
end