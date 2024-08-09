"""
Robin Shindelman
2024-08-06

A brief interpolation script, meant to be used in conjunction
with the cross-validation script, interpolation_cv.py.
"""
module Interp 
    
using Interpolations


function build_grids_2(data)
    coarse_y_dim = size(data, 1)
    coarse_x_dim = size(data, 2)

    min_y, max_y = minimum(data[:, 1]), maximum(data[:, 1])
    min_x, max_x = minimum(data[1, :]), maximum(data[1, :])
    coarse_ys = LinRange(min_y, max_y, coarse_y_dim)
    coarse_xs = LinRange(min_x, max_x, coarse_x_dim)

    fine_yi_dim = (coarse_y_dim * 2) + 1
    fine_xi_dim = (coarse_x_dim * 2) + 1
    fine_YI = LinRange(min_y, max_y, fine_yi_dim)
    fine_XI = LinRange(min_x, max_x, fine_xi_dim)

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end

function regular_interp(samples, xs, ys, XI, YI)
    # samples = convert(AbstractArray{Float64, 2}, samples)
    xs = convert(Vector{Float64}, xs)
    ys = convert(Vector{Float64}, ys)
    XI = convert(Vector{Float64}, XI)
    YI = convert(Vector{Float64}, YI)

    @show size(xs) size(ys) size(samples)
    
    itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
    # itp_results = zeros(Float64, (length(YI), length(XI)))
    # itp_results = [itp(yi, xi) for yi in YI, xi in XI]
    # return itp_results
end

function test_interp(data)
    Y, X, YI, XI = build_grids_2(data)
    itpd_data = regular_interp(data, X, Y, XI, YI)
    return itpd_data
end

end