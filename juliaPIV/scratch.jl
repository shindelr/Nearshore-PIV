using Plots
using DelimitedFiles
using Interpolations

function regular_interp(samples, xs, ys, XI, YI)
    itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
    # extp = Interpolations.extrapolate(itp, Line())

    itp_results = zeros(Float64, (length(YI), length(XI)))

    for (mi, yi) in enumerate(YI)
        for (ni, xi) in enumerate(XI)
            # Interpolate the interior of the matrix
            # try 
                itp_results[mi, ni] = itp(yi, xi)
                # samples[mi, ni] = itp(yi, xi)
            # Otherwise, extrapolate
            # catch BoundsError
            #     println("Extrapolated at $mi, $ni\n")
            #     itp_results[mi, ni] = extp(yi, xi)
            # end
        end
    end
    return itp_results
end

function build_grids(wins, overlap, sx, sy, i)
    next_win_x = wins[i + 1, 1]
    next_win_y = wins[i + 1, 2]

    # Final window size is duplicated, so check for equality.
    if wins[i, 1] != next_win_x
        X = (1:((1 - overlap) * 2 * next_win_x):
                sx - 2 * next_win_x + 1) .+ next_win_x
        XI = (1:((1 - overlap) * next_win_x):
                sx - next_win_x + 1) .+ (next_win_x / 2)
    else
        X = (1:((1 - overlap) * next_win_x):
                sx - next_win_x + 1) .+ (next_win_x / 2)
        XI = (1:((1 - overlap) * next_win_x):
                sx - next_win_x + 1) .+ (next_win_x / 2)
        X = copy(XI)
    end

    if wins[i, 2] != next_win_y
        Y = (1:((1 - overlap) * 2 * next_win_y): 
                sy - 2 * next_win_y + 1) .+ next_win_y
        YI = (1:((1 - overlap) * next_win_y):
                sy - next_win_y + 1) .+ (next_win_y / 2)
    else
        Y = (1:((1 - overlap) * next_win_y):
                sy - next_win_y + 1) .+ (next_win_y / 2)
        YI = (1:((1 - overlap) * next_win_y):
                sy - next_win_y + 1) .+ (next_win_y / 2)
        Y = copy(YI)
    end

return X, Y, XI, YI
end

function peaks(x, y)
    z = 3 * (1 - x)^2 * exp(-x^2 - (y + 1)^2) - 10 * (x / 5 - x^3 - y^5) * exp(-x^2 - y^2) - 1 / 3 * exp(-(x + 1)^2 - y^2)
    return z
end

function build_grids_2(data)
    coarse_y_dim = size(data, 1)
    coarse_x_dim = size(data, 2)
    coarse_ys = LinRange(0, 1, coarse_y_dim)
    coarse_xs = LinRange(0, 1, coarse_x_dim)

    fine_yi_dim = (coarse_y_dim * 2) + 1
    fine_xi_dim = (coarse_x_dim * 2) + 1
    fine_YI = LinRange(0, 1, fine_yi_dim)
    fine_XI = LinRange(0, 1, fine_xi_dim)

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end

function intpeak(x1, y1, R, Rxm1, Rxp1, Rym1, Ryp1, N)
    if length(N) == 2
        M = N[1]; N = N[2]
    else
        M = N
    end
    x01 = x1 + ((log(Complex(Rxm1)) - log(Complex(Rxp1))) / ((2 * log(complex(Rxm1))) - (4 * log(R)) + (2 * log(complex(Rxp1)))))
    y01 = y1 + ((log(Complex(Rym1)) - log(Complex(Ryp1))) / ((2 * log(complex(Rym1))) - (4 * log(R)) + (2 * log(complex(Ryp1)))))
    x0 = x01 - M
    y0 = y01 - N

    x0 = real(x0)
    y0 = real(y0)

    return x0, y0
end

max_x1, max_y1 = 16, 16
R = 0.1684079708669066
Rxm1 = -8.42042795691657e-17
Rxp1 = 0.022767921593265242
Rym1 = 0.013138210493162958
Ryp1 = -0.020155209279284103
N = 16.0

intpeak(max_x1, max_y1, R, Rxm1, Rxp1, Rym1, Ryp1, N)
