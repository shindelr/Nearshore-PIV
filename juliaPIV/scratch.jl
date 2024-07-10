using ScatteredInterpolation
using Interpolations
using DelimitedFiles

samples = [0.0; 0.5; 0.5; 0.5; 1.0]
points = [0.0 0.0; 0.0 1.0; 1.0 0.0; 0.5 0.5; 1.0 1.0]'

n = 5
x = range(0, stop = 1, length = n)
y = range(0, stop = 1, length = n)
X = repeat(x, n)[:]
Y = repeat(y', n)[:]
gridPoints = [X Y]'

itp = interpolate(Multiquadratic(), points, samples)
interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
gridded = reshape(interpolated, n, n)
# ---------------------------------------------------------------------

samples = readdlm("tests/juliaOut/naninterp_testing/jtest_DATAX.csv", ',', Float64)
xs = collect(33:32:3041)
ys = collect(33:32:2017)
itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
extp = Interpolations.extrapolate(itp, Line())

XI = collect(17.0:16.0:3057.0)
YI = collect(17.0:16.0:2033.0)

# Interpolate interior first
# interior_itp = zeros(Float64, (length(ys), length(xs)))
# extrapolated = zeros(Float64, (length(YI), length(XI)))
# interior_itp = copy(extrapolated)
itp_results = zeros(Float64, (length(YI), length(XI)))

for (mi, yi) in enumerate(YI)
    for (ni, xi) in enumerate(XI)
        if (33 < xi < 3041) && (33 < yi < 2017)
            # interior_itp[mi, ni] = itp(yi, xi)
            itp_results[mi, ni] = itp(yi, xi)
        else
            # extrapolated[mi, ni] = extp(yi, xi)
            itp_results[mi, ni] = extp(yi, xi)
        end
    end
end

samples = round.(Int, itp_results)
writedlm("tests/juliaOut/normal_interp/scratch_samples.csv", samples,',')

# ====================== ScatteredInterpolation.jl ==========================
                # # Collect coarse points to interpolate on
                # points = collect.(Iterators.product(X, Y))
                # points_coords = hcat([Tuple(p) for p in points]...)
                # points_coords_matrix = hcat([i[1] for i in points_coords], [i[2] for i in points_coords])   

                # display(points_coords_matrix)

                # # Collect values to interpolate using coords above
                # points_vals_datax = vec(datax)
                # points_vals_datay = vec(datay)

                # # Interpolate!
                # itp_x = interpolate(InverseMultiquadratic(), points_coords_matrix', points_vals_datax)
                # itp_y = interpolate(InverseMultiquadratic(), points_coords_matrix', points_vals_datay)

                # # Create fine grid
                # fine_grid_points = collect.(Iterators.product(XI, YI))
                # fine_coords = hcat([Tuple(p) for p in fine_grid_points]...)
                # fine_coords_m = hcat([i[1] for i in fine_coords], [i[2] for i in fine_coords])   

                # interpolated_datax = ScatteredInterpolation.evaluate(itp_x, fine_coords_m')


                # for p in fine_grid_points
                #         itp_val_vec_x = ScatteredInterpolation.evaluate(itp_x, p)
                #         itp_val_vec_y = ScatteredInterpolation.evaluate(itp_y, p)
                #         # itp_val = itp_val_vec[1]
                #         # sample[c] = itp_val
                # end
# ====================== ScatteredInterpolation.jl ==========================

# ====================== Interpolations.jl ==========================
                # Creating new NaN matrices to interpolate into
                # m = length(YI)
                # n = length(XI)
                # itp_datax = fill(NaN, m+2, n+2)
                # itp_datay = fill(NaN, m+2, n+2)

                # Build interpolate func, layer onto datax, copy into NaN matrix
                # itp_x = Interpolations.interpolate((Y, X), datax, Gridded(Linear()))
                # datax = [itp_x(yi, xi) for yi in YI, xi in XI]
                
                # itp_y = Interpolations.interpolate((Y, X), datay, Gridded(Linear()))
                # datay = [itp_y(yi, xi) for yi in YI, xi in XI]
                
                # # itp_datax[2:end-1, 2:end-1] = round.(Int, datax)
                # # itp_datay[2:end-1, 2:end-1] = round.(Int, datay)

                # # Interpolate out the NaN's we put in
                # # datax, datay = linear_naninterp(itp_datax, itp_datay)

                # datax = round.(Int, datax)
                # datay = round.(Int, datay)
# ====================== Interpolations.jl ==========================