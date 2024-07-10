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