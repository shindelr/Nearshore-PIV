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