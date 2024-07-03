using ScatteredInterpolation

sample = [23 51 10 NaN 8;
          19 20 NaN NaN 3; 
          103 7 7 113 23]

coords = findall(x -> isnan(x), sample)
coords_M = hcat([i[1] for i in coords], [i[2] for i in coords])

itp = interpolate(Multiquadratic(), coords_M', sample) 

for c in coords
    c = [c[1]; c[2]]
    @show sample[c]
    @show c
    sample[c] = evaluate(itp, c)
    # @show sample[c]
end
