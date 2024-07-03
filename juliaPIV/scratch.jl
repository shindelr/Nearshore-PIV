using ScatteredInterpolation

sample = [23 51 10 NaN 8;
          19 20 NaN NaN 3; 
          103 7 7 113 23]

nan_coords = findall(x -> isnan(x), sample)

non_nan_coords = findall(x -> !isnan(x), sample)
non_nan_coords_M = hcat([i[1] for i in non_nan_coords], [i[2] for i in non_nan_coords])
non_nan_vals = [sample[c] for c in non_nan_coords]

itp = interpolate(Multiquadratic(), non_nan_coords_M', non_nan_vals) 

val_vec = evaluate(itp, [2,4])
val = val_vec[1]

for c in nan_coords
    c_m = [c[1]; c[2]]
    itp_val_vec = evaluate(itp, c_m)
    itp_val = itp_val_vec[1]
    @show itp_val
    sample[c] = itp_val
end

@show sample