# Not currently used/functioning
function naninterp(sample, pass)
    nan_coords = findall(x -> isnan(x), sample)
    non_nan_coords= findall(x -> !isnan(x), sample)
    non_nan_coords_matrix = hcat([i[1] for i in non_nan_coords], [i[2] for i in non_nan_coords])
    non_nan_vals= [sample[c] for c in non_nan_coords]
    println("made it: 479") # Last place it gets past

    # Debugging
    if pass == 2
        # writedlm("tests/juliaOut/second_pass_interp/stuck_sample.csv", sample, ',')
        display(non_nan_coords_matrix)
        display(non_nan_vals)
    end

    itp = ScatteredInterpolation.interpolate(InverseMultiquadratic(), non_nan_coords_matrix', non_nan_vals)
    for c in nan_coords
        c_extracted = [c[1]; c[2]]
        itp_val_vec = ScatteredInterpolation.evaluate(itp, c_extracted)
        itp_val = itp_val_vec[1]
        sample[c] = itp_val
    end

    return sample
end
