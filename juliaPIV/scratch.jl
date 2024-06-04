if size(R, 1) == (N - 1)
    max_val, max_coords = findmax(R)
    max_y1, max_x1 = Tuple(max_coords)
else
    max_val, max_coords = findmax(
                    R[Int(floor(.5*N+2)):Int(floor(1.5*N-3)), 
                    Int(floor(.5*M+2)):Int(floor(1.5*M-3))])
    max_y1, max_x1 = Tuple(max_coords)
end




# Handle multiple or no maximum values
if length(max_coords) > 1
    max_x1 = round(
        Int, sum([c[2] * i for (i, c) in enumerate(max_coords)]) / sum(c[2] for c in max_coords)
        )


    max_y1 = round(Int, sum([c[1] * i for (i, c) in enumerate(max_coords)]) / sum(c[1] for c in max_coords))
elseif isempty(max_coords)
    idx[cj, ci] = NaN
    idy[cj, ci] = NaN
    max_x1 = NaN
    max_y1 = NaN
else
    max_y1, max_x1 = max_coords[1]  # Unpack the single tuple to max_y1 and max_x1
end