"""
        build_grids(wins, overlap, sx, sy, i)

    Build grids for interpolation.

    Arguments
    ---------
    - `wins::Array{Float64,2}`: A matrix containing the window sizes for each iteration.
    - `overlap::Float64`: The overlap between adjacent windows.
    - `sx::Int`: The size of the x-axis.
    - `sy::Int`: The size of the y-axis.
    - `i::Int`: The index of the current window.

    Returns
    ---------
    - `X::Array{Float64,1}`: The x-coordinates of the coarse grid points.
    - `Y::Array{Float64,1}`: The y-coordinates of the coarse grid points.
    - `XI::Array{Float64,1}`: The x-coordinates of the fine interpolated grid points.
    - `YI::Array{Float64,1}`: The y-coordinates of the fine interpolated grid points.
"""
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