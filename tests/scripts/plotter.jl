using Plots, DelimitedFiles

function get_grid(data)
    x = 1:size(data, 1)
    y = 1:size(data, 2)
    xgrid, ygrid = collect.((x, y))
    return xgrid, ygrid
end

function read_in_data(path)
    # Path must be to a comma delimited csv.
    data = readdlm(path, ',', Float64)
    return data
end

function display_heatmap(path, title)
    data = read_in_data(path)
    xgrid, ygrid = get_grid(data)
    heatmap(ygrid, xgrid, data, title=title, c=:viridis, clims=(-25., 25.), ylims=(1, 127), xlims=(1, 191))
end

function display_3d(path, title)
    data = read_in_data(path)
    xgrid, ygrid = get_grid(data)
    plot(ygrid, xgrid, data, st=:surface, c=:viridis, title=title)
end


path = "tests/juliaOut/datax1.csv"
title = "Julia interp_datax1"

display_heatmap(path, title)

path = "tests/mlabOut/datax1.csv"
title = "Mlab datax1"

display_heatmap(path, title)
# display_3d(path, title)

