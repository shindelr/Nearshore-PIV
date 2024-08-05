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
    heatmap(ygrid, xgrid, data, title=title, c=:viridis)
end

function display_3d(path, title)
    data = read_in_data(path)
    xgrid, ygrid = get_grid(data)
    plot(ygrid, xgrid, data, st=:surface, c=:viridis, title=title)
end


path = "tests/juliaOut/multipass_loop/penultimate_datax.csv"
title = "Julia penultimate_datax w/Modified grids"

path = "tests/mlabOut/multipass_loop/reginterp_datax.csv"
title = "Mlab reginterp_datax"

display_heatmap(path, title)
display_3d(path, title)
