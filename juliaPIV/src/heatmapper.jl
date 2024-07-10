using Plots, DelimitedFiles

j_array = readdlm("tests/juliaOut/normal_interp/scratch_samples.csv", ',', Int64)
m_array = readdlm("tests/mlabOut/normal_interp/m_fully_itpd.csv", ',', Float64)
unmodified_m_array = readdlm("tests/mlabOut/naninterp_testing/mtest_unmodifiedDATAX.csv", ',', Float64)


function get_grid(data)
    x = 1:size(data, 1)
    y = 1:size(data, 2)
    xgrid, ygrid = collect.((x, y))
    return xgrid, ygrid
end

xgrid_m, ygrid_m = get_grid(m_array)
xgrid_j, ygrid_j = get_grid(j_array)
xgrid_m_unmod, ygrid_m_unmod = get_grid(unmodified_m_array)




mlab = heatmap(ygrid_m, xgrid_m, m_array, title="Fully Interpolated Matlab data", c=:viridis)
jl = heatmap(ygrid_j, xgrid_j, j_array, title="Julia Scratch Samples", c=:viridis)
unmod_m = heatmap(ygrid_m_unmod, xgrid_m_unmod, unmodified_m_array, title="Unmodified Matlab Data", c=:viridis)

savefig(mlab, "tests/heatmaps/m_fully_itpd.png")
savefig(jl, "tests/heatmaps/scratch_samples.png")
