using Plots, DelimitedFiles

j_array = readdlm("tests/juliaOut/naninterp_testing/jtest_roundedInverseMultQDATAX.csv", ',', Float64,)
m_array = readdlm("tests/mlabOut/naninterp_testing/mtest_1stNaNinterpDATAX.csv", ',', Float64)
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




mlab = heatmap(ygrid_m, xgrid_m, m_array, title="First NaNInterp Matlab data", c=:viridis)
jl = heatmap(ygrid_j, xgrid_j, j_array, title="Julia Inverser Multiquadratic", c=:viridis)
unmod_m = heatmap(ygrid_m_unmod, xgrid_m_unmod, unmodified_m_array, title="Unmodified Matlab Data", c=:viridis)

savefig(mlab, "tests/heatmaps/1stNaNInterpMlab.png")
savefig(jl, "tests/heatmaps/InvQMultJulia.png")
