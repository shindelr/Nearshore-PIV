using Plots
using DelimitedFiles
using Interpolations

function regular_interp(samples, xs, ys, XI, YI)
    itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
    # extp = Interpolations.extrapolate(itp, Line())

    itp_results = zeros(Float64, (length(YI), length(XI)))

    for (mi, yi) in enumerate(YI)
        for (ni, xi) in enumerate(XI)
            # Interpolate the interior of the matrix
            # try 
                itp_results[mi, ni] = itp(yi, xi)
                # samples[mi, ni] = itp(yi, xi)
            # Otherwise, extrapolate
            # catch BoundsError
            #     println("Extrapolated at $mi, $ni\n")
            #     itp_results[mi, ni] = extp(yi, xi)
            # end
        end
    end
    return itp_results
end

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

function peaks(x, y)
    z = 3 * (1 - x)^2 * exp(-x^2 - (y + 1)^2) - 10 * (x / 5 - x^3 - y^5) * exp(-x^2 - y^2) - 1 / 3 * exp(-(x + 1)^2 - y^2)
    return z
end

function build_grids_2(data)
    coarse_y_dim = size(data, 1)
    coarse_x_dim = size(data, 2)
    coarse_ys = LinRange(0, 1, coarse_y_dim)
    coarse_xs = LinRange(0, 1, coarse_x_dim)

    fine_yi_dim = (coarse_y_dim * 2) + 1
    fine_xi_dim = (coarse_x_dim * 2) + 1
    fine_YI = LinRange(0, 1, fine_yi_dim)
    fine_XI = LinRange(0, 1, fine_xi_dim)

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end

subset = [  0.0463474    0.0227534   -0.00837549  -0.0266841   0.00601272    0.00298008    0.00456666  -0.0509503     0.00785566  -0.0130272   -0.0123018   -0.042495    -0.0558764;
0.0293588    0.0432518    0.0587956    0.0153987   0.0384511     0.0552958     0.073269         0.0633309     0.0208217    0.0346896    0.00275334  -0.0184351    0.00800629;
0.0330418    0.0488336    0.0386758    0.0088719   0.0112397     0.0730536     0.065421         0.0692089     0.0231118    0.0320262    0.00790876  -0.00892043   0.0419694;
0.00507062  -0.0177193    0.00156691   0.0117704   0.0196212     0.0536875     0.0216615        0.0940214     0.0159455   -0.00724363   0.0169934    0.0269343    0.0245413;
0.00815702  -0.0112511    0.0203176    0.0944224   0.0550089     0.0515141     0.0341976        0.145847      0.00290203   0.0393136    0.00961366   0.00208987  -0.0160093;
0.0441218    0.0282753    0.0860296    0.0597418   0.122801    NaN           NaN           NaN           NaN           -0.0002102   -0.0230424    0.0117914   -0.0705229;
0.0527618   -0.00293387   0.0364239    0.126171    0.147663    NaN           NaN              NaN           NaN            0.0738228    0.0260172    0.0376218   -0.0690405;
0.0511356    0.0461903    0.0927345    0.164809    0.166869    NaN           NaN              NaN           NaN            0.106338     0.0836875    0.0579849   -0.070148;
0.0316338    0.0662742    0.0770998    0.179601    0.192132    NaN           NaN              NaN           NaN            0.173396     0.130567     0.0410564   -0.0158512;
-0.012834     0.132887     0.117992     0.157399    0.150163    NaN           NaN              NaN           NaN            0.040521     0.0799092    0.0637829    0.0440173;
0.00538126   0.0227474    0.101663     0.0855544   0.0424873   NaN           NaN           NaN           NaN            0.0240775    0.0551653   -0.0187337   -0.00780226;
-0.0921788   -0.0169682    0.0558944   -0.0239361  -0.0279447   NaN           NaN              NaN           NaN            0.0139158   -0.0114896   -0.060729    -0.0109905;
-0.0806385   -0.0676329   -0.0374586   -0.0938637  -0.078249      0.0206674     0.0228627        0.0737724     0.0684696    0.00655972   0.00428526  -0.014801    -0.0224422;
-0.0496796   -0.0943678   -0.10319     -0.0948814  -0.0810961    -0.0152199     0.00673767      -0.00118545   -0.0577585   -0.056953    -0.0548346   -0.0634352   -0.0335631;
-0.0102227   -0.109281    -0.117695    -0.0841187  -0.0909032    -0.0604942    -0.0637831       -0.0247401    -0.0615724   -0.0804332   -0.0886972   -0.0481862   -0.0176157;
0.00496655  -0.0424937   -0.101847    -0.0973482  -0.0540484    -0.0562604     0.00845846  -0.0348084    -0.0262642   -0.0505313   -0.0660822   -0.0235215   -0.0130509]

maximum(subset)