function firstpass(A, B, N, overlap, idx, idy)
    M = N[1]; N = N[2]
    sy, sx = size(A)
    xx = zeros(ceil((size(A, 1) - N) / ((1 - overlap) * N)) + 1, ceil((size(A, 2) - M) / ((1 - overlap) * M)) + 1)
    yy = similar(xx)
    datax = similar(xx)
    datay = similar(xx)
    IN = zeros(size(A))

    cj = 1
    for jj in 1:((1 - overlap) * N):sy - N + 1
        ci = 1
        for ii in 1:((1 - overlap) * M):sx - M + 1
            if IN[jj + div(N, 2), ii + div(M, 2)] != 1
                if isnan(idx[cj, ci])
                    idx[cj, ci] = 0
                end
                if isnan(idy[cj, ci])
                    idy[cj, ci] = 0
                end
                if jj + idy[cj, ci] < 1
                    idy[cj, ci] = 1 - jj
                elseif jj + idy[cj, ci] > sy - N + 1
                    idy[cj, ci] = sy - N + 1 - jj
                end
                if ii + idx[cj, ci] < 1
                    idx[cj, ci] = 1 - ii
                elseif ii + idx[cj, ci] > sx - M + 1
                    idx[cj, ci] = sx - M + 1 - ii
                end

                C = A[jj:jj + N - 1, ii:ii + M - 1]
                D = B[jj + idy[cj, ci]:jj + N - 1 + idy[cj, ci], ii + idx[cj, ci]:ii + M - 1 + idx[cj, ci]]
                C .-= mean(C)
                D .-= mean(D)
                stad1 = std(C)
                stad2 = std(D)

                if stad1 == 0
                    stad1 = NaN
                end
                if stad2 == 0
                    stad2 = NaN
                end

                R = xcorrf2(C, D) / (N * M * stad1 * stad2)

                if size(R, 1) == N - 1
                    max_y1, max_x1 = findmax(R)
                else
                    max_y1, max_x1 = findmax(R[div(N, 2) + 2:div(3 * N, 2) - 3, div(M, 2) + 2:div(3 * M, 2) - 3])
                end

                if length(max_x1) > 1
                    max_x1 = round(Int, sum(max_x1 .* (1:length(max_x1))) / sum(max_x1))
                    max_y1 = round(Int, sum(max_y1 .* (1:length(max_y1))) / sum(max_y1))
                elseif isempty(max_x1)
                    idx[cj, ci] = NaN
                    idy[cj, ci] = NaN
                    max_x1 = NaN
                    max_y1 = NaN
                end

                datax[cj, ci] = -(max_x1 - M) + idx[cj, ci]
                datay[cj, ci] = -(max_y1 - N) + idy[cj, ci]
                xx[cj, ci] = ii + div(M, 2)
                yy[cj, ci] = jj + div(N, 2)
                ci += 1
            else
                xx[cj, ci] = ii + div(M, 2)
                yy[cj, ci] = jj + div(N, 2)
                datax[cj, ci] = NaN
                datay[cj, ci] = NaN
                ci += 1
            end
        end
        cj += 1
    end
    return xx, yy, datax, datay
end

