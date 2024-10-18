using CUDA
using Images
using FileIO
using StaticArrays
using FFTW


function main()
    # Get images, upload to GPU
    im_pair = load_images()
    im1_cu = CuArray{Float32}(undef, size(im_pair[1]))
    im2_cu = CuArray{Float32}(undef, size(im_pair[2]))
    copyto!(im1_cu, im_pair[1])
    copyto!(im2_cu, im_pair[2])

    test_win = CUDA.zeros(Float32, size(im_pair[1]))

    # Other params
    win_size = Int32(64)
    ol = 0.5f0

    # Kernel config here:
    gpu_pass_launch(im1_cu, im2_cu, win_size, ol, test_win)
end

function load_images()
    im1 = load("../data/im1.jpg")
    im2 = load("../data/im2.jpg")
    crops = (24, 2424, 1, 2048)
    im1 = im1[crops[3]:crops[4], crops[1]:crops[2]]
    im2 = im2[crops[3]:crops[4], crops[1]:crops[2]]
    im_pair = (cu(Gray.(im1)), cu(Gray.(im2)))

    return im_pair
end

function gpu_pass_launch(im1::CuArray{Float32}, im2::CuArray{Float32},
                         win_size::Int32, ol::Float32, test_win::CuArray{Float32})
    
    kernel = @cuda launch=false gpu_pass!(im1, im2, win_size, ol, test_win)
    config = launch_configuration(kernel.fun)

    threads = min(length(im1), config.threads)
    blocks = cld(length(im1), threads)

    CUDA.@sync begin
        kernel(im1, im2, win_size, ol, test_win; threads, blocks)
    end
end

function gpu_pass!(im1::CuDeviceMatrix{Float32}, im2::CuDeviceMatrix{Float32},
                    win_size::Int32, ol::Float32, test_win::CuDeviceMatrix{Float32})

    # Thread for window center
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    # Map thread to window center in images
    M, N = size(im1)
    num_windows = M / win_size
    row = i / num_windows
    col = i % num_windows

    half_win = win_size / 32
    center_x = row * win_size + half_win
    center_y = col * win_size + half_win

    if center_x >= half_win && center_x <= (M - half_win) &&
        center_y >= half_win && center_y <= (N - half_win)

        @inbounds begin
            C = @view im1[center_x - half_win: center_x + half_win,
                          center_y - half_win: center_y + half_win] 
            D = @view im2[center_x - half_win: center_x + half_win,
                          center_y - half_win: center_y + half_win]
        end

        if i == 32
            @cushow center_x center_y
            # for i in 1:size(C, 1)
            #     for j in 1:size(C, 2)
            #         test_win[i, j] = C[i, j]
            #     end
            # end
        end
        


    end
        

    return
end

main()
