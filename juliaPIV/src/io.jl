using VideoIO
using Base.Threads
using Images
using Plots
include("main.jl")
using .JuliaPIV

function vid(f)
    frame_pairs = []
    im1 = read(f)
    i = 0
    while !eof(f)
        im2 = read(f)
        
        g_im1 = Gray.(im1)
        g_im2 = Gray.(im2)

        if i % 8 == 0 && i > 0 
            push!(frame_pairs, (g_im1, g_im2))
            im1 = im2
        end

        i += 1
    end
    return frame_pairs
end

function thread_it()
    io = VideoIO.open("tests/test_data/IMG_6566.avi")
    f = VideoIO.openvideo(io)
    pairs = vid(f)

    println("Got pairs, starting threads")

    Threads.@threads for (frame1, frame2) in pairs
        JuliaPIV.main(frame1, frame2)
    end
end

@show Threads.nthreads()
thread_it()