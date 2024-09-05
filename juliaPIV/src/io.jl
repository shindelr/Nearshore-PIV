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

        if i % 6 == 0 && i > 0 
            push!(frame_pairs, (g_im1, g_im2))
            im1 = im2
        end

        i += 1
    end
    return frame_pairs
end

# function thread_it()
#     io = VideoIO.open("tests/test_data/IMG_6566.avi")
#     f = VideoIO.openvideo(io)
#     pairs = vid(f)

#     println("Got pairs, starting threads")

#     i = 0
#     chunks = Iterators.partition(pairs, length(pairs) ÷ Threads.nthreads())
#     tasks = map((chunk1, chunk2)) do chunk1, chunk2
#         Threads.@spawn JuliaPIV.main(chunk1, chunk2)
#     end
#     chunk_frames = fetch.(tasks)
#     return JuliaPIV.main(chunk_frames...)
# end

# @show Threads.nthreads()
# thread_it()

function process_chunk(chunk, start_idx)
    for (i, (frame1, frame2)) in enumerate(chunk)
        fig = JuliaPIV.main(frame1, frame2)
        savefig(fig, "$(start_idx + i).png")
    end
end

function thread_it()
    io = VideoIO.open("tests/test_data/IMG_6566.avi")
    f = VideoIO.openvideo(io)
    pairs = vid(f)

    println("Got pairs, starting threads")

    num_threads = Threads.nthreads()
    chunk_size = ceil(Int, length(pairs) / num_threads)
    chunks = [
            pairs[i:min(i + chunk_size - 1, end)] 
            for i in 1:chunk_size:length(pairs)
        ]

    # Launch threads to process each chunk
    tasks = [
            Threads.@spawn process_chunk(chunk, start_idx) 
            for (chunk, start_idx) in zip(chunks, 1:chunk_size:length(pairs))
        ]
    fetch.(tasks)
end

thread_it()