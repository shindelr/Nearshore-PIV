using CUDA
CUDA.allowscalar(false)
using BenchmarkTools

N = 2050
a = CUDA.rand(N)
b = CUDA.rand(N)
c = similar(a)

function my_add!(c::AbstractArray, a::AbstractArray, b::AbstractArray)
    for i in eachindex(c)
        c[i] = a[i] + b[i]
    end
    nothing
end

function _my_add_kernel(c, a, b)
    i = CUDA.threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i <= length(c)
        c[i] = a[i] + b[i]
    end

    return nothing
end

function my_add!(c::CuArray, a::CuArray, b::CuArray)
    @cuda blocks=cld(length(c), 1024) threads=1024 _my_add_kernel(c, a, b)
    nothing
end    

my_add!(c, a, b)

isapprox(Array(c), Array(a) .+ Array(b))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

function _est_pi_kernel(global_count)
    # Create som shared mem for each thread in the block
    hits = CUDA.@cuStaticSharedMem(UInt16, 256)

    # Throw a dart and calculate whether it has hit
    x = rand(Float32) * 2 - 1
    y = rand(Float32) * 2 - 1
    is_hit = (x^2+y^2<=1)

    idx = threadIdx().x
    # Record hit in shared mem
    hits[idx] = UInt16(is_hit)

    # Perform reduction on shared memory
    step_size = Int32(div(256, 2))
    while step_size != 0
        CUDA.sync_threads()
        if idx <= step_size
            hits[idx] += hits[idx + step_size]
        end
        step_size = Int32(div(step_size, 2))
    end
    total = hits[1]  # By the end of the reduction, this is only one elem

    # add the count from the block into the global count
    if idx == 1
        # Atomic makes it happen sequentially between blocks, eliminating race cond
        CUDA.@atomic global_count[] += hits[1]
        # Emtpy brackets indexes the first element
    end

    return nothing
end

function est_pi(n)
    # calculate number of threads blocks
    threads = 256
    blocks = cld(n, threads)

    # create some memory to store the count
    total_count = CUDA.zeros(UInt32, 1)
    # run kernel
    @cuda blocks=blocks threads=threads _est_pi_kernel(total_count)
    # Transfer the finished count from the GPU
    count = UInt32(0)
    # In this case, just getting a single value out of the array, so no need to
    # create a new array to transfer to like Array(total_count). Instead, just
    # allow scalar idx and grab the first elem
    CUDA.@allowscalar begin
        count = total_count[]
    end

    # Free the memory from the GPU if you know you won't use it again
    CUDA.unsafe_free!(total_count) 

    return 4 * count / (blocks * threads)
end


# Check for type stability:
# @device_code_warntype est_pi(2^20)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

function matrix_add!()
    d_a = CUDA.rand(Float32, 100, 100)
    d_b = CUDA.rand(Float32, 100, 100)

    # Mem to store the matrix addition
    d_c = CUDA.zeros(Float32, 100, 100)

    threads = (32, 32)
    blocks = (cld(size(d_a, 1), threads[1]), cld(size(d_a, 2), threads[2]))
    @cuda blocks=blocks threads=threads _matrix_add(d_a, d_b, d_c)

    isapprox(Array(d_c), Array(d_a) .+ Array(d_b))
    return nothing
end

function _matrix_add(d_a, d_b, d_c)
    i = Int32(threadIdx().x + (blockIdx().x - 1) * blockDim().x)
    j = Int32(threadIdx().y + (blockIdx().y - 1) * blockDim().y)

    if i <= size(d_c, 1) && j <= size(d_c, 2)
        d_c[i, j] = d_a[i, j] + d_b[i, j]
    end

    return nothing
end

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

function scale_picture!(n, m)
    d_pic_in = CUDA.rand(Float32, n, m)
    d_pic_out = CUDA.zeros(Float32, n, m)

    threads = (32, 32)  # Threads per block
    dim_grid = (Int32(ceil(n/threads[1])), Int32(ceil(m/threads[2])))  # Blocks in grid
    @cuda threads=threads blocks=dim_grid _scale_picture(d_pic_out, d_pic_in, n, m)

    # isapprox(Array(d_pic_out), Array(d_pic_in) .* 2)
    return nothing
end

function _scale_picture(d_pic_out::T, d_pic_in::T, n::F, m::F) where {T, F}

    # Calculate row, col
    # Always minus 1 from blockIdx bc julia is 1-indexed
    row = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    col = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    # Scale each pixel by 2 if w/in boundaries
    # Store in the out picture
    # <= bc julia is 1-indexed
    if row <= m && col <= n
        # -1 from row bc julia is 1-indexed
        @inbounds d_pic_out[(row - 1) * n + col] = d_pic_in[(row - 1) * n + col] * 2
    end

    return nothing
end

function tiled_matrix_mul!(n, m)
    matrix_a = CUDA.rand(Float32, n, m)
    matrix_b = CUDA.rand(Float32, n, m)
    matrix_c = CUDA.zeros(Float32, n, m)

    threads = (2, 2)
    dim_grid = (Int32(ceil(n/threads[1])), Int32(ceil(m/threads[2])))
    @cuda threads=threads blocks=dim_grid shmem=sizeof(matrix_a) _matrix_mul(matrix_c, matrix_a, matrix_b, n)

    return isapprox(Array(matrix_c), Array(matrix_a) .* Array(matrix_b))
end

function _matrix_mul(matrix_c, matrix_a, matrix_b, global_width)

    # Details about this particular thread
    bx = blockIdx().x - 1
    by = blockIdx().y - 1
    tx = threadIdx().x
    ty = threadIdx().y

    # Global indices of the whole matrix
    i = blockDim().y * by + ty
    j = blockDim().x * bx + tx

    # Shared memory tile allocation 
    shared_a = CuDynamicSharedArray(Float32, 4)
    shared_b = CuDynamicSharedArray(Float32, 4)

    num_tiles = global_width / 2
    val = 0f0
    phase = 0
    while phase < num_tiles
        @cushow ((i) * global_width + phase * 2 + tx)
        if i < global_width && j < global_width
            shared_a[ty, tx] = matrix_a[(i) * global_width + phase * 2 + tx]
            shared_b[ty, tx] = matrix_b[(phase * global_width + ty) * global_width + j]
        end
        
        phase += 1
    end

    return nothing
end