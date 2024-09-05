function xcorrf2!(R::Matrix{T}, A::Matrix{T}, B::Matrix{T}, plan!, iplan!, 
    pad_matrix_a::Matrix{ComplexF32}, pad_matrix_b::Matrix{ComplexF32}) where {T}

# Unpack size() return tuple into appropriate variables
mb, nb = size(B)

# Reverse conjugate
B = conj(B[mb:-1:1, nb:-1:1])

# Transfer data from og matrix to optimized sized ones
pad_matrix_a[1:size(A,1), 1:size(A,2)] = A[1:size(A,1), 1:size(A,2)]
pad_matrix_b[1:size(B,1), 1:size(B,2)] = B[1:size(B,1), 1:size(B,2)]

# FFT in place
plan! * pad_matrix_a
plan! * pad_matrix_b

# Mult transforms and invert
pad_matrix_a .*= pad_matrix_b    
iplan! * pad_matrix_a

# Make all real and store in R
R .= real.(pad_matrix_a)

# Reset pads
fill!(pad_matrix_a, 0)
fill!(pad_matrix_b, 0)
end