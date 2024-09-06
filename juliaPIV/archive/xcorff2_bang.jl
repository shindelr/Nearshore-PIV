function xcorrf2!(R::Matrix{T}, A::Matrix{T}, B::Matrix{T}, plan!, iplan!, M::Int32) where {T}
    M = N
    # Unpack size() return tuple into appropriate variables
    mb, nb = size(B)

    # Reverse conjugate
    B = conj(B[mb:-1:1, nb:-1:1])

    pad_matrix_A = pad_for_xcorr(A[1:M, 1:N])
    pad_matrix_B = pad_for_xcorr(B[1:M, 1:N])

    # Transfer data from og matrix to optimized sized ones
    pad_matrix_A[1:size(A,1), 1:size(A,2)] = A[1:size(A,1), 1:size(A,2)]
    pad_matrix_B[1:size(B,1), 1:size(B,2)] = B[1:size(B,1), 1:size(B,2)]

    R = real(iplan! * ((plan! * pad_matrix_b) .* (plan! * pad_matrix_a))
            )[1:ma + mb - 1, 1:na + nb - 1]
end