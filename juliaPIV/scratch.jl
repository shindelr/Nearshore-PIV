function pad_matrix(matrix::Matrix{T}, pad_size::Tuple{Int,Int}) where T
    padded_matrix = zeros(
                    T, 
                    size(matrix, 1) + 2*pad_size[1], 
                    size(matrix, 2) + 2*pad_size[2]
                    )
    
    padded_matrix[
        (pad_size[1]+1):(end-pad_size[1]), 
        (pad_size[2]+1):(end-pad_size[2])
        ] = matrix

    return padded_matrix
end

function pad(og_matrix, pad_rows, pad_cols)
    og_m_rows = size(og_matrix, 1)
    og_m_cols = size(og_matrix, 2)
    m_zeros = zeros(
        eltype(og_matrix), 
        og_m_rows + 2*pad_rows, 
        og_m_cols + 2*pad_cols
        )
    m_zeros[(pad_rows + 1):(end - pad_rows), (pad_cols + 1):(end - pad_cols)] = og_matrix
    return m_zeros 
end

# ------ TEST ZONE ------

a = [16 2 3 13;
5 11 10 8;  
9 7 6 12;
4 14 15 1]
b = [1 2 3; 
4 5 6; 
7 8 9]
ma, na = size(a)
mb, nb = size(b)
mf = Base._nextpow2(ma + mb)
nf = Base._nextpow2(na + nb)
c = pad(b, mf, nf)
println(c)
# ------ TEST ZONE ------