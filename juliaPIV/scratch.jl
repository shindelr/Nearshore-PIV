A = [1 2 3 4 5;
16 2 3 13 2;
5 11 10 8 9;  
9 7 6 12 4;
4 14 15 1 7
]

B = [3 4 5 6 7;
8 9 10 11 12;
13 14 15 16 17;
18 19 20 21 22;
23 24 25 26 27
]

pivwin = 16
log2pivwin = log2(pivwin)
if log2pivwin - round(log2pivwin) != 0
    error("pivwin must be factor of 2")
end
pass_sizes = 2 .^ collect(6:-1:log2pivwin)
push!(pass_sizes, pass_sizes[end]) # Duplicate final element
pass_sizes = [pass_sizes pass_sizes] # Convert vector to matrix
dt = 1; overlap = 0.5; validvec = 3
A = convert(Matrix{Float64}, A)
B = convert(Matrix{Float64}, B)
sy, sx = size(A)
iter = size(pass_sizes, 1)

data_dim_1 = floor(Int64, (sy/(pass_sizes[1,1] * 1-overlap)))
data_dim_2 = floor(Int64, (sx/(pass_sizes[1,2]*(1-overlap)))) 
datax = zeros(eltype(A), (data_dim_1, data_dim_2))
datay = zeros(eltype(A), (data_dim_1, data_dim_2))
