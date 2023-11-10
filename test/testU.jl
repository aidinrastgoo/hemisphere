Amplitude = 1.0
ε_r1 = 1.0
ε_r2 = 1.0
n_max = 10
k_max = 10
radius = 1.0  
U = generate_U_matrix(n_max, k_max)
A = calculate_A(Amplitude, ε_r1, ε_r2, k_max)
M = calculate_M(radius, ε_r1, ε_r2, n_max, k_max)  # Calculate the M matrix
@show B = inv(M)*A  



@test sum(B) == 0