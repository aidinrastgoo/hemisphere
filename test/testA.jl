Amplitude = rand()
ε_r1 = ε_r2 = 1.0
k_max =n_max= 11# rand(10:342)
radius = 1.0
B_h = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
A_h = calculate_A(Amplitude, ε_r1, ε_r2, k_max)

@test sum(B_h) + sum(A_h) == 0  
D = calculate_D(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)'
C = calculate_C(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)'
S = calculate_in(Amplitude, radius, ε_r1, n_max, k_max)
@test isapprox(D, C, atol=1e-10)
@test isapprox(S,D)
