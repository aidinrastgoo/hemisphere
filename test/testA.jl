Amplitude = rand()
ε_r1 = ε_r2 = 1.0
k_max =n_max= rand(10:342)
radius = rand()
B_h = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
A_h = calculate_A(Amplitude, ε_r1, ε_r2, k_max)

@test sum(B_h) + sum(A_h) == 0  