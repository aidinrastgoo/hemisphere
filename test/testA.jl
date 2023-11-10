Amplitude = 99999999999999999999999.0
ε_r1 = 1.0
ε_r2 = 1.0
k_max = 342
A = calculate_A(Amplitude, ε_r1, ε_r2, k_max)
@test sum(A) == 0