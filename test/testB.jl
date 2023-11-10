Amplitude = 99999999999999999999999.0
ε_r1 = 1.0
ε_r2 = 1.0
n_max = 342 
k_max = 342
radius = 7.9  
B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
@test sum(B) == 0