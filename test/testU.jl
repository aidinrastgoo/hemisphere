Amplitude = 1.0
ε_r1 = 1.0
ε_r2 = 1.0
n_max = 6
k_max = 6
radius = 1.0  

A = calculate_A(Amplitude, ε_r1, ε_r2, k_max)
M = calculate_M(radius, ε_r1, ε_r2, n_max, k_max)  # Calculate the M matrix
B = A \ M  





@test sum(B)==3.757952410879924