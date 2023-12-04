
Amplitude = rand()
ε_r1 = ε_r2 = rand()
n_max = k_max = rand(10:342)
radius = rand()

B_h = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
B_s =calculate_BSphere(Amplitude, radius, ε_r1, n_max, k_max)

diff = sqrt(sum((B_h .- B_s) .^ 2))
#@test  diff < 1.0e-16

 D = calculate_D(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)'
 C = calculate_C(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)'
 S = calculate_in(Amplitude, radius, ε_r1, n_max, k_max)
@test isapprox(D, C , atol=1e-10) || isapprox(D, S , atol=1e-10)