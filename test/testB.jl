Amplitude = rand()
ε_r1 = rand()
ε_r2 = ε_r1
n_max = k_max = rand(10:342)
radius = rand()

B_h = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
B_s =calculate_BSphere(Amplitude, radius, ε_r1, n_max, k_max)


function euclidean_distance(vec1, vec2)
    return sqrt(sum((vec1 .- vec2) .^ 2))
end

distance = euclidean_distance(B_h, B_s)
@show distance

@test distance < 1.0e-16
