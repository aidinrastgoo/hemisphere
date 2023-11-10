include("A.jl")
include("M.jl")
function calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    M = calculate_M(radius, ε_r1, ε_r2, n_max, k_max)
    A = calculate_A(Amplitude, ε_r1, ε_r2, k_max)
    B = inv(M)*A
    return @show B
end