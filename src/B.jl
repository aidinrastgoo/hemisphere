include("A.jl")
include("M.jl")
function calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    M_h = calculate_M(radius, ε_r1, ε_r2, n_max, k_max)
    A_h = calculate_A(Amplitude, ε_r1, ε_r2, k_max)
    
    # Transpose A_h to make it a column vector
    B_h = A_h'/M_h
    return B_h
end

