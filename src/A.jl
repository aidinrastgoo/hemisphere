using SpecialFunctions
include("M.jl")  # Include the M.jl file for the U matrix and η_k function

function calculate_A(Amplitude, ε_r1, ε_r2, k_max)
    U = generate_U_matrix(1, k_max)  # Generate U matrix for n=1
    A = zeros(Float64, k_max)
    
    for k in 1:k_max
        η_k_val = η_k(k, ε_r1, ε_r2)
        term1 = η_k_val * k * ε_r1
        term2 = η_k_val - (-1)^(1 + k) * k * ε_r2 - (-1)^(1 + k)
        A[k] = Amplitude * (term1 - term2) * U[1, k]
    end
    return A
end


