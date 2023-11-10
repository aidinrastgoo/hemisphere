using SpecialFunctions
include("M.jl")  # Include the M.jl file for the U matrix and η_k function

function calculate_A(Amplitude, ε_r1, ε_r2, k_max)
    U = generate_U_matrix( 10 , k_max)  # Generate U matrix . Hier ist wert n_max=10 .In Generell soll der Wert  n > 2 , weil der zweite Reihe von Matrix verwendet wird, was ist n=1
    A = zeros(Float64, k_max)

    for k in 0:k_max-1
        η_k_val = η_k(k, ε_r1, ε_r2)
        term1 = η_k_val * k * ε_r1
        term2 = ((-1)^(1 + k)) * k * ε_r2 
        term3 = (-1)^(1 + k)
        A[k+1] = Amplitude * (term1 - η_k_val + term2 - term3) * U[2, k+1]
    end
    return @show A
end


