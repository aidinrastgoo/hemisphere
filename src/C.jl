include("M.jl")
include("B.jl")
include("D.jl")
function calculate_C(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D = calculate_D(Amplitude, radius, ε_r1, ε_r2, n_max, k_max) 
    η = calculate_η(k_max, ε_r1, ε_r2)
    return η.* D
end


function calculate_C_My(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D = calculate_D_My(Amplitude, radius, ε_r1, ε_r2, n_max, k_max) 
    η = calculate_η(k_max, ε_r1, ε_r2)
    return η.* D
end

function calculate_C_test(Amplitude, radius, ε_r1, ε_r2, n_max, k_max) # it is calculate only by one Potential_eq. Thus is not valid.
     B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
     U = generate_U_matrix(n_max, k_max)
     C = zeros(Float64, n_max)
     term_1 = zeros(Float64, n_max)
     term_2 = zeros(Float64, k_max)
     term_3 = zeros(Float64, n_max, k_max)
       for n in 0:n_max-1  
            for k in 0:k_max-1
                term_1[n+1] = radius^(-n-1) * B[n+1] * U[k+1, n+1]
                term_2[k+1] = radius * Amplitude * U[2, k+1]
                term_3[k+1, n+1] = U[k+1, n+1] * (radius^k)
            end
        end
        return inv(term_3) * (term_1 - term_2) 
end

