include("M.jl")
include("B.jl")
function calculate_LHP(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)  #lower Hiemspere Potential-eq. 
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    term1 = zeros(Float64, n_max)
    term2 = zeros(Float64, k_max)
    X = zeros(Float64, n_max, k_max)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1[n+1] = (-1)^(k + n) * radius^(-n-1) * B[n+1] * U[k+1, n+1]
            term2[k+1] = (-1)^(k + 1) * radius * Amplitude * U[2, k+1]
            X[k+1, n+1] = U[k+1, n+1] * (radius^k)*(-1)^(k+n)
        end
    end
    return inv(X) * (term1 - term2)
end

using LinearAlgebra   # it is not possible to invers Matrix of X , therefore is used Pseudoinverse. Thus we need the package of LinearAlgebra

function calculate_LHE(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)   #lower Hiemspere Electric-Field-eq.
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    term1_LHE = zeros(Float64, n_max)
    term2_LHE = zeros(Float64, k_max)
    X_LHE = zeros(Float64, n_max, k_max)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1_LHE[n+1] = (-n-1)*(-1)^(k + n) * radius^(-n-2) * B[n+1] * U[k+1, n+1]
            term2_LHE[k+1] = (-1)^(k + 1) * Amplitude * U[2, k+1]
            X_LHE[k+1, n+1] += U[k+1, n+1]  * k * ε_r2* (radius^(k-1))*((-1)^(k+n))
        end
    end
    return pinv(X_LHE) * (term1_LHE - term2_LHE)   #Pseudoinverse
end

function calculate_D_test(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_1 = calculate_LHP(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_2 = calculate_LHE(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    return D_1 + D_2
end 



function calculate_UHP(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)  #Upper Hiemspere Potential-eq. 
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    η = calculate_η(k_max, ε_r1, ε_r2)
    term1 = zeros(Float64, n_max)
    term2 = zeros(Float64, k_max)
    X = zeros(Float64, n_max, k_max)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1[n+1] = radius^(-n-1) * B[n+1] * U[k+1, n+1]
            term2[k+1] = radius * Amplitude * U[2, k+1]
            X[k+1, n+1] = η[k+1] * U[k+1, n+1] * (radius^k)
        end
    end
    return inv(X) * (term1 - term2)
end

function calculate_UHE(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)   #upper Hiemspere Electric-Field-eq.
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    term1_UHE = zeros(Float64, n_max)
    term2_UHE = zeros(Float64, k_max)
    X_UHE = zeros(Float64, n_max, k_max)
    η = calculate_η(k_max, ε_r1, ε_r2)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1_UHE[n+1] = (-n-1) * radius^(-n-2) * B[n+1] * U[k+1, n+1]
            term2_UHE[k+1] = Amplitude * U[2, k+1]
            X_UHE[k+1, n+1] += U[k+1, n+1]  * k * η[k+1] * ε_r1 * (radius^(k-1))
        end
    end
    return pinv(X_UHE) * (term1_UHE - term2_UHE)   #Pseudoinverse
end

function calculate_D_test2(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_1 = calculate_UHP(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_2 = calculate_UHE(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_3 = calculate_LHP(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    D_4 = calculate_LHE(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    return D_1 - D_2 + D_3 - D_4
end 


function calculate_D_test3(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)   #Lower Hiemspere Potantial-eq + Electric-Field-eq.
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    term1_LHE = zeros(Float64, n_max)
    term2_LHE = zeros(Float64, n_max, k_max)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1_LHE[n+1] =((-1)^(k+1)) * (-n-2)  * Amplitude * U[2, k+1]
            term2_LHE[k+1, n+1] = ((-1)^(k+1))* (radius^(k-1)) *U[k+1, n+1] * ((k * ε_r2) + (n+1))
        end
    end
    return inv(term2_LHE) * (term1_LHE)  
end


function calculate_D(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)   # Papier
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U =  generate_U_matrix(n_max, k_max)
    term1_LHE = zeros(Float64, n_max)
    term2_LHE = zeros(Float64, n_max, k_max)
     η = calculate_η(k_max, ε_r1, ε_r2)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1_LHE[n+1] = -Amplitude * radius * ( (1/(n+1)) + 1 + (((-1)^(1+n)) / (n+1)) + ((-1)^(1+n)) ) * U[2, n+1]
            term2_LHE[k+1, n+1] = radius^(k) * (η[k+1] + ((η[k+1] * k * ε_r1 )/(n+1)) +  ((-1)^(n+k)) +  ( (((-1)^(n+k)) * k * ε_r2)  /  (n+1))  ) * U[k+1, n+1] 
        end
    end
     term1_LHE 
     term2_LHE 
    return inv(term2_LHE) * (term1_LHE)  
end




function calculate_D_My(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)   # My calculation
    B = calculate_B(Amplitude, radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    term1_LHE = zeros(Float64, n_max)
    term2_LHE = zeros(Float64, n_max, k_max)
    η = calculate_η(k_max, ε_r1, ε_r2)
    for n in 0:n_max-1  
        for k in 0:k_max-1
            term1_LHE[n+1] = -Amplitude * ( ((-1)^(n+1)) + ((n+1) * (-1)^(n+1)) + ((n+1) * (-1)^(n+n)) + ((-1)^(n+n)) ) * U[2, n+1]
            term2_LHE[k+1, n+1] = (radius^(k-1)) * ((-1)^(k+n)) * ( (k * ε_r2) + (n+1) + ((n+1) * η[k+1]) + (ε_r1 * k * η[k+1])    ) * U[k+1, n+1] 
        end
    end
     term1_LHE 
     term2_LHE 
    return inv(term2_LHE) * (term1_LHE)  

end
