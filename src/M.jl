using SpecialFunctions

function A_kn(k, n)
    numerator_k = gamma(k/2 + 1) / gamma(k/2 + 1/2)
    denominator_n = gamma(n/2 + 1/2) / gamma(n/2 + 1)
    return numerator_k * denominator_n
end

function create_matrix(n, k)
    if n == k
        return 1.0 / (2n + 1)
    elseif (n + k) % 2 == 0
        return 0.0
    elseif n == 0 && k >= 1
        t1 = cos(π / (2k)) * (1 / A_kn(k, n))
        t2 = -k * (k + 1)
        t3 = sin(π / (2k)) * A_kn(k, n)
        return (2 / π) * ((t1 / t2) - (t3 / t2))
    elseif k == 0 && n >= 1
        t11 = sin(π / (2n)) * (1 / A_kn(k, n))
        t12 = n * (n + 1)
        t13 = cos(π / (2n)) * A_kn(k, n)
        return (2 / π) * ((t11 / t12) - (t13 / t12))
    else
        t21 = (sin(π / (2n)) * cos(π / (2k))) / A_kn(k, n)
        t22 = n * (n + 1) - k * (k + 1)
        t23 = sin(π / (2k)) * cos(π / (2n)) * A_kn(k, n)
        return (2 / π) * ((t21 / t22) - (t23 / t22))
    end
end



function generate_U_matrix(n_max, k_max)
    U = zeros(Float64, n_max, k_max)
    for n in 0:n_max-1
        for k in 0:k_max-1
            U[n+1, k+1] = create_matrix(n, k)
        end
    end
    return U
end



function η_k(k, ε_r1, ε_r2)
    if k % 2 == 0
        return 1.0
    else
        return ε_r2 / ε_r1
    end
end

function calculate_M(radius, ε_r1, ε_r2, n_max, k_max)
    U = generate_U_matrix(n_max, k_max)
    M = zeros(Float64, n_max, k_max)
    
    for n in 0:n_max-1
        for k in 0:k_max-1
            η_k_val = η_k(k, ε_r1, ε_r2)
            term1 = η_k_val * (n + 1)
            term2 = η_k_val * k * ε_r1
            term3 = ((-1)^(n + k)) * (n + 1)
            term4 = ((-1)^(n + k)) * k * ε_r2
            M[n+1, k+1] = (radius^(-(n + 2))) * (term1 + term2 + term3 + term4) * U[n+1, k+1]
        end
    end
    return @show M
end