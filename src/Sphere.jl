function create_σ(n, k)
    if n == k
        return 1
    else
        return 0
    end
end

function generate_σ(n_max, k_max)
    σ = zeros(Int, n_max, k_max)
    for n in 1:n_max
        for k in 1:k_max
            σ[n, k] = create_σ(n, k)
        end
    end
    return σ
end

function calculate_MSph(radius, ε_r1, n_max, k_max)
    σ = generate_σ(n_max, k_max)
    M_s = zeros(Float64, n_max, k_max)
    
    for n in 0:n_max-1
        for k in 0:k_max-1
            term1 = - n * ε_r1 
            term2 = -(n + 1)
            M_s[n+1, k+1] = (radius^(-(n + 2))) * (term1 + term2) * σ[n+1, k+1]
        end
    end
    
    return @show M_s
end

function calculate_ASph(Amplitude, ε_r1, k_max)
    σ = generate_σ(k_max, k_max)
    A_s = zeros(Float64, 1, k_max)  # Adjusted dimensions to be a row vector
    for k in 0:k_max-1
        A_s[1, k+1] = Amplitude * (1-(k * ε_r1)) * σ[2, k+1]
    end
    
    return A_s
end
