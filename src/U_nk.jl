

function create_matrix(n, k)
    if n == k
        return 1.0 / (2n + 1)
    elseif (n + k) % 2 == 0
        return 0.0
    else
        A_nk = 1
        return (2 / π) * ((sin(π / (2n)) * cos(π / (2k)) * A_nk) / (n * (n + 1) - k * (k + 1)) - (sin(π / (2k)) * cos(π / (2n)) * A_nk) / (n * (n + 1) - k * (k + 1)))
    end
end

function U(n::Int, k::Int)
    U = zeros(Float64, n, k)
    for i in 1:n
        for j in 1:k 
            U[i, j] = create_matrix(i, j)
        end
    end
    return U
end



