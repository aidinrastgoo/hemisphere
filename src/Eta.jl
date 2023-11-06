function Eta(k::Integer, permittivity1::Float64 , permittivity2::Float64)
if k % 2==0
    return permittivity2/permittivity1
    else
        return 1
    end

end