module hemisphere

import Pkg
Pkg.add("SpecialFunctions")

include("A.jl")
include("M.jl")
include("B.jl")
include("Sphere.jl")


export calculate_A
export A_kn
export create_matrix
export η_k
export generate_U_matrix
export calculate_M
export calculate_B
export calculate_MSph
export calculate_ASph
export calculate_BSphere

end # module hemisphere
