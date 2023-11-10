module hemisphere

import Pkg
Pkg.add("SpecialFunctions")

include("A.jl")
include("M.jl")
include("B.jl")


export calculate_A
export A_kn
export create_matrix
export η_k
export generate_U_matrix
export calculate_M
export calculate_B

end # module hemisphere
