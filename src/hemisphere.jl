module hemisphere

import Pkg
Pkg.add("SpecialFunctions")

include("A.jl")
include("M.jl")


export calculate_A
export A_kn
export create_matrix
export η_k
export generate_U_matrix

end # module hemisphere
