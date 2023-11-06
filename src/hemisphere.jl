module hemisphere

import Pkg
Pkg.add("SpecialFunctions")

include("U_nk.jl")

export A_kn
export create_matrix
export U

end # module hemisphere
