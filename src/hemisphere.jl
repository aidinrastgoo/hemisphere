module hemisphere
using LinearAlgebra
import Pkg
Pkg.add("SpecialFunctions")

include("A.jl")
include("M.jl")
include("B.jl")
include("Sphere.jl")
include("D.jl")
include("C.jl")


export calculate_A
export A_kn
export create_matrix
export calculate_η
export generate_U_matrix
export calculate_M
export calculate_B
export calculate_MSph
export calculate_ASph
export calculate_BSphere
export calculate_D
export calculate_C
export calculate_D_My
export calculate_C_My
export calculate_in

end # module hemisphere

