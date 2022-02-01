__precompile__()

module LA

using Random
using LinearAlgebra

export magnitude
export normalize_minmax
export outer
export rms
export vec_angle

include("magnitude.jl")
include("normalize_minmax.jl")
include("outer.jl")
include("rms.jl")
include("vec_angle.jl")

end
