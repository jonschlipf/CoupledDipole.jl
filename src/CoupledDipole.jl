module CoupledDipole

greet() = print("Hello World!")

end # module

ε0=8.85e-3 #C/(V*nm)
include("singleparticle.jl")
include("array_factor.jl")
include("fresnel.jl")
include("modified_fresnel.jl")
include("crosssection.jl")
