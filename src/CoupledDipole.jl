module CoupledDipole

greet() = print("Hello World!")

end # module

ε0=8.85e-3
include("singleparticle.jl")
include("array_factor.jl")
include("fresnel.jl")
function α_to_array(α,S)
	return 1 ./((1 ./α)-S)
end
function modified_r_s(α,pitch,n1,n2,θ1,ε0=ε0)
	rho=pitch^-2
	θ2=real.(asin.(n1./n2.*sin.(θ1)))-1im*abs.(imag.(asin.(n1./n2*sin.(θ1))))
	denom=n1.*cos.(θ1).+n2.*cos.(θ2).-1im*k.*α*rho/ε0
	r=(n1.*cos.(θ1).-n2.*cos.(θ2).+1im*k.*α*rho/ε0)./denom
	t=2n1.*cos.(θ1)./denom
	return r,t
end
function modified_r_p(α,pitch,n1,n2,θ1,ε0=ε0)
	rho=pitch^-2
	θ2=real.(asin.(n1./n2.*sin.(θ1)))-1im*abs.(imag.(asin.(n1./n2*sin.(θ1))))
	denom=n1.*cos.(θ2).+n2.*cos.(θ1).-1im*k.*α*rho/ε0.*cos.(θ2).*cos.(θ1)
	r=(n1.*cos.(θ2).-n2.*cos.(θ1).+1im*k.*α*rho/ε0*cos.(θ1).*cos.(θ2))./denom
	t=2n1.*cos.(θ1)./denom
	return r,t
end
