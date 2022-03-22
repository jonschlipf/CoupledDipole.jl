module CoupledDipole

greet() = print("Hello World!")

end # module

ε0=8.85e-3
include("singleparticle.jl")
include("array_factor.jl")
include("fresnel.jl")
#incorporate dipole polarizability into array
function α_to_array(α,S)
	return 1 ./((1 ./α)-S)
end
#modified fresnel coefficients for dipoles on interface
function modified_r_s(α,pitch,n1,n2,θ1,k,ε0=ε0)
    θ1=θ1+0im
	rho=pitch.^-2 #dipole area density
	θ2=real.(asin.(n1./n2.*sin.(θ1)))-1im*abs.(imag.(asin.(n1./n2*sin.(θ1)))) #Snell
	denom=n1.*cos.(θ1).+n2.*cos.(θ2).-1im*k.*α.*rho/ε0 #denominator
	r=(n1.*cos.(θ1).-n2.*cos.(θ2).+1im*k.*α.*rho/ε0)./denom#forward
	rb=(n2.*cos.(θ2).-n1.*cos.(θ1).+1im*k.*α.*rho/ε0)./denom#backward
	t=2n1.*cos.(θ1)./denom #forward t
	tb=2n2.*cos.(θ2)./denom #backward t
	return r,t,rb,tb
end
#same for p polarization
function modified_r_p(α,pitch,n1,n2,θ1,k,ε0=ε0)
    θ1=θ1+0im
	rho=pitch.^-2
	θ2=real.(asin.(n1./n2.*sin.(θ1)))-1im*abs.(imag.(asin.(n1./n2*sin.(θ1))))
	denom=n1.*cos.(θ2).+n2.*cos.(θ1).-1im*k.*α.*rho/ε0.*cos.(θ2).*cos.(θ1)
	r=(n1.*cos.(θ2).-n2.*cos.(θ1).+1im*k.*α.*rho/ε0*cos.(θ1).*cos.(θ2))./denom
	rb=(n2.*cos.(θ1).-n1.*cos.(θ2).+1im*k.*α.*rho/ε0*cos.(θ1).*cos.(θ2))./denom
	t=2n1.*cos.(θ1)./denom
	tb=2n2.*cos.(θ2)./denom
	return r,t,rb,tb
end
