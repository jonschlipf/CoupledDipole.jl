#https://pubs.acs.org/doi/10.1021/nn102166t
#Gold, Platinum, and Aluminum Nanodisk Plasmons: Material Independence, Subradiance, and Damping Mechanisms
#Bohren-Huffman

"""
    α0(a,b,ε,εout=1,ε0=ε0)
Computes the polarizability of a spheroid particle
# Arguments
* `a` : x and y radius
* `b` : z radius
* `ε` : inside relative permittivity
* `εout` : optional outside relative permittivity
* `ε0` : optional permittivity of free space
# Output
* `α0`: polarizability of the particle
"""
function α0(a,b,ε,εout=1,ε0=ε0)
	if a==b #spherical case
		L=1/3
	else #ellipsoid, L is depolarization factor
		esquared=1-b^2*a^-2
		g=√((1-esquared)/esquared)
		L=g/(2*esquared)*(π/2-atan(g))-.5*g^2
	end
	return ε0*4π*a^2*b/3*(ε .-εout)./(εout .+L*(ε.-εout))
end
"""
    αMLWA(α0,k,a,ε0=ε0)
    αMLWA(a,b,ε,k,εout=1,ε0=ε0)
Computes the polarizability of a spheroid particle in the modified long wavelength approximation
# Arguments
* `α0` : non-MLWA polarizability
* `k` : wavenumber of propagation outside the particle
* `a` : x and y radius
* `b` : z radius
* `ε` : inside relative permittivity
* `εout` : optional outside relative permittivity
* `ε0` : optional permittivity of free space
# Output
* `αMLWA`: polarizability of the particle
"""
function αMLWA(α0,k,a,ε0=ε0)
	α0=α0/ε0 #convert α to different convention
	return ε0*α0./(1 .-k.^2 ./(4π*a).*α0-1im*k.^3 ./(6π).*α0)
end
function αMLWA(a::Real,b::Real,ε,k,εout=1,ε0=ε0)
	return αMLWA(α0(a,b,ε,εout,ε0),k,a,ε0)
end
