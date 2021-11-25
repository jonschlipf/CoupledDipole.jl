#MIE theory
function α0(a,b,ε,εout=1,ε0=ε0)
	if a==b #spherical case
		L=1/3
	else #L is depolarization factor
		e=1-b^2*a^-2
		g=√((1-e^2)*e^-2)
		L=g/(2*e^2)*(π/2-atan(g))-.5*g^2
	end
	return ε0*4π*a^2*b/3*(ε .-εout)./(εout .+L*(ε.-εout))
end
#MLWA approximation
function αMLWA(α0,k,a,ε0=ε0)
	α0=α0/ε0 #convert α to different convention
	return ε0*α0./(1 .-k.^2 ./(4π*a).*α0-1im*k.^3 ./(6π).*α0)
end
#both in one
function αMLWA(a::Real,b::Real,ε,k,εout=1,ε0=ε0)
	return αMLWA(α0(a,b,ε,εout,ε0),k,a,ε0)
end
