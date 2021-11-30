

#Snell's law
function fresnel_angles(n,θ0)
	θ=real.(asin.(n[1]./n*sin(θ0)))-1im*abs.(imag.(asin.(n[1]./n*sin(θ0))))
	return θ
end

function fresnel_coefs_s(n,theta)
	#s polarization
	denom=n[1:end-1].*cos.(theta[1:end-1])+n[2:end].*cos.(theta[2:end])
	r=(n[1:end-1].*cos.(theta[1:end-1])-n[2:end].*cos.(theta[2:end]))./denom
	t=2n[1:end-1].*cos.(theta[1:end-1])./denom
	rb=(n[2:end].*cos.(theta[2:end])-n[1:end-1].*cos.(theta[1:end-1]))./denom
	tb=2n[2:end].*cos.(theta[2:end])./denom
	return r,t,rb,tb
end
#p polarization
function fresnel_coefs_p(n,theta)
	denom=n[1:end-1].*cos.(theta[2:end])+n[2:end].*cos.(theta[1:end-1])
	r=(n[1:end-1].*cos.(theta[2:end])-n[2:end].*cos.(theta[1:end-1]))./denom
	rb=(n[2:end].*cos.(theta[1:end-1])-n[1:end-1].*cos.(theta[2:end]))./denom
	t=2n[1:end-1].*cos.(theta[1:end-1])./denom
	tb=2n[2:end].*cos.(theta[2:end])./denom
	return r,t,rb,tb
end
#phase shift when passing through interface
function fresnel_phase(n,theta,h,k)
	delta=h.*k.*n[2:end-1].*cos.(theta[2:end-1])
	return delta
end
#t-matrix method
function fresnel_transfermatrix(r,t,delta,rb=r,tb=t)
	#preset
	M=[1 0;0 1]
	#first interfaces and propagation
	for i=1:length(r)-1
		#M i generalized for the r12!=-r21 case
		M=M*(1/t[i])*[1 -rb[i];r[i] t[i]*tb[i]-r[i]*rb[i]]*
		[exp(-1im*delta[i]) 0;0 exp(1im*delta[i])]
	end
	#last interface
	M=M*(1/t[end])*[1 -rb[end];r[end] t[end]*tb[end]-r[end]*rb[end]]
	return M
end
#get total reflection and transmission coefficient
function fresnel_total(r,t,delta,rb=r,tb=t)
	M=fresnel_transfermatrix(r,t,delta,rb,tb)
	rtot=M[2,1]/M[1,1]
	ttot=1/M[1,1]
	if length(r)==1
		rtot=r[1]
		ttot=t[1]
	end
	return rtot,ttot
end
#total reflection and transmission
function fresnel_intensity(r,t,delta,n,θ,rb=-r,tb=t)
	r,t=fresnel_total(r,t,delta,rb,tb)
	R=abs.(r).^2
	#normalize T by wave vector z component
	T=abs.(t).^2*real(n[end]*cos(θ[end]))/real(n[1]*cos(θ[1]))
	return R,T
end
