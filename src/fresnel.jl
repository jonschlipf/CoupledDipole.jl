


function fresnel_angles(n,θ0)
	θ=real.(asin.(n[1]./n*sin(θ0)))-1im*abs.(imag.(asin.(n[1]./n*sin(θ0))))
	return θ
end

function fresnel_coefs_s(n,theta)
	r=(n[1:end-1].*cos.(theta[1:end-1])-n[2:end].*cos.(theta[2:end]))./(n[1:end-1].*cos.(theta[1:end-1])+n[2:end].*cos.(theta[2:end]))
	t=2n[1:end-1].*cos.(theta[1:end-1])./(n[1:end-1].*cos.(theta[1:end-1])+n[2:end].*cos.(theta[2:end]))
	return r,t
end
function fresnel_coefs_p(n,theta)
	r=(n[1:end-1].*cos.(theta[2:end])-n[2:end].*cos.(theta[1:end-1]))./(n[1:end-1].*cos.(theta[2:end])+n[2:end].*cos.(theta[1:end-1]))
	t=2n[1:end-1].*cos.(theta[1:end-1])./(n[1:end-1].*cos.(theta[2:end])+n[2:end].*cos.(theta[1:end-1]))
	return r,t
end
function fresnel_phase(n,theta,h,k)
	delta=h.*k.*n[2:end-1].*cos.(theta[2:end-1])
	return delta
end
function fresnel_transfermatrix(r,t,delta)
	M=[1 0;0 1]
	for i=1:length(r)-1
		M=M*(1/t[i])*[1 r[i];r[i] 1]*
		[exp(-1im*delta[i]) 0;0 exp(1im*delta[i])]
	end
	M=M*(1/t[end])*[1 r[end];r[end] 1]
	return M
end
function fresnel_total(r,t,delta)
	M=fresnel_transfermatrix(r,t,delta)
	rtot=M[2,1]/M[1,1]
	ttot=1/M[1,1]
	if length(r)==1
		rtot=r[1]
		ttot=t[1]
	end
	return rtot,ttot
end
function fresnel_intensity(r,t,delta,n,θ)
	r,t=fresnel_total(r,t,delta)
	R=abs.(r).^2
	T=abs.(t).^2*real(n[end]*cos(θ[end]))/real(n[1]*cos(θ[1]))
	return R,T
end
