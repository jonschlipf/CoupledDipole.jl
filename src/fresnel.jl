
"""
    fresnel_angles(n,θ0)           
Computes the angles of ray propagation through a multilayer stack according to Snell's law                                  
# Arguments
* `n` : array of refractive indices 
* `θ0` :  angle of ray in first layer                     
# Output
* `θ`: array of angles in all layers                      
"""
function fresnel_angles(n,θ0)
    #complex Snell's law
	θ=real.(asin.(n[1]./n*sin(θ0)))-1im*abs.(imag.(asin.(n[1]./n*sin(θ0))))
	return θ
end

"""
    fresnel_coefs_s(n,theta)           
Computes the fresnel coefficients for a multilayer stack for s/TE polarization
# Arguments
* `n` : array of refractive indices 
* `theta` : array of angles                     
# Output
* `r`: array of fresnel reflection coefficients for forward propagation
* `t`: array of fresnel transmission coefficients for forward propagation
* `rb`: array of fresnel reflection coefficients for backward propagation
* `tb`: array of fresnel transmission coefficients for backward propagation
"""
function fresnel_coefs_s(n,theta)
	#s polarization
    #get denominator
	denom=n[1:end-1].*cos.(theta[1:end-1])+n[2:end].*cos.(theta[2:end])
    #reflection
	r=(n[1:end-1].*cos.(theta[1:end-1])-n[2:end].*cos.(theta[2:end]))./denom
    #transmission
	t=2n[1:end-1].*cos.(theta[1:end-1])./denom
    #backward reflection
	rb=(n[2:end].*cos.(theta[2:end])-n[1:end-1].*cos.(theta[1:end-1]))./denom
    #backward transmission
	tb=2n[2:end].*cos.(theta[2:end])./denom
	return r,t,rb,tb
end
"""
    fresnel_coefs_p(n,theta)           
Computes the fresnel coefficients for a multilayer stack for p/TM polarization
# Arguments
* `n` : array of refractive indices 
* `theta` : array of angles                     
# Output
* `r`: array of fresnel reflection coefficients for forward propagation
* `t`: array of fresnel transmission coefficients for forward propagation
* `rb`: array of fresnel reflection coefficients for backward propagation
* `tb`: array of fresnel transmission coefficients for backward propagation
"""
function fresnel_coefs_p(n,theta)
    #see s polarization
	denom=n[1:end-1].*cos.(theta[2:end])+n[2:end].*cos.(theta[1:end-1])
	r=(n[1:end-1].*cos.(theta[2:end])-n[2:end].*cos.(theta[1:end-1]))./denom
	rb=(n[2:end].*cos.(theta[1:end-1])-n[1:end-1].*cos.(theta[2:end]))./denom
	t=2n[1:end-1].*cos.(theta[1:end-1])./denom
	tb=2n[2:end].*cos.(theta[2:end])./denom
	return r,t,rb,tb
end
"""
    fresnel_phase(n,theta,h,k)           
Computes the phase shift of propagation through the layers of a multilayer stack
# Arguments
* `n` : array of refractive indices 
* `theta` : array of angles                     
* `h` : array of  layer thicknesses                    
* `k` :  free-space wavenumber                    
# Output
* `delta`: array of phase shift angles
"""
function fresnel_phase(n,theta,h,k)
	delta=h.*k.*n[2:end-1].*cos.(theta[2:end-1])
	return delta
end
"""
    fresnel_transfermatrix(r,t,delta,rb=-r,tb=t)           
Computes the transfer matrix of the stack from fresnel coefficients
# Arguments
* `r`: array of fresnel reflection coefficients for forward propagation
* `t`: array of fresnel transmission coefficients for forward propagation
* `delta`: array of phase shift angles
* `rb`: optional array of fresnel reflection coefficients for backward propagation
* `tb`: optional array of fresnel transmission coefficients for backward propagation
# Output
* `M`: transfer matrix
"""
function fresnel_transfermatrix(r,t,delta,rb=-r,tb=t)
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
"""
    fresnel_total(r,t,delta,rb=-r,tb=t)           
Computes the total amplitude reflection and transmission coefficient of the layer stack from fresnel coefficients
# Arguments
* `r`: array of fresnel reflection coefficients for forward propagation
* `t`: array of fresnel transmission coefficients for forward propagation
* `delta`: array of phase shift angles
* `rb`: optional array of fresnel reflection coefficients for backward propagation
* `tb`: optional array of fresnel transmission coefficients for backward propagation
# Output
* `rtot`: total amplitude reflection coefficient
* `ttot`: total amplitude transmission coefficient
"""
function fresnel_total(r,t,delta,rb=-r,tb=t)
    #build transfer matrix
	M=fresnel_transfermatrix(r,t,delta,rb,tb)
    #reflection
	rtot=M[2,1]/M[1,1]
    #transmission
	ttot=1/M[1,1]
    #special case only one layer
	if length(r)==1
		rtot=r[1]
		ttot=t[1]
	end
	return rtot,ttot
end
"""
    fresnel_intensity(r,t,delta,n,θ,rb=r,tb=t)           
Computes the total power reflection and transmission coefficient of the layer stack from fresnel coefficients
# Arguments
* `r`: array of fresnel reflection coefficients for forward propagation
* `t`: array of fresnel transmission coefficients for forward propagation
* `delta`: array of phase shift angles
* `n`: array of refractive indices
* `θ`: array of ray angles
* `rb`: optional array of fresnel reflection coefficients for backward propagation
* `tb`: optional array of fresnel transmission coefficients for backward propagation
# Output
* `R`: total power reflection coefficient
* `T`: total power transmission coefficient
"""
function fresnel_intensity(r,t,delta,n,θ,rb=-r,tb=t)
    #get amplitudes
	r,t=fresnel_total(r,t,delta,rb,tb)
    #convert to power
	R=abs.(r).^2
	#normalize T by wave vector z component
	T=abs.(t).^2*real(n[end]*cos(θ[end]))/real(n[1]*cos(θ[1]))
	return R,T
end
