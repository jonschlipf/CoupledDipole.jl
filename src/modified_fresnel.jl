#from: https://pubs.acs.org/doi/full/10.1021/nn302879j equation 3/4
"""
    modified_r_s(α,pitch,n1,n2,θ1,k,ε0=ε0)
Computes modified fresnel coefficient for a surface with dipoles for s/TE polarization                                          
# Arguments
* `α` : polarizability  
* `pitch` : (mean) spacing of particles on array  
* `n1` : refractive index before interface  
* `n2` : refractive index after interface  
* `θ1` : incidence angle before interface  
* `k` : wavenumber of propagation outside the particle
* `ε0` : optional permittivity of free space
# Output
* `r`: forward reflection coefficient                      
* `t`: forward transmission coefficient                      
* `rb`: backward reflection coefficient                      
* `tb`: backward transmission coefficient                      
"""
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

"""
    modified_r_p(α,pitch,n1,n2,θ1,k,ε0=ε0)
Computes modified fresnel coefficient for a surface with dipoles for p/TM polarization                                          
# Arguments
* `α` : polarizability  
* `pitch` : (mean) spacing of particles on array  
* `n1` : refractive index before interface  
* `n2` : refractive index after interface  
* `θ1` : incidence angle before interface  
* `k` : wavenumber of propagation outside the particle
* `ε0` : optional permittivity of free space
# Output
* `r`: forward reflection coefficient                      
* `t`: forward transmission coefficient                      
* `rb`: backward reflection coefficient                      
* `tb`: backward transmission coefficient                      
"""
function modified_r_p(α,pitch,n1,n2,θ1,k,ε0=ε0)
    #see p polarization version
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

