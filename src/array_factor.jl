using ProgressMeter
#see p version for comments
#this is slow, there should be more efficient algorithms that exploit lattice symmetry
#https://pure.tue.nl/ws/files/3897048/728024624764265.pdf
#Citation for published version (APA):Zijlstra, P., Orrit, M., & Koenderink, A. F. (2014). Metal nanoparticles for microscopy and spectroscopy. In C.Mello Donegá, de (Ed.),Nanoparticles : Workhorses of Nanoscience Springer. https://doi.org/10.1007/978-3-662-44823-6_3
#Equation 3.31
#oblique implementation was not tested fully

"""
    Ssquare(nshells,pitch,k,θ=0,α=0,ε0=ε0)           
Computes the array factor in the coupled dipole approximation for a square lattice
# Arguments
* `nshells` : number of shells of nearest neighbors in the lattice to consider 
* `pitch` : lattice pitch 
* `k` : wavenumber of propagation outside the particle
* `θ` : optional incidence angle
* `α` : optional azimuth angle 
* `ε0` : optional permittivity of free space
# Output
* `S`: array factor
"""


function Ssquare(nshells,pitch,k,θ=0,α=0,ε0=ε0)
	#preallocate
	S=0im*k
    kx=k*sin(θ)*cos(α) #x component of k
    ky=k*sin(θ)*sin(α) #y component of k
	function Spart(i,j) #call for each element
		r=pitch*sqrt(i^2+j^2) #dipole distance
		φ=atan(j,i) #angle
        phase=k*r+kx*pitch*i+ky*pitch*j
		S1=exp.(1im*phase).*(3*cos(φ)^2-1).*(1 .-1im*k*r)*r^-3 #first part
        S2=exp.(1im*phase).*(sin(φ)^2).*k.^2*r^-1 #second part
		return S1+S2
	end
	    @showprogress for i=-nshells:nshells #iterate over rows
        for j=-nshells:nshells #iterate over columns
            if (i != 0)||(j!=0) #exclude the element itself
                S=S+Spart(i,j) 
            end
        end
    end
    return S/(4π*ε0) #unit system conversion
end
"""
    α_to_array(α,S)
Computes array polarizability from particle polarizability and array factor
# Arguments
* `α` : particle polarizability  
* `S` : array factor
# Output
* `p`: array polarizability                      
"""
function α_to_array(α,S)
    return 1 ./((1 ./α)-S)
end

