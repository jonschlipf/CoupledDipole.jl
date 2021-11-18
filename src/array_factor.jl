#no oblique working yet
function Ssquare_s(nshells,pitch,k,θ=0,ε0=ε0) #only vertical incidence for now
S=0im*k # initialize complex array
kin=k*sin(θ)
for i=1:nshells
	#nearest neighbours horizontal and vertical
    r=pitch*i
    S=S+2exp.(1im*k*r+1im*r*kin).*(2*(1 .-1im*k*r)*r^-3) 
    S=S+2exp.(1im*k*r).*((-1)*(1 .-1im*k*r)*r^-3+1*k.^2/r)
	#diagonal
    r=pitch*i*√2
    S=S+4exp.(1im*k*r+1im*r*kin).*(.5*(1 .-1im*k*r)*r^-3+.5*k.^2/r) 
	#others
    for j=1:i-1
        r=pitch*sqrt(i^2+j^2)
        θ1=atan(j/i)
        θ2=atan(i/j)
        S=S+4exp.(1im*k*r+1im*kin*i*pitch).*((3*cos(θ1).^2-1)*r^-3*(1 .-1im*k*r)+k.^2 .*sin(θ1)^2/r) #S4
        S=S+4exp.(1im*k*r+1im*kin*j*pitch).*((3*cos(θ2).^2-1)*r^-3*(1 .-1im*k*r)+k.^2 .*sin(θ2)^2/r) #S5
    end
end
return S/(4π*ε0)
end
function Ssquare_p(nshells,pitch,k,θ=0,ε0=ε0) #only vertical incidence for now
S=0im*k # initialize complex array
kin=k*sin(θ)
for i=1:nshells
	#nearest neighbours horizontal and vertical
    r=pitch*i
    S=S+2exp.(1im*k*r).*(2*(1 .-1im*k*r)*r^-3) 
    S=S+2exp.(1im*k*r+1im*r*kin).*((-1)*(1 .-1im*k*r)*r^-3+1*k.^2/r)
	#diagonal
    r=pitch*i*√2
    S=S+4exp.(1im*k*r+1im*r*kin).*(.5*(1 .-1im*k*r)*r^-3+.5*k.^2/r) 
	#others
    for j=1:i-1
        r=pitch*sqrt(i^2+j^2)
        θ1=atan(j/i)
        θ2=atan(i/j)
        S=S+4exp.(1im*k*r+1im*kin*j*pitch).*((3*cos(θ1).^2-1)*r^-3*(1 .-1im*k*r)+k.^2 .*sin(θ1)^2/r) #S4
        S=S+4exp.(1im*k*r+1im*kin*i*pitch).*((3*cos(θ2).^2-1)*r^-3*(1 .-1im*k*r)+k.^2 .*sin(θ2)^2/r) #S5
    end
end
return S/(4π*ε0)
end
