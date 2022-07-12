using Pkg
Pkg.activate(".")
using ProgressMeter

#spectral range and material
include("../src/fresnel.jl")
E=1.4:.01:2.4
λ=1240 ./E
Ep=1.3e16*1240/(2π*3e8*1e9)
gamma=2.5e13*1240/(2π*3e8*1e9)
drudemod(E)=1 .-Ep^2 ./(E*(E+1im*gamma))
ε=[drudemod(e) for e in E]

#parameters of nanodisk
θ0=0
k=2π./λ
a=60
b=10
include("../src/CoupledDipole.jl")

#refractive indices and effective permittivity
n1=1.52
n2=1
εout=1.9
#actually computing L vs taking from Svedendahl
#esquared=1-b^2*a^-2
#g=√((1-esquared)/esquared)
#L=g/(2*esquared)*(π/2-atan(g))-.5*g^2
L=.09

#compute α as given by Svedendahl
keff=k.*sqrt(εout)
V=4/3*π*a^2*b
Elspr=sqrt.(Ep.^2*L./(εout*(1-L)+L))
c=1240/(2π)
#F=wlspr^2*a^2*2*b/(9c^3*L)      #blub
F=Elspr^2*a^2*2*b/(9c^3*L)      
α=V*Elspr.^2/L ./(Elspr^2 .-E.^2 .-1im*E.*(gamma.+F.*E.^2))*ε0

#compute modified fresnel of surface
rho=12.6e-6
pitch=rho^-.5
rint,tint,rintb,tintb=modified_r_s(α,pitch,n1,n2,θ0,keff)

#data caputred from literature
using CSV,DataFrames
sveddata=CSV.read("test/Svedendahl_normal.csv",DataFrame)

using Plots
plot(xlabel="E_phot / eV",ylabel="T")
#literature data
plot!(sveddata[!,:x],sveddata[!,:Curve1],label="Svedendahl")
#computed transmission
plot!(E,vec(abs.(tint.^2))*n2/n1,label="CDA, ε_eff=1.9")
using Optim
ESved=sveddata[!,:x]
λSved=1240 ./ESved
kSved=2π./λSved
εSved=[drudemod(l) for l in λSved]
#fit effective permittivity to svedendahl
#parametric transmission computation
function get_T(ts)
    keff=kSved.*real(sqrt(ts[1]+0im))
    V=4/3*π*a^2*b
    Elspr=sqrt.(Ep.^2*L./(ts[1]*(1-L)+L))
    c=1240/(2π)
    F=Elspr^2*a^2*2*b/(9c^3*L)      #blub
    α=V*Elspr.^2/L ./(Elspr^2 .-ESved.^2 .-1im*ESved.*(gamma.+F.*ESved.^2))*ε0
    rint,tint,rintb,tintb=modified_r_s(α,rho^-.5,n1,n2,θ0,keff)
    return vec(abs.(tint.^2))*n2/n1
end
#cost function
function err(ts)
    return sqrt(sum((get_T(ts)-sveddata[!,:Curve1]).^2))
end
#start value
ts0=[εout]
#fit
results=Optim.optimize(err,ts0,Optim.Options(time_limit=120,show_trace=true,x_tol=1e-3))
#show fit result
println(results.minimizer)
#plot data with fitted εout
plot!(ESved,get_T(results.minimizer),label=string("CDA, fitted ε_eff=",round(results.minimizer[1],digits=2)))
plot!()
savefig("test/svedendahl.png")
