using Pkg
Pkg.activate(".")
using ProgressMeter

include("../src/CoupledDipole.jl")
λ=600:10:2000
#get Al permittivity from RCWA library
using RigorousCoupledWaveAnalysis
Al=ModelPerm(RigorousCoupledWaveAnalysis.al_rakic)
ε=[get_permittivity(Al,l) for l in λ]


θ0=0
n=1
k0=2π./λ
k=n*2π./λ
h=[0,0]
a=50
b=50
pitch=1000
nshells=1000
α=αMLWA(a,b,ε,k,n^2,ε0)
S=Ssquare(nshells,pitch,k,θ0)
p=α_to_array(α,S)
rint,tint,rintb,tintb=modified_r_s(p,pitch,n,n,θ0,k)
println("array computed")
R=zeros(length(λ))
T=zeros(length(λ))
@showprogress for i=1:length(λ)
    narr=[n,n,n,n]
    θ=fresnel_angles(narr,θ0)
    r,t,rb,tb=fresnel_coefs_s(narr,θ)
    r[2]=rint[i]
    t[2]=tint[i]
    rb[2]=rintb[i]
    tb[2]=tintb[i]
    delta=fresnel_phase(narr,θ,h,k0[i])
    R[i],T[i]=fresnel_intensity(r,t,delta,n,θ,rb,tb)
end
using CSV,DataFrames
lumdata=CSV.read("test/npa_scpml_smol_sphere.txt",DataFrame)
using Plots
plot(λ,R,label="R")

plot!(lumdata[!,:lambda]*1e9,lumdata[!,:R],label="FDTD")
rho=pitch^-2
#p=p*correction
plot!(λ,abs.(-1im*k.*p.*rho/ε0./(2 .-1im*k.*p*rho/ε0)).^2)
#p=p*correction
plot!(λ,abs.(-1im*k.*p.*rho/ε0./(2 .+1im*k.*p*rho/ε0)).^2)

