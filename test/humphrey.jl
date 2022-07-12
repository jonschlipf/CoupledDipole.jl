#https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075404
#Plasmonic surface lattice resonances on arrays of different lattice symmetry
#Reproduces figure 1

#virtual environment for packages
using Pkg
Pkg.activate(".")
#parameters
diameter=120
height=30
pitch=480
n=1.515
λ=500:2:900
k=n*2π./λ
nshells=2000

#material from Palik
λAl=[494.9,516.6,539.1,563.6,590.4,619.9,652.6,688.8,729.3,774.9,826.6]
εAl=([.130,.130,.129,.120,.121,.131,.14,.14,.148,.143,.145].+1im*[2.88,3.07,3.25,3.45,3.66,3.88,4.15,4.44,4.74,5.09,5.5]).^2
using Interpolations
Al_palik=extrapolate(interpolate((λAl,),εAl,Gridded(Linear())),Interpolations.Flat())
ε=Al_palik(λ)

#load coupled dipole source code
include("../src/CoupledDipole.jl")
#single particle polarizability
α=αMLWA(diameter/2,height/2,ε,k,n^2,ε0)
#array factor
S=Ssquare(nshells,pitch,k)
#array polarizability
p=α_to_array(α,S)
#extinction cross-sections for both
C_single=Cext(k,α)
C_multi=Cext(k,p)
#show in plot
using Plots
plot(xlabel="λ / nm",ylabel="Cext / μm²",framestyle=:box)
plot!(λ,C_single*1e-6,label="single particle")
plot!(λ,C_multi*1e-6,label="particle in array")
plot!(leg=:topleft)

#plot in the paper
using CSV,DataFrames
barnesdata=CSV.read("test/humphrey.csv",DataFrame,delim=",",decimal=',')
plot!(barnesdata[!,:x],barnesdata[!,:Single],label="Single, Humphrey")
plot!(barnesdata[!,:x],barnesdata[!,:Array],label="Array, Humphrey")



savefig("test/humphrey.png")
