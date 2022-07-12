#https://edoc.ub.uni-muenchen.de/2367/1/Soennichsen_Carsten.pdf
#Plasmons in metal nanostructures, Carsten Sönnichsen, Doktorarbeit LMU
"""
    Csca(k,α)           
Computes the scattering cross-section of a particle or array from the polarizability                                                
# Arguments
* `k` : wavenumber of propagation outside the particle
* `α` : polarizability   
# Output
* `Csca`: scattering cross-section                      
"""

function Csca(k,α)
    return k.^4 .*abs.(α.^2)/(6π*ε0^2)
end
"""
    Csca(k,α)           
Computes the scattering cross-section of a particle or array from the polarizability                                                
# Arguments
* `k` : wavenumber of propagation outside the particle
* `α` : polarizability   
# Output
* `Cabs`: absorption cross-section                      
"""
function Cabs(k,α)
    return k.*imag.(α/ε0)
end
"""
    Cext(k,α)           
Computes the extinction cross-section of a particle or array from the polarizability                                                
# Arguments
* `k` : wavenumber of propagation outside the particle
* `α` : polarizability   
# Output
* `Cext`: extinction cross-section                      
"""
function Cext(k,α)
    return Cabs(k,α)+Csca(k,α)
end
