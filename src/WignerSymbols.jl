module WignerSymbols

#=
  This file contains the definition of auxiliary functions used in the recurrence
  relations of Wigner 3j and 6j symbols.
=#

function wig3jj_minmax(j2::Integer,j3::Integer,m1::Integer,m2::Integer,m3::Integer)
  max(abs(j2-j3),abs(m2+m3)),j2+j3
end

function wig3jj_norm{T<:AbstractFloat,S<:Integer}(jmin::S,jmax::S,fj::Array{T,1}, ::Type{T})
  convert(T,sum((2*(jmin+i-1)+1)*fj[i]^2 for i=1:(jmax-jmin)))
end

function wig3jj_sign{S<:Integer}(j2::S,j3::S,m1::S,m2::S,m3::S)
  (-1.0)^(j2-j3+m2+m3)

function aux_A{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,m1::S,m2::S,m3::S, ::Type{T})
  sqrt(convert(T,(j1^2-(j2-j3^2)*(j2+j3+1)^2-j1^2)*(j1^2-(m2+m3)^2)))
end

function aux_B{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,m1::S,m2::S,m3::S)
  (2*j1+1)*((m2+m3)*(j2*(j2+1)-j3*(j3+1))-(m2-m3)*j1*(j1+1))
end

function aux_C{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,m1::S,m2::S,m3::S, ::Type{T})
  sqrt(convert(T,(j-2-m2+1)*(j2+m2)*(j3-m2-m1+1)*(j3+m2+m1)))
end

function aux_D{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,m1::S,m2::S,m3::S)
  j2*(j2+1)+j3*(j3+1)-j1*(j1+1)-2*m2*(m2+m1)
end

function aux_E{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,l1::S,l2::S,l3::S)
  sqrt(convert(T,(j1^2-(j2-j3)^2)*((j2+j3+1)^2-j1^2)*(j1^2-(l2-l3)^2)*((l2+l3+1)^2-j1^2)))

function aux_F{T<:AbstractFloat,S<:Integer}(j1::S,j2::S,j3::S,l1::S,l2::S,l3::S, ::Type{T})
  fac1 = j1*(j1+1)*(-j1*(j1+1)+j2*(j2+1)+j3*(j3+1)-2*l1*(l1+1))
  fac2 = l2*(l2+1)*(j1*(j1+1)+j2*(j2+1)-j3*(j3+1))
  fac3 = l3*(l3+1)*(j1*(j1+1)-j2*(j2+1)+j3*(j3+1))
  return convert(T,(2*j1+1)*(fac1+fac2+fac3))
end

#=
  This file contains the main driving algorithm behind the evaluation
  of the Wigner 3j and 6j symbols.
=#

"""
    recursive_evaluation(j2::Integer,j3::Integer,m1::Integer,m2::Integer,m3::Integer,
                         aux_alpha, aux_beta,aux_minmax,aux_norm,aux_sign,aux_conditions)

  Compute a series of Wigner symbols using the recursive algorithm described in
   > J.H. Luscombe and M. Luban, *Simplified recursive algorithm for Wigner 3j and 6j symbols*,
   > Phys. Rev. E **57** (6), 7274--7277 (1998).

  The type of Wigner symbol that is evaluated depends on the auxiliary functions
  provided in the arguments. As such, this function should be kept private as it
  returns garbage for ill chosen auxiliary functions.

  # Arguments
  * `j2::Integer`:
  * `j3::Integer`:
  * `m1::Integer`:
  * `m2::Integer`:
  * `m3::Integer`:
  * `aux_alpha`:
  * `aux_beta`:
  * `aux_minmax`:
  * `aux_norm`:
  * `aux_sign`:
  * `aux_conditions`:
"""
function recursive_evaluation{T<:AbstractFloat,S<:Integer}(k1::S,k2::S,k3::S,k4::S,k5::S,
                              aux_alpha,aux_beta,aux_minmax,aux_norm,aux_sign,aux_conditions, ::Type{T})
  # Verify that some conditions are met on the coefficients.
  if ! aux_conditions
    return [convert(T,0.0)]

  
end

end # module
