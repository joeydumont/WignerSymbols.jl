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