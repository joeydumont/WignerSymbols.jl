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
                              aux_alpha,aux_beta,aux_minmax,aux_norm,aux_sign,aux_conditions)
  # Verify that some conditions are met on the coefficients.
  if ! aux_conditions
    return [0.0]
end
