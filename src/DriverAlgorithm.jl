#=
  This file contains the main driving algorithm behind the evaluation
  of the Wigner 3j and 6j symbols.
=#

"""
    recursive_evaluation(j2::Integer,j3::Integer,m1::Integer,m2::Integer,m3::Integer,
                         aux_alpha, aux_beta,aux_minmax,aux_norm,aux_sign)

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
"""
function recursive_evaluation(k1::Integer,k2::Integer,k3::Integer,k4::Integer,k5::Integer,
                              aux_alpha,aux_beta,aux_minmax,aux_norm,aux_sign)
  # Evaluation here.
end
