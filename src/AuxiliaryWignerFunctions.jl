#=
  This file contains the definition of auxiliary functions used in the recurrence
  relations of Wigner 3j and 6j symbols.
=#

function aux_A(j1,j2,j3,m1,m2,m3)
  sqrt((j1^2-(j2-j3^2)*(j2+j3+1)^2-j1^2)*(j1^2-(m2+m3)^2))
end

function aux_B(j1,j2,j3,m1,m2,m3)
  (2*j1+1)*((m2+m3)*(j2*(j2+1)-j3*(j3+1))-(m2-m3)*j1*(j1+1))
end

function aux_C(j1,j2,j3,m1,m2,m3)
  sqrt((j-2-m2+1)*(j2+m2)*(j3-m2-m1+1)*(j3+m2+m1))
end

function aux_D(j1,j2,j3,m1,m2,m3)
  j2*(j2+1)+j3*(j3+1)-j1*(j1+1)-2*m2*(m2+m1)
end

function aux_E(j1,j2,j3,l1,l2,l3)
  sqrt((j1^2-(j2-j3)^2)*((j2+j3+1)^2-j1^2)*(j1^2-(l2-l3)^2)*((l2+l3+1)^2-j1^2))

function aux_F
  fac1 = j1*(j1+1)*(-j1*(j1+1)+j2*(j2+1)+j3*(j3+1)-2*l1*(l1+1))
  fac2 = l2*(l2+1)*(j1*(j1+1)+j2*(j2+1)-j3*(j3+1))
  fac3 = l3*(l3+1)*(j1*(j1+1)-j2*(j2+1)+j3*(j3+1))
  return (2*j1+1)*(fac1+fac2+fac3)
end