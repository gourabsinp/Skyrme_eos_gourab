import numpy as np

import matplotlib.pyplot as plt

import h5py

from tov_new_general import TOV_SOLVER

q=np.linspace(0.01,1.92,500)

nb=[]

M=[]

R=[]

P_c=[]

e_c=[]

rho_c=[]

for i in range(len(q)):
   M_max,R_max,P_max,E_max,rho_max=TOV_SOLVER(q[i])
   nb.append(q[i])
   M.append(M_max)
   R.append(R_max)
   P_c.append(P_max)
   e_c.append(E_max)
   rho_c.append(rho_max)
   
   
with h5py.File('SLY4_EOS_M_R_general.h5','w') as hdf:
      G1=hdf.create_group('parameter')
      G1.create_dataset('nb',data=nb)
      G1.create_dataset('M',data=M)
      G1.create_dataset('R',data=R)
      G1.create_dataset('P',data=P_c)
      G1.create_dataset('e',data=e_c)
      G1.create_dataset('rho',data=rho_c)


plt.scatter(R,M)
plt.xlabel("R[Km]")
plt.ylabel(r"M$[M_0]$")
plt.savefig("M_R_curve_TOV3_general.png")
plt.show()
plt.close()
plt.scatter(rho_c,M)
plt.xlabel(r"$\rho_c[gm/cm^3]$")
plt.ylabel(r"M$[M_0]$")
plt.savefig("M_RHO_curve_TOV3_general.png")
plt.show()
plt.close()
plt.scatter(rho_c,R)
plt.xlabel(r"$\rho_c[gm/cm^3]$")
plt.ylabel("R[Km]")
plt.savefig("R_RHO_curve_TOV3_general.png")
plt.show() 
plt.close()
plt.scatter(e_c,M)
plt.xlabel(r"$\epsilon_c[erg/cm^3]$")
plt.ylabel(r"M$[M_0]$")
plt.savefig("M_epsilon_curve_TOV3_general.png")
plt.show()
plt.close()
plt.scatter(e_c,R)
plt.xlabel(r"$\epsilon_c[erg/cm^3]$")
plt.ylabel("R[Km]")
plt.savefig("R_epsilon_curve_TOV3_general.png")
plt.show()   

