import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib as mpl
import matplotlib.lines as mlines
plt.rcParams['xtick.labelsize']=10
plt.rcParams['ytick.labelsize']=10
plt.rcParams['xtick.major.pad']='5.0'
plt.rcParams['ytick.major.pad']='5.0'
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams["axes.labelsize"] = 10
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.labelcolor"] = "black"
mpl.rcParams["xtick.color"] = "black"
mpl.rcParams["ytick.color"] = "black"
plt.rcParams["xtick.labelcolor"] = "black"
plt.rcParams["ytick.labelcolor"] = "black"


import skyrm2

#from skyrm2 import rho_sat

from skyrm2 import beta_eos

from skyrm2 import E_A_pnm

from skyrm2 import E_A_snm

from skyrm2 import E_sym_A_def

from skyrm2 import E_A_sym_2

from skyrm2 import E_A_sym_4

from skyrm2 import E_A_sym_6

from skyrm2 import P_snm

#from skyrm2 import P_bary

from skyrm2 import P_A_pnm

from skyrm2 import P_A_sym_2

from skyrm2 import P_A_sym_4

from skyrm2 import P_A_sym_6

#rho_sat_new=rho_sat(0.1)

rho_sat_new=0.16

print(rho_sat_new)

rho=np.linspace(0.0001,1.92,1000)

rho1=[]
for i in range(len(rho)):
  rho1.append(rho[i])

I_new=[]

rho_p_new=[]

rho_n_new=[]

rho_e_new=[]

rho_mu_new=[]

P_e_new=[]

P_mu_new=[]

epsilon_e_new=[]

epsilon_mu_new=[]

mu_e_new=[]

mu_mu_new=[]

P_beta_new=[]

a_sym2=[]

E_A_pnm_new=[]

E_A_snm_new=[]

E_A_sym_def_new=[]

Press_snm=[]

Press_pnm=[]

Press_sym=[]

P_baryon=[]

rho_sat_new1=[]

for i in range(len(rho)):
    I,rhop,rhon,rhoe,rhomu,Pe,Pmu,epsilone,epsilonmu,mue,mumu,Pb,Pbeta=beta_eos(rho[i],1.0)
    I_new.append(I)
    rho_p_new.append(rhop)
    rho_n_new.append(rhon)
    rho_e_new.append(rhoe)
    rho_mu_new.append(rhomu)
    P_e_new.append(Pe)
    P_mu_new.append(Pmu)
    epsilon_e_new.append(epsilone)
    epsilon_mu_new.append(epsilonmu)
    mu_e_new.append(mue)
    mu_mu_new.append(mumu)
    P_beta_new.append(Pbeta)
    P_baryon.append(Pb)
    
P_e_b=[]
P_mu_b=[]
P_bary_b=[]
for i in range(len(rho)):
    P1=P_e_new[i]/P_beta_new[i]
    P2=P_mu_new[i]/P_beta_new[i]
    P3=P_baryon[i]/P_beta_new[i]
    P_e_b.append(P1)
    P_mu_b.append(P2)
    P_bary_b.append(P3)
    
    
    
frac_e_b=[]
frac_mu_b=[]
frac_proton_b=[]
frac_neutron_b=[]
for i in range(len(rho)):
    P1=rho_e_new[i]/rho[i]
    P2=rho_mu_new[i]/rho[i]
    P3=rho_p_new[i]/rho[i]
    P4=rho_n_new[i]/rho[i]
    frac_e_b.append(P1)
    frac_mu_b.append(P2)
    frac_proton_b.append(P3)
    frac_neutron_b.append(P4)
    
    
for i in range(len(rho)):
    #a_sym2.append(E_A_sym_2(rho[i])+E_A_sym_4(rho[i])+E_A_sym_6(rho[i]))
    #E_A_sym_4_corr.append(E_A_sym_4(rho[i],I_new[i]))
    #E_A_sym_6_corr.append(E_A_sym_6(rho[i],I_new[i]))
    E_A_pnm_new.append(E_A_pnm(rho[i],I_new[i]))
    E_A_snm_new.append(E_A_snm(rho[i],I_new[i]))
    E_A_sym_def_new.append(E_sym_A_def(rho[i],I_new[i]))
    #E_A_sym_approx_new.append(E_sym_A_approx(rho[i],I_new[i]))
    #rho_sat_new1.append(rho_sat_new)
    Press_snm.append(P_snm(rho[i]))
    #Press_sym.append(P_A_sym_2(rho[i])+P_A_sym_4(rho[i])+P_A_sym_6(rho[i]))
    Press_pnm.append(P_A_pnm(rho[i]))
    
#relative_ratio=[]
#for i in range(len(rho)):
    #A=1-(a_sym2[i]/E_A_sym_def_new[i])
    #relative_ratio.append(A)
    
#with h5py.File('SLY4_EOS_SKYRM_SYM_ENER_upto__general_corr_new.h5','w') as hdf:
      #G1=hdf.create_group('parameter')
      #G1.create_dataset('I',data=I_new)
      #G1.create_dataset('rhop',data=rho_p_new)
      #G1.create_dataset('rhon',data=rho_n_new)
      #G1.create_dataset('rhoe',data=rho_e_new)
      #G1.create_dataset('rhomu',data=rho_mu_new)
      #G1.create_dataset('Pe',data=P_e_new)
      #G1.create_dataset('Pmu',data=P_mu_new)
      #G1.create_dataset('epsilone',data=epsilon_e_new)
      #G1.create_dataset('epsilonmu',data=epsilon_mu_new)
      #G1.create_dataset('mue',data=mu_e_new)
      #G1.create_dataset('mumu',data=mu_mu_new)
      #G1.create_dataset('Pressure_beta_equi',data=P_beta_new)
      #G1.create_dataset('press_baryon',data=P_baryon)
      #G1.create_dataset('total_baryon_den',data=rho1)
      #G1.create_dataset('sym_energy_per_A',data=a_sym2)
      #G1.create_dataset('pnm_minus_snm_per_A',data=E_A_sym_def_new)
      #G1.create_dataset('pnm_energy_per_A',data=E_A_pnm_new)
      #G1.create_dataset('snm_energy_per_A',data=E_A_snm_new)
      #G1.create_dataset('snm_press',data=Press_snm)
      #G1.create_dataset('sym_press',data=Press_sym)
      #G1.create_dataset('pnm_press',data=Press_pnm)
      #G1.create_dataset('frac_e_press',data=P_e_b)
      #G1.create_dataset('frac_mu_press',data=P_mu_b)
      #G1.create_dataset('frac_bary_press',data=P_bary_b)
      #G1.create_dataset('frac_e_number',data=frac_e_b)
      #G1.create_dataset('frac_mu_number',data=frac_mu_b)
      #G1.create_dataset('frac_n_number',data=frac_neutron_b)
      #G1.create_dataset('frac_p_number',data=frac_proton_b)
      #G1.create_dataset('sym_comparison',data=relative_ratio)


I_old=[]
for i in range(len(rho)):
  r1=rho_n_new[i]-rho_p_new[i]
  r2=rho_n_new[i]+rho_p_new[i]
  I_old.append(r1/r2)


print(rho_n_new[0],rho_n_new[-1])

print(rho_p_new[0],rho_p_new[-1])


#print(len(I_new))
#print(len(rho))
print((P_baryon))

#print(P_mu_new)

#print(rho)


plt.plot(rho/rho_sat_new,Press_snm,color='r',label='c=snm')
plt.plot(rho/rho_sat_new,Press_pnm,color='g',label='c=pnm')
#plt.plot(rho/rho_sat_new,Press_sym,color='b',label='c=sym')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$P_{c}$")
plt.legend()
#plt.savefig("pnm_snm_press.png")
plt.show()




#plt.plot(rho/rho_sat_new,E_A_snm_new,color='r',label='c=snm')
#plt.plot(rho/rho_sat_new,E_A_pnm_new,color='g',label='c=pnm')
#plt.plot(rho/rho_sat_new,a_sym2,color='b',label='c=sym')
#plt.xlabel(r"$\rho/\rho_0$")
#plt.ylabel(r"$\frac{E_{c}}{A}$")
#plt.legend()
#plt.savefig("I_new.png")
#plt.show()



#plt.plot(rho/rho_sat_new,a_sym2,color='r',label='$a_s^{(2)}(\rh)$')
#plt.plot(rho/rho_sat_new,E_A_sym_def_new,color='g',label='Pnm-Snm')
#plt.xlabel(r"$\rho/\rho_0$")
#plt.ylabel(r"$E_{sym}$")
#plt.legend()
#plt.savefig("I_new.png")
#plt.show()





#plt.plot(rho/rho_sat_new,relative_ratio,color='r')
#plt.xlabel(r"$\rho/\rho_0$")
#plt.ylabel(r"$(\frac{a_s^{(2)}}{Pnm-Snm})-1$")
#plt.legend()
#plt.savefig("I_new.png")
#plt.show()




#plt.plot(I_new,mun_mup_new,color='r',label='s=mun-mup')
#plt.plot(rho/rho_sat_new,mu_e_new,color='g',label='s=e')
#plt.xlabel(r"$\rho/\rho_0$")
#plt.ylabel(r"$\mu_n-\mu_p$")
#plt.legend()
#plt.savefig("I_new.png")
#plt.show()


plt.plot(rho/rho_sat_new,P_e_b,color='r',label='s=e')
plt.plot(rho/rho_sat_new,P_mu_b,color='g',label='s=mu')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$P_s/P_{\beta}$")
plt.legend()
#plt.savefig("I_new.png")
plt.show()



plt.plot(rho/rho_sat_new,I_new,color='r',label='s=I_new')
#plt.plot(rho/rho_sat_new,I_old,color='g',label='s=I_old')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$\mathtt{I_s}$")
plt.legend()
#plt.savefig("I_new.png")
plt.show()

plt.plot(rho,rho_p_new,color='r',label='s=p')
plt.plot(rho,rho_n_new,color='g',label='s=n')
plt.xlabel(r"$\rho$")
plt.ylabel(r"$\rho_s$")
plt.legend()
#plt.savefig("I_new.png")
plt.show()

#plt.plot(rho)



plt.plot(rho/rho_sat_new,P_beta_new,color='r',label='s=beta-pressure')
plt.plot(rho/rho_sat_new,P_baryon,color='g',label='s=baryon-pressure')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$P_{\beta}$")
plt.legend()
#plt.savefig("I_new.png")
plt.show()

#plt.plot(rho/rho_sat_new,E_A_snm_new,color='r')
#plt.xlabel(r"$\rho/\rho_0$")
#plt.ylabel(r"$E_{snm}$")
#plt.legend()
#plt.savefig("I_new.png")
#plt.show()

plt.plot(rho/rho_sat_new,Press_snm,color='r')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$P_{snm}$")
#plt.legend()
#plt.savefig("I_new.png")
plt.show()


plt.plot(rho/rho_sat_new,P_baryon,color='r')
plt.xlabel(r"$\rho/\rho_0$")
plt.ylabel(r"$P_{baryon}$")
#plt.legend()
#plt.savefig("I_new.png")
plt.show()





