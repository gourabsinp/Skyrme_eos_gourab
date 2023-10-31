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
 
h_c=197.32698                                                         # Mev-fm unit(it is actually planck constant times speed of light ).
m_n_c_2=939.0                                                         # Mev unit(it is equivalent mass energy of nucleon).
pi=3.141592653589793
pi2=pi**2
m_e_c_2=0.511                                                         # Mev unit(it is equivalent mass energy of electron)      
m_mu_c_2=105.6583755                                                  # Mev unit(it is equivalent mass energy of muon)


#SLY4 EOS
t0=-2488.91 
t1=486.82 
t2=-546.39
t3=13777.0
x0=0.834
x1=-0.344
x2=-1.000
x3=1.354
sigma=1.0/6.0


# Define First order derivative of isospin asymmetry factor(I) with respect to baryon density(rho):

def DI(rho,I):
    d=-(I/rho)
    return(d)
    
    
# Define asymmetry factor(F):

def F(m,I):
    a=((1.0+I)**m)+((1.0-I)**m)
    X=a/2.0
    return(X)
    
    
# Define First order derivative of asymmetry factor(F) with respect to Isospin asymmetry(I):


def D1F(m,I):
    F1=m*(((1.0+I)**(m-1.0))-((1.0-I)**(m-1.0)))
    X=F1/2.0
    return(X)
    
    
# Now the energy density, e=T_k+H0+H3+Heff, 
#  where T_k is the kinetic energy density,
#  H0 is a zero range term, 
#  H3 is the density dependent term,
#  Heff is the effective-mass term.
# Now I will define each term of energy density.

 
# First I define kinetic energy density term(T_k):

def T_k(rho,I):
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    m=5/3
    c=(rho**(5/3))*F(m,I)
    d=(3*a*(b**(2/3))*c)/5
    return(d)
    

# Now I define the zero range term(H0):

def H0(rho,I):
    m=2
    a=(2*(2+x0))-((2*x0+1)*F(m,I))
    b=(t0*(rho**2)*a)/8
    return(b)
   
   
# Now I define the density dependent term(H3):


def H3(rho,I):
    m=2
    a=(2*(2+x3))-((2*x3+1)*F(m,I))
    b=(t3*(rho**(sigma+2))*a)/48
    return(b)
    
    
# Now I define the effective mass term(Heff):


def Heff(rho,I):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**(8/3))*func3)/40
    return(d)
    
    

# Now above I define all components of energy density T_k, H0, H3, Heff. So now I define total energy density e:


def e(rho,I):
    add=T_k(rho,I)+H0(rho,I)+H3(rho,I)+Heff(rho,I)+rho*(m_n_c_2)
    return(add)
    
    
# Now I define energy per nucleon term(E/A):


def E_A(rho,I):
    E_A=e(rho,I)/rho
    return(E_A)
    
    
# Now I define the symmetry nuclear matter energy per nucleon(E_A_snm):


def E_A_snm(rho,I):
        E_snm=E_A(rho,0.0)
        return(E_snm)
        
        
# Now I define the pure neutron matter energy per nucleon(E_A_pnm):


def E_A_pnm(rho,I):
        E_pnm=E_A(rho,1.0)
        return(E_pnm)       
        
        
        
#First order derivative of E_A_snm with respect to rho(D1E_A_snm):


def D1E_A_snm(rho):
    A=(h_c**2)/(5.0*m_n_c_2)
    B=(3*pi2)/2.0
    func=(3*t1)+((5+(4*x2))*t2)
    first_term=(A*(B**(2/3)))*(rho**(-1/3))
    second_term=(3*t0)/8.0
    third_term=(t3*(sigma+1)*(rho**sigma))/16.0
    fourth_term=(func*(B**(2/3))*(rho**(2/3)))/16.0
    X=first_term+second_term+third_term+fourth_term
    return(X)
    
    
# Symmetric nuclear matter pressure:
 
 
def P_snm(rho):
    P=(rho**2)*D1E_A_snm(rho)
    return(P)
    
    
# Now I define the symmetry energy per nucleon using its basic defination ((E_A)_sym=(E_A)_pnm-(E_A)_snm).


def E_sym_A_def(rho,I):
    E_sym=E_A_pnm(rho,I)-E_A_snm(rho,I)
    return(E_sym)        
        


# Now I define the symmetry energy per nucleon for second order correction a_s^(2)=1/2!(d^2(E/A)/dI^2)|I=0. 

    
def E_A_sym_2(rho):
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    func1=((2*x0)+1)
    func2=((2*x3)+1)
    func3=(3*t1*x1)-(4*t2)-(5*t2*x2)
    term1=A*(B**(2/3))*(rho**(2.0/3.0))
    term2=(t0*func1*rho)/8.0
    term3=(t3*func2*(rho**(sigma+1)))/48.0
    term4=((B**(2.0/3.0))*func3*(rho**(5.0/3.0)))/24.0
    return(term1-term2-term3-term4)



def dE_A_sym_2(rho):    
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    func1=((2*x0)+1)
    func2=((2*x3)+1)
    func3=(3*t1*x1)-(4*t2)-(5*t2*x2)
    term1=(2*A*(B**(2.0/3.0))*(rho**(-1/3)))/3.0
    term2=(t0*func1)/8
    term3=(t3*func2*(sigma+1.0)*(rho**(sigma)))/48.0
    term4=(5*(B**(2.0/3.0))*func3*(rho**(2.0/3.0)))/72.0
    return(term1-term2-term3-term4)


def P_A_sym_2(rho):
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    func1=((2*x0)+1)
    func2=((2*x3)+1)
    func3=(3*t1*x1)-(4*t2)-(5*t2*x2)
    term1=(2*A*(B**(2.0/3.0))*(rho**(5.0/3.0)))/3.0
    term2=(t0*(rho**2.0)*func1)/8.0
    term3=(t3*func2*(sigma+1.0)*(rho**(sigma+2)))/48.0
    term4=(5*(B**(2.0/3.0))*func3*(rho**(8.0/3.0)))/72.0
    return(term1-term2-term3-term4)
    
    
    
    
def E_A_sym_4(rho):
    m1=5.0/3.0
    D4F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)
    m2=8.0/3.0
    D4F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(2.0/3.0))*D4F53
    term1=(3*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D4F53)+((func2*D4F83)/2)
    term2=(3*(b**(2/3))*(rho**(5/3))*func3)/40.0
    t=(term1+term2)/24.0
    return(t)
    
    
def dE_A_sym_4(rho):
    m1=5.0/3.0
    D4F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)
    m2=8.0/3.0
    D4F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(-1.0/3.0))*D4F53
    term1=(2*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D4F53)+((func2*D4F83)/2)
    term2=((b**(2/3))*(rho**(2.0/3.0))*func3)/8.0
    t=(term1+term2)/24.0
    return(t)
    
    
def P_A_sym_4(rho):
    m1=5.0/3.0
    D4F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)
    m2=8.0/3.0
    D4F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(5.0/3.0))*D4F53
    term1=(2*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D4F53)+((func2*D4F83)/2)
    term2=((b**(2/3))*(rho**(8.0/3.0))*func3)/8.0
    t=(term1+term2)/24.0
    return(t)
    
    
    
def E_A_sym_6(rho):
    m1=5.0/3.0
    D6F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)*(m1-4.0)*(m1-5.0)
    m2=8.0/3.0
    D6F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)*(m2-4.0)*(m2-5.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(2.0/3.0))*D6F53
    term1=(3*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D6F53)+((func2*D6F83)/2)
    term2=(3*(b**(2/3))*(rho**(5/3))*func3)/40.0
    t=(term1+term2)/720.0
    return(t)
    
    

def dE_A_sym_6(rho):
    m1=5.0/3.0
    D6F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)*(m1-4.0)*(m1-5.0)
    m2=8.0/3.0
    D6F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)*(m2-4.0)*(m2-5.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(-1.0/3.0))*D6F53
    term1=(2*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D6F53)+((func2*D6F83)/2)
    term2=((b**(2/3))*(rho**(2.0/3.0))*func3)/8.0
    t=(term1+term2)/720.0
    return(t)
    
    
def P_A_sym_6(rho):
    m1=5.0/3.0
    D6F53=m1*(m1-1.0)*(m1-2.0)*(m1-3.0)*(m1-4.0)*(m1-5.0)
    m2=8.0/3.0
    D6F83=m2*(m2-1.0)*(m2-2.0)*(m2-3.0)*(m2-4.0)*(m2-5.0)
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    c=(rho**(5.0/3.0))*D6F53
    term1=(2*a*(b**(2.0/3.0))*c)/5
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*D6F53)+((func2*D6F83)/2)
    term2=((b**(2/3))*(rho**(8.0/3.0))*func3)/8.0
    t=(term1+term2)/720.0
    return(t)
    
    
def D1E_A_pnm(rho):
    A=(h_c**2)/(5.0*m_n_c_2)
    B=(3*pi2)/2.0
    func=(3*t2*(1+x2))+(t1*(1-x1))
    first_term=(A*(B**(2/3)))*(rho**(-1/3))*(2**(2.0/3.0))
    second_term=((2-(2*x0))*t0)/8.0
    third_term=(t3*(2-(2*x3))*(sigma+1)*(rho**sigma))/48.0
    fourth_term=(func*(B**(2/3))*(rho**(2/3))*(2**(2.0/3.0)))/8.0
    X=first_term+second_term+third_term+fourth_term
    return(X)
    
    
    
def P_A_pnm(rho):
    A=(rho**2)*D1E_A_pnm(rho)
    return(A)


# Calculation of the saturation density using newton raphson method:


def rho_sat(rho0):
    rho_select=np.zeros([1000001])
    rho_select[0]=rho0
    rho_sat=[]
    for i in range(1,len(rho_select)):
          rho_select[i]=rho_select[i-1]-(D1E_A_snm(rho_select[i-1])/D2E_A_snm(rho_select[i-1]))
          if(abs(rho_select[i]-rho_select[i-1])<0.0000001):
               return(rho_select[i-1])
            
            
# Chemical potential of electron:

def mu_e(rho,I,order_approx):
  if order_approx==1.0:
       X=(4*E_A_sym_2(rho)*(I))
       return(X)
  if order_approx==2.0:
       X=(4*E_A_sym_2(rho)*I)+(8*E_A_sym_4(rho)*(I**3))
       return(X)
  if order_approx==3.0:
       X=(4*E_A_sym_2(rho)*I)+(8*E_A_sym_4(rho)*(I**3))+(12*E_A_sym_6(rho)*(I**5))
       return(X)
  if order_approx==0.0:
       a=(6*(h_c**2))/(10*m_n_c_2)
       b=(3*pi2)/2
       DF53=(5*(((1+I)**(2.0/3.0))-((1-I)**(2.0/3.0))))/6.0
       DF2=2*I
       DF83=(8*(((1+I)**(5.0/3.0))-((1-I)**(5.0/3.0))))/6.0
       c=(rho**(-1/3))*DF53
       term1=a*(b**(2.0/3.0))*c
       term2=(t0*((2*x0)+1)*DF2)/4.0
       term3=(t3*((2*x3)+1)*DF2*(rho**sigma))/24.0
       f=(t1*(x1+2))+(t2*(x2+2))
       g=((t2*((2*x2)+1))-(t1*((2*x1)+1)))/2.0
       term4=(3*(b**(2.0/3.0))*(rho**(2/3))*((f*DF53)+(g*DF83)))/20
       term=rho*(term1-term2-term3+term4)
       return(term)    
    
    
# Electron number density calculation:
# mu_e^2=h_c^2 k_f_e^2 + m_e_c_2^2 , 
# also k_f_e^3=3 pi2 rho_e, 
# using these two relation we get rho_e=[(mu_e^2-m_e_c_2)/(h_c^2)]^(3/2)/(3 pi2)
    
    
def rho_e(rho,I,order_approx):
    X1=mu_e(rho,I,order_approx)**2-m_e_c_2**2
    Y1=X1/(h_c**2)
    Z1=(Y1**(3.0/2.0))/(3*pi2)
    return(Z1)
    
    
def k_f_e(rho,I,order_approx):
    X1=mu_e(rho,I,order_approx)**2-m_e_c_2**2
    Y1=X1/(h_c**2)
    Z1=(Y1**(1.0/2.0))
    return(Z1)
    
    
def epsilon_e(rho,I,order_approx):
    x=(h_c*k_f_e(rho,I,order_approx))/m_e_c_2
    A=(m_e_c_2**4)/(8*pi2*(h_c**3))
    e1_term1=(1+(x**2))**(1.0/2.0)
    e1_term2=1+(2*(x**2))
    e1=x*e1_term1*e1_term2
    #e2=np.log(x+e1_term1)
    e2=np.arcsinh(x)
    epsilon=A*(e1-e2)
    return(epsilon)
    
def P_e(rho,I,order_approx):
    P=(rho_e(rho,I,order_approx)*mu_e(rho,I,order_approx))-epsilon_e(rho,I,order_approx)
    return(P)




# Muon number density calculation:
# mu_mu^2=h_c^2 k_f_mu^2 + m_mu_c_2^2 , 
# also k_f_mu^3=3 pi2 rho_mu, 
# using these two relation we get rho_mu=[(mu_mu^2-m_mu_c_2)/(h_c^2)]^(3/2)/(3 pi2)
# Now according to the beta equilibrium condition mu_n=mu_p+mu_e and mu_e=mu_mu, so we can 
# also write the muon number density as, rho_mu=[(mu_e^2-m_mu_c_2)/(h_c^2)]^(3/2)/(3 pi2)


def rho_mu(rho,I,order_approx):
    X2=mu_e(rho,I,order_approx)**2-m_mu_c_2**2
    Y2=X2/(h_c**2)
    Z2=(Y2**(3.0/2.0))/(3*pi2)
    return(Z2)
    

def k_f_mu(rho,I,order_approx):
    X2=mu_e(rho,I,order_approx)**2-m_mu_c_2**2
    Y2=X2/(h_c**2)
    Z2=(Y2**(1.0/2.0))
    return(Z2)
    

def epsilon_mu(rho,I,order_approx):
    x=(h_c*k_f_mu(rho,I,order_approx))/m_mu_c_2
    A=(m_mu_c_2**4)/(8*pi2*(h_c**3))
    e1_term1=(1+(x**2))**(1.0/2.0)
    e1_term2=1+(2*(x**2))
    e1=x*e1_term1*e1_term2
    #e2=np.log(x+e1_term1)
    e2=np.arcsinh(x)
    epsilon=A*(e1-e2)
    return(epsilon)
    
def P_mu(rho,I,order_approx):
    P=(rho_mu(rho,I,order_approx)*mu_e(rho,I,order_approx))-epsilon_mu(rho,I,order_approx)
    return(P)
    
    
    



#Now I calculate the pressure using the formula P=rho^2 d(E/A)/d(rho) or P=rho^2 d(e/rho)/d(rho)
#Now we know that (E/A)=e/rho and in above I define all components of 
#energy density. So now I define the pressure component for each component
#of energy density and after it adding all component I calculate the total pressure.


# Define the pressure due to kinetic energy density term(T_k):


def P_T_k(rho,I):
    A=(3*(h_c**2))/(10.0*m_n_c_2)
    B=(3*pi2)/2.0
    m1=5.0/3.0
    term1=(2*(rho**(-1/3))*F(m1,I))/3
    term=(rho**2)*A*(B**(2/3))*term1
    return(term)
    
    
    
# Define the pressure due to zero range term(H0):


def P_H0(rho,I):
    m1=2
    term1=((2*(x0+2))-(((2*x0)+1)*F(m1,I)))
    term=((rho**2)*t0*term1)/8
    return(term)




# Define the pressure due to density dependent term(H3):


def P_H3(rho,I):
    m1=2
    term1=(sigma+1)*(rho**sigma)*((2*(x3+2))-(((2*x3)+1)*F(m1,I)))
    term=((rho**2)*t3*term1)/48
    return(term)
    
    
    

# Define the pressure due to effective mass term(Heff):



def P_Heff(rho,I):
    m1=5/3
    m2=8/3
    B=(3*pi2)/2.0
    func1=((t1*(x1+2))+(t2*(x2+2)))*F(m1,I)
    func2=((t2*((2*x2)+1))-(t1*((2*x1)+1)))*F(m2,I)
    term1=(5*(rho**(2/3))*(func1+(func2/2)))/3
    term=((rho**2)*3*(B**(2/3))*term1)/40
    return(term)
    
    


# Now define the total baryonic pressure(P_bary) without any approximation:


def P_general(rho,I):
    P=P_T_k(rho,I)+P_H0(rho,I)+P_H3(rho,I)+P_Heff(rho,I)
    return(P)
    



def P_bary(rho,I,order_approx):
    if order_approx==1.0: 
        P=P_snm(rho)+((I**2)*P_A_sym_2(rho))
        return(P)
    if order_approx==2.0:
        P=P_snm(rho)+((I**2)*P_A_sym_2(rho))+((I**4)*P_A_sym_4(rho))
        return(P)
    if order_approx==3.0:
        P=P_snm(rho)+((I**2)*P_A_sym_2(rho))+((I**4)*P_A_sym_4(rho))+((I**6)*P_A_sym_6(rho))
        return(P)
    if order_approx==0.0:
        P=P_general(rho,I)
        return(P)
    
    



# Now define the beta equilibrium pressure with the presence of muon(P_beta_with_muon):

def P_beta_with_muon(rho,I,order_approx):
    P=P_bary(rho,I,order_approx)+P_e(rho,I,order_approx)+P_mu(rho,I,order_approx) 
    return(P)
    
    



# Now define the beta equilibrium pressure in the absence of muon(P_beta_without_muon):  

def P_beta_without_muon(rho,I,order_approx):
    P=P_bary(rho,I,order_approx)+P_e(rho,I,order_approx) 
    return(P)   





# Calculation of isospin asymmetry(I):

def beta_eos(rho,order_approx):
    I_start=1.0
    rho_old=rho
    for j in range(100000000):
        if mu_e(rho_old,I_start,order_approx)>m_mu_c_2:
            mue=mu_e(rho_old,I_start,order_approx)
            mumu=mue
            Pe=P_e(rho_old,I_start,order_approx)
            Pmu=P_mu(rho_old,I_start,order_approx)
            rhoe=rho_e(rho_old,I_start,order_approx)
            rhomu=rho_mu(rho_old,I_start,order_approx)
            epsilone=epsilon_e(rho_old,I_start,order_approx)
            epsilonmu=epsilon_mu(rho_old,I_start,order_approx)
            rhop=rhoe+rhomu
            rhon=rho_old-rhop
            I_next=(rhon-rhop)/(rhon+rhop)
            Pb=P_bary(rho_old,I_start,order_approx)
            Pbeta=P_beta_with_muon(rho_old,I_start,order_approx)
            if abs(I_next-I_start)<1.e-6:
                return(I_start,rhop,rhon,rhoe,rhomu,Pe,Pmu,epsilone,epsilonmu,mue,mumu,Pb,Pbeta)
                break
            else:
                I_start=(0.8*I_start)+(0.2*I_next)
                j=j+1
        else:
            I0=I_start
            mue=mu_e(rho_old,I_start,order_approx)
            mumu=0.0
            Pe=P_e(rho_old,I_start,order_approx)
            Pmu=0.0
            rhoe=rho_e(rho_old,I_start,order_approx)
            rhomu=0.0
            epsilone=epsilon_e(rho_old,I_start,order_approx)
            epsilonmu=0.0
            rhop=rho_e(rho_old,I_start,order_approx)
            Pbeta=P_beta_without_muon(rho_old,I_start,order_approx)
            Pb=P_bary(rho_old,I_start,order_approx)
            rhon=rho_old-rhop
            I_next=(rhon-rhop)/(rhon+rhop)
            if abs(I_next-I_start)<1.e-6:
                return(I_start,rhop,rhon,rhoe,rhomu,Pe,Pmu,epsilone,epsilonmu,mue,mumu,Pb,Pbeta)
                break
            else:
                I_start=(0.8*I_start)+(0.2*I_next)
                j=j+1                
