import numpy as np
import matplotlib.pyplot as plt
import h5py


h_c=197.32698                                                         # Mev-fm unit(it is actually planck constant times speed of light ).
m_n_c_2=939.0                                                         # Mev unit(it is equivalent mass energy of nucleon).
pi=3.141592653589793
pi2=pi**2
m_e_c_2=0.511                                                         # Mev unit(it is equivalent mass energy of electron)      
m_mu_c_2=105.6583755                                                  # Mev unit(it is equivalent mass energy of muon)



def F(m,I):
    a=((1.0+I)**m)+((1.0-I)**m)
    X=a/2.0
    return(X)


def T_k(rho,I):
    a=(h_c**2)/(2*m_n_c_2)
    b=(3*pi2)/2
    m=5/3
    c=(rho**(5/3))*F(m,I)
    d=(3*a*(b**(2/3))*c)/5
    return(d)
    

# Now I define the zero range term(H0):

def H0(rho,I,x0,t0):
    m=2
    a=(2*(2+x0))-((2*x0+1)*F(m,I))
    b=(t0*(rho**2)*a)/8
    return(b)
   
   
# Now I define the density dependent term(H31):


def H31(rho,I,x31,t31,sigma1):
    m=2
    a=(2*(2+x31))-((2*x31+1)*F(m,I))
    b=(t31*(rho**(sigma1+2))*a)/48
    return(b)
    
    
def H32(rho,I,x32,t32,sigma2):
    m=2
    a=(2*(2+x32))-((2*x32+1)*F(m,I))
    b=(t32*(rho**(sigma2+2))*a)/48
    return(b)
    
def H33(rho,I,x33,t33,sigma3):
    m=2
    a=(2*(2+x33))-((2*x33+1)*F(m,I))
    b=(t33*(rho**(sigma3+2))*a)/48
    return(b)
    
    
# Now I define the effective mass term(Heff):


def Heff(rho,I,t1,t2,x1,x2):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    func1=(t1*(2+x1))+(t2*(2+x2))
    func2=(t2*(2*x2+1))-(t1*(2*x1+1))
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**(8/3))*func3)/40
    return(d)
    
    
def Heff2(rho,I,t4,x4,delta):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(8.0/3.0)+delta
    func1=t4*(x4+2)
    func2=t4*((2*x4)+1)
    func3=(func1*F(m1,I))-((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**(n))*func3)/40
    return(d)
    
def Heff3(rho,I,t5,x5,gamma):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(8.0/3.0)+gamma
    func1=t5*(x5+2)
    func2=t5*((2*x5)+1)
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**(n))*func3)/40
    return(d)    

# Now above I define all components of energy density T_k, H0, H3, Heff. So now I define total energy density e:


def e(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=T_k(rho,I)
    a2=H0(rho,I,x0,t0)
    a3=H31(rho,I,x31,t31,sigma1)
    a4=H32(rho,I,x32,t32,sigma2)
    a5=H33(rho,I,x33,t33,sigma3)
    a6=Heff(rho,I,t1,t2,x1,x2)
    a7=Heff2(rho,I,t4,x4,delta)
    a8=Heff3(rho,I,t5,x5,gamma)
    add=a1+a2+a3+a4+a5+a6+a7+a8
    return(add)
    
    
# Now I define energy per nucleon term(E/A):


def E_A(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    E_A=e(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)/rho
    return(E_A)
    
def Y_p(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    y=(1.0-I)/2.0
    return(y)
    
    
def debyrho_drho_T_k(rho,I):
    A=(3*(h_c**2))/(10.0*m_n_c_2)
    B=(3*pi2)/2.0
    m1=5.0/3.0
    term1=(2*(rho**(-1/3))*F(m1,I))/3
    term=A*(B**(2/3))*term1
    return(term)
    
    
def debyrho_drho_H0(rho,I,x0,t0):
    m1=2
    term1=((2*(x0+2))-(((2*x0)+1)*F(m1,I)))
    term=(t0*term1)/8
    return(term)
    

def debyrho_drho_H31(rho,I,x31,t31,sigma1):
    m1=2
    term1=(sigma1+1)*(rho**sigma1)*((2*(x31+2))-(((2*x31)+1)*F(m1,I)))
    term=(t31*term1)/48
    return(term)   
    

def debyrho_drho_H32(rho,I,x32,t32,sigma2):
    m1=2
    term1=(sigma2+1)*(rho**sigma2)*((2*(x32+2))-(((2*x32)+1)*F(m1,I)))
    term=(t32*term1)/48
    return(term)   
    
    
def debyrho_drho_H33(rho,I,x33,t33,sigma3):
    m1=2
    term1=(sigma3+1)*(rho**sigma3)*((2*(x33+2))-(((2*x33)+1)*F(m1,I)))
    term=(t33*term1)/48
    return(term) 
    
    
def debyrho_drho_Heff(rho,I,t1,t2,x1,x2):
    m1=5/3
    m2=8/3
    B=(3*pi2)/2.0
    func1=((t1*(x1+2))+(t2*(x2+2)))*F(m1,I)
    func2=((t2*((2*x2)+1))-(t1*((2*x1)+1)))*F(m2,I)
    term1=(5*(rho**(2/3))*(func1+(func2/2)))/3
    term=(3*(B**(2/3))*term1)/40
    return(term)
    
def debyrho_drho_Heff2(rho,I,t4,x4,delta):    
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+delta
    func1=t4*(x4+2)
    func2=t4*((2*x4)+1)
    func3=(func1*F(m1,I))-((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(rho**(n-1.0))*func3)/40
    return(d)


def debyrho_drho_Heff3(rho,I,t5,x5,gamma):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+gamma
    func1=t5*(x5+2)
    func2=t5*((2*x5)+1)
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(rho**(n-1.0))*func3)/40
    return(d)
    

def debyrho_drho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=debyrho_drho_T_k(rho,I)
    a2=debyrho_drho_H0(rho,I,x0,t0)
    a3=debyrho_drho_H31(rho,I,x31,t31,sigma1)
    a4=debyrho_drho_H32(rho,I,x32,t32,sigma2)
    a5=debyrho_drho_H33(rho,I,x33,t33,sigma3)
    a6=debyrho_drho_Heff(rho,I,t1,t2,x1,x2)
    a7=debyrho_drho_Heff2(rho,I,t4,x4,delta)
    a8=debyrho_drho_Heff3(rho,I,t5,x5,gamma)
    add=a1+a2+a3+a4+a5+a6+a7+a8
    return(add)  

    
def P_T_k(rho,I):
    A=(3*(h_c**2))/(10.0*m_n_c_2)
    B=(3*pi2)/2.0
    m1=5.0/3.0
    term1=(2*(rho**(-1/3))*F(m1,I))/3
    term=(rho**2)*A*(B**(2/3))*term1
    return(term)
    
    
def P_H0(rho,I,x0,t0):
    m1=2
    term1=((2*(x0+2))-(((2*x0)+1)*F(m1,I)))
    term=((rho**2)*t0*term1)/8
    return(term)
    
    
def P_H31(rho,I,x31,t31,sigma1):
    m1=2
    term1=(sigma1+1)*(rho**sigma1)*((2*(x31+2))-(((2*x31)+1)*F(m1,I)))
    term=((rho**2)*t31*term1)/48
    return(term)
    
def P_H32(rho,I,x32,t32,sigma2):
    m1=2
    term1=(sigma2+1)*(rho**sigma2)*((2*(x32+2))-(((2*x32)+1)*F(m1,I)))
    term=((rho**2)*t32*term1)/48
    return(term)   
    
    
def P_H33(rho,I,x33,t33,sigma3):
    m1=2
    term1=(sigma3+1)*(rho**sigma3)*((2*(x33+2))-(((2*x33)+1)*F(m1,I)))
    term=((rho**2)*t33*term1)/48
    return(term) 
    
def P_Heff(rho,I,t1,t2,x1,x2):
    m1=5/3
    m2=8/3
    B=(3*pi2)/2.0
    func1=((t1*(x1+2))+(t2*(x2+2)))*F(m1,I)
    func2=((t2*((2*x2)+1))-(t1*((2*x1)+1)))*F(m2,I)
    term1=(5*(rho**(2/3))*(func1+(func2/2)))/3
    term=((rho**2)*3*(B**(2/3))*term1)/40
    return(term)
    
def P_Heff2(rho,I,t4,x4,delta):    
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+delta
    func1=t4*(x4+2)
    func2=t4*((2*x4)+1)
    func3=(func1*F(m1,I))-((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**2)*n*(rho**(n-1.0))*func3)/40
    return(d)


def P_Heff3(rho,I,t5,x5,gamma):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+gamma
    func1=t5*(x5+2)
    func2=t5*((2*x5)+1)
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*(rho**2)*n*(rho**(n-1.0))*func3)/40
    return(d)  
    
def P(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=P_T_k(rho,I)
    a2=P_H0(rho,I,x0,t0)
    a3=P_H31(rho,I,x31,t31,sigma1)
    a4=P_H32(rho,I,x32,t32,sigma2)
    a5=P_H33(rho,I,x33,t33,sigma3)
    a6=P_Heff(rho,I,t1,t2,x1,x2)
    a7=P_Heff2(rho,I,t4,x4,delta)
    a8=P_Heff3(rho,I,t5,x5,gamma)
    add=a1+a2+a3+a4+a5+a6+a7+a8
    return(add)
    
def d2ebyrho_d2rho_T_k(rho,I):
    A=(3*(h_c**2))/(10.0*m_n_c_2)
    B=(3*pi2)/2.0
    m1=5.0/3.0
    term1=-(2*(rho**(-4/3))*F(m1,I))/9
    term=A*(B**(2/3))*term1
    return(term)
    
def d2ebyrho_d2rho_H31(rho,I,x31,t31,sigma1):
    m1=2
    term1=(sigma1+1)*sigma1*(rho**(sigma1-1))*((2*(x31+2))-(((2*x31)+1)*F(m1,I)))
    term=(t31*term1)/48
    return(term) 
    
def d2ebyrho_d2rho_H32(rho,I,x32,t32,sigma2):   
    m1=2
    term1=(sigma2+1)*sigma2*(rho**(sigma2-1))*((2*(x32+2))-(((2*x32)+1)*F(m1,I)))
    term=(t32*term1)/48
    return(term)
    
def d2ebyrho_d2rho_H33(rho,I,x33,t33,sigma3):
    m1=2
    term1=(sigma3+1)*sigma3*(rho**(sigma3-1))*((2*(x33+2))-(((2*x33)+1)*F(m1,I)))
    term=(t33*term1)/48
    return(term)

def d2ebyrho_d2rho_Heff(rho,I,t1,t2,x1,x2):
    m1=5/3
    m2=8/3
    B=(3*pi2)/2.0
    func1=((t1*(x1+2))+(t2*(x2+2)))*F(m1,I)
    func2=((t2*((2*x2)+1))-(t1*((2*x1)+1)))*F(m2,I)
    term1=(5*2*(rho**(-1/3))*(func1+(func2/2)))/9
    term=(3*(B**(2/3))*term1)/40
    return(term)
    
def d2ebyrho_d2rho_Heff2(rho,I,t4,x4,delta):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+delta
    func1=t4*(x4+2)
    func2=t4*((2*x4)+1)
    func3=(func1*F(m1,I))-((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(n-1.0)*(rho**(n-2.0))*func3)/40
    return(d)
    

def d2ebyrho_d2rho_Heff3(rho,I,t5,x5,gamma):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+gamma
    func1=t5*(x5+2)
    func2=t5*((2*x5)+1)
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(n-1.0)*(rho**(n-2.0))*func3)/40
    return(d) 

def d2ebyrho_d2rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=d2ebyrho_d2rho_T_k(rho,I)
    a3=d2ebyrho_d2rho_H31(rho,I,x31,t31,sigma1)
    a4=d2ebyrho_d2rho_H32(rho,I,x32,t32,sigma2)
    a5=d2ebyrho_d2rho_H33(rho,I,x33,t33,sigma3)
    a6=d2ebyrho_d2rho_Heff(rho,I,t1,t2,x1,x2)
    a7=d2ebyrho_d2rho_Heff2(rho,I,t4,x4,delta)
    a8=d2ebyrho_d2rho_Heff3(rho,I,t5,x5,gamma)
    add=a1+a3+a4+a5+a6+a7+a8
    return(add)  
    
def K(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    add=9*(rho**2)*d2ebyrho_d2rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)
    return(add)

def rho_sat(rho0,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    rho_select=np.zeros([1000001])
    rho_select[0]=rho0
    rho_sat=[]
    for i in range(1,len(rho_select)):
          rho_select[i]=rho_select[i-1]-(debyrho_drho(rho_select[i-1],0.0,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)/d2ebyrho_d2rho(rho_select[i-1],0.0,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma))
          if(abs(rho_select[i]-rho_select[i-1])<0.0000001):
               return(rho_select[i-1])
               

def d3ebyrho_d3rho_T_k(rho,I):
    A=(3*(h_c**2))/(10.0*m_n_c_2)
    B=(3*pi2)/2.0
    m1=5.0/3.0
    term1=(2*4*(rho**(-7/3))*F(m1,I))/(9*3)
    term=A*(B**(2/3))*term1
    return(term) 
    
def d3ebyrho_d3rho_H31(rho,I,x31,t31,sigma1):
    m1=2
    term1=(sigma1+1)*sigma1*(sigma1-1)*(rho**(sigma1-2))*((2*(x31+2))-(((2*x31)+1)*F(m1,I)))
    term=(t31*term1)/48
    return(term) 
    
def d3ebyrho_d3rho_H32(rho,I,x32,t32,sigma2):   
    m1=2
    term1=(sigma2+1)*sigma2*(sigma2-1)*(rho**(sigma2-2))*((2*(x32+2))-(((2*x32)+1)*F(m1,I)))
    term=(t32*term1)/48
    return(term)
    
def d3ebyrho_d3rho_H33(rho,I,x33,t33,sigma3):
    m1=2
    term1=(sigma3+1)*sigma3*(sigma3-1)*(rho**(sigma3-2))*((2*(x33+2))-(((2*x33)+1)*F(m1,I)))
    term=(t33*term1)/48
    return(term)    
    
def d3ebyrho_d3rho_Heff(rho,I,t1,t2,x1,x2):
    m1=5/3
    m2=8/3
    B=(3*pi2)/2.0
    func1=((t1*(x1+2))+(t2*(x2+2)))*F(m1,I)
    func2=((t2*((2*x2)+1))-(t1*((2*x1)+1)))*F(m2,I)
    term1=-(5*2*(rho**(-4/3))*(func1+(func2/2)))/(9*3)
    term=(3*(B**(2/3))*term1)/40
    return(term)
    
def d3ebyrho_d3rho_Heff2(rho,I,t4,x4,delta):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+delta
    func1=t4*(x4+2)
    func2=t4*((2*x4)+1)
    func3=(func1*F(m1,I))-((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(n-1.0)*(n-2.0)*(rho**(n-3.0))*func3)/40
    return(d)
    
def d3ebyrho_d3rho_Heff3(rho,I,t5,x5,gamma):
    m1=5/3
    m2=8/3
    b=(3*pi2)/2
    n=(5.0/3.0)+gamma
    func1=t5*(x5+2)
    func2=t5*((2*x5)+1)
    func3=(func1*F(m1,I))+((func2*F(m2,I))/2)
    d=(3*(b**(2/3))*n*(n-1.0)*(n-2.0)*(rho**(n-3.0))*func3)/40
    return(d) 
    
def d3ebyrho_d3rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=d3ebyrho_d3rho_T_k(rho,I)
    a3=d3ebyrho_d3rho_H31(rho,I,x31,t31,sigma1)
    a4=d3ebyrho_d3rho_H32(rho,I,x32,t32,sigma2)
    a5=d3ebyrho_d3rho_H33(rho,I,x33,t33,sigma3)
    a6=d3ebyrho_d3rho_Heff(rho,I,t1,t2,x1,x2)
    a7=d3ebyrho_d3rho_Heff2(rho,I,t4,x4,delta)
    a8=d3ebyrho_d3rho_Heff3(rho,I,t5,x5,gamma)
    add=a1+a3+a4+a5+a6+a7+a8
    return(add) 
    
def Q(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    add=27.0*(rho**3)*d3ebyrho_d3rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)
    return(add)  
    
def E_sym_T_k(rho,I):  
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    term1=A*(B**(2/3))*(rho**(2.0/3.0))
    return(term1)
    
def E_sym_H0(rho,I,x0,t0):
    func1=((2*x0)+1)
    term2=(t0*func1*rho)/8.0
    return(term2)
    
def E_sym_H31(rho,I,x31,t31,sigma1):
    func2=((2*x31)+1)
    term3=(t31*func2*(rho**(sigma1+1)))/48.0
    return(term3)
    
    
def E_sym_H32(rho,I,x32,t32,sigma2):
    func2=((2*x32)+1)
    term3=(t32*func2*(rho**(sigma2+1)))/48.0
    return(term3)
    
def E_sym_H33(rho,I,x33,t33,sigma3):
    func2=((2*x33)+1)
    term3=(t33*func2*(rho**(sigma3+1)))/48.0
    return(term3)
    
def E_sym_Heff(rho,I,t1,t2,x1,x2):
    B=(3*pi2)/2.0
    func3=-(3*t1*x1)+(4*t2)+(5*t2*x2)
    term4=((B**(2.0/3.0))*func3*(rho**(5.0/3.0)))/24.0
    return(term4)
    
def E_sym_Heff2(rho,I,t4,x4,delta):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+delta
    a=t4*x4*(rho**n)*(B**(2.0/3.0))
    term5=a/8.0
    return(term5)
    
def E_sym_Heff3(rho,I,t5,x5,gamma):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+gamma
    a=t5*((5*x5)+4)*(rho**n)*(B**(2.0/3.0))
    term6=a/24.0
    return(term6)
    
    
def E_sym(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=E_sym_T_k(rho,I)
    a2=E_sym_H0(rho,I,x0,t0)
    a3=E_sym_H31(rho,I,x31,t31,sigma1)
    a4=E_sym_H32(rho,I,x32,t32,sigma2)
    a5=E_sym_H33(rho,I,x33,t33,sigma3)
    a6=E_sym_Heff(rho,I,t1,t2,x1,x2)
    a7=E_sym_Heff2(rho,I,t4,x4,delta)
    a8=E_sym_Heff3(rho,I,t5,x5,gamma)
    add=a1-a2-a3-a4-a5+a6-a7-a8
    return(add)
    
    
def dE_sym_drho_T_k(rho,I):
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    term1=(2*A*(B**(2/3))*(rho**(-1.0/3.0)))/3.0
    return(term1)
    
  
def dE_sym_drho_H0(rho,I,x0,t0):
    func1=((2*x0)+1)
    term2=(t0*func1)/8.0
    return(term2)
    
def dE_sym_drho_H31(rho,I,x31,t31,sigma1):
    func2=((2*x31)+1)
    term3=(t31*func2*(sigma1+1)*(rho**(sigma1)))/48.0
    return(term3)
    
def dE_sym_drho_H32(rho,I,x32,t32,sigma2):
    func2=((2*x32)+1)
    term3=(t32*func2*(sigma2+1)*(rho**(sigma2)))/48.0
    return(term3)
    
def dE_sym_drho_H33(rho,I,x33,t33,sigma3):
    func2=((2*x33)+1)
    term3=(t33*func2*(sigma3+1)*(rho**(sigma3)))/48.0
    return(term3) 
    
def dE_sym_drho_Heff(rho,I,t1,t2,x1,x2):
    B=(3*pi2)/2.0
    func3=-(3*t1*x1)+(4*t2)+(5*t2*x2)
    term4=(5*(B**(2.0/3.0))*func3*(rho**(2.0/3.0)))/72.0
    return(term4)  
    
def dE_sym_drho_Heff2(rho,I,t4,x4,delta):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+delta
    a=t4*x4*n*(rho**(n-1))*(B**(2.0/3.0))
    term5=a/8.0
    return(term5) 
    
def dE_sym_drho_Heff3(rho,I,t5,x5,gamma):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+gamma
    a=t5*((5*x5)+4)*n*(rho**(n-1))*(B**(2.0/3.0))
    term6=a/24.0
    return(term6)

def dE_sym_drho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=dE_sym_drho_T_k(rho,I)
    a2=dE_sym_drho_H0(rho,I,x0,t0)
    a3=dE_sym_drho_H31(rho,I,x31,t31,sigma1)
    a4=dE_sym_drho_H32(rho,I,x32,t32,sigma2)
    a5=dE_sym_drho_H33(rho,I,x33,t33,sigma3)
    a6=dE_sym_drho_Heff(rho,I,t1,t2,x1,x2)
    a7=dE_sym_drho_Heff2(rho,I,t4,x4,delta)
    a8=dE_sym_drho_Heff3(rho,I,t5,x5,gamma)
    add=a1-a2-a3-a4-a5+a6-a7-a8
    return(add)
    

def L(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a=3*rho*dE_sym_drho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)
    return(a)
    
    
def d2E_sym_d2rho_T_k(rho,I):
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    term1=-(2*A*(B**(2/3))*(rho**(-4.0/3.0)))/9.0
    return(term1)
    
def d2E_sym_d2rho_H31(rho,I,x31,t31,sigma1):
    func2=((2*x31)+1)
    term3=(t31*func2*(sigma1+1)*sigma1*(rho**(sigma1-1)))/48.0
    return(term3)
    
def d2E_sym_d2rho_H32(rho,I,x32,t32,sigma2):
    func2=((2*x32)+1)
    term3=(t32*func2*(sigma2+1)*sigma2*(rho**(sigma2-1)))/48.0
    return(term3)
    
def d2E_sym_d2rho_H33(rho,I,x33,t33,sigma3):
    func2=((2*x33)+1)
    term3=(t33*func2*(sigma3+1)*sigma3*(rho**(sigma3-1)))/48.0
    return(term3)
    
def d2E_sym_d2rho_Heff(rho,I,t1,t2,x1,x2):
    B=(3*pi2)/2.0
    func3=-(3*t1*x1)+(4*t2)+(5*t2*x2)
    term4=(5*2*(B**(2.0/3.0))*func3*(rho**(-1.0/3.0)))/(72.0*3)
    return(term4)
    
def d2E_sym_d2rho_Heff2(rho,I,t4,x4,delta):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+delta
    a=t4*x4*n*(n-1.0)*(rho**(n-2))*(B**(2.0/3.0))
    term5=a/8.0
    return(term5) 
    
def d2E_sym_d2rho_Heff3(rho,I,t5,x5,gamma):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+gamma
    a=t5*((5*x5)+4)*n*(n-1.0)*(rho**(n-2))*(B**(2.0/3.0))
    term6=a/24.0
    return(term6)

def d2E_sym_d2rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=d2E_sym_d2rho_T_k(rho,I)
    a3=d2E_sym_d2rho_H31(rho,I,x31,t31,sigma1)
    a4=d2E_sym_d2rho_H32(rho,I,x32,t32,sigma2)
    a5=d2E_sym_d2rho_H33(rho,I,x33,t33,sigma3)
    a6=d2E_sym_d2rho_Heff(rho,I,t1,t2,x1,x2)
    a7=d2E_sym_d2rho_Heff2(rho,I,t4,x4,delta)
    a8=d2E_sym_d2rho_Heff3(rho,I,t5,x5,gamma)
    add=a1-a3-a4-a5+a6-a7-a8
    return(add)
    
    
def Ksym(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a=9*(rho**2)*d2E_sym_d2rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)
    return(a)
    
def d3E_sym_d3rho_T_k(rho,I):
    A=(h_c**2.0)/(6.0*m_n_c_2)
    B=(3*pi2)/2.0
    term1=(2*4*A*(B**(2/3))*(rho**(-7.0/3.0)))/(9.0*3.0)
    return(term1)
    
def d3E_sym_d3rho_H31(rho,I,x31,t31,sigma1):
    func2=((2*x31)+1)
    term3=(t31*func2*(sigma1+1)*sigma1*(sigma1-1)*(rho**(sigma1-2)))/48.0
    return(term3)
    
def d3E_sym_d3rho_H32(rho,I,x32,t32,sigma2):
    func2=((2*x32)+1)
    term3=(t32*func2*(sigma2+1)*sigma2*(sigma2-1)*(rho**(sigma2-2)))/48.0
    return(term3)
    
def d3E_sym_d3rho_H33(rho,I,x33,t33,sigma3):
    func2=((2*x33)+1)
    term3=(t33*func2*(sigma3+1)*sigma3*(sigma3-1)*(rho**(sigma3-2)))/48.0
    return(term3)
    
def d3E_sym_d3rho_Heff(rho,I,t1,t2,x1,x2):
    B=(3*pi2)/2.0
    func3=-(3*t1*x1)+(4*t2)+(5*t2*x2)
    term4=-(5*2*(B**(2.0/3.0))*func3*(rho**(-4.0/3.0)))/(72.0*3*3)
    return(term4)
    
def d3E_sym_d3rho_Heff2(rho,I,t4,x4,delta):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+delta
    a=t4*x4*n*(n-1.0)*(n-2.0)*(rho**(n-3))*(B**(2.0/3.0))
    term5=a/8.0
    return(term5) 
    
def d3E_sym_d3rho_Heff3(rho,I,t5,x5,gamma):
    B=(3*pi2)/2.0
    n=(5.0/3.0)+gamma
    a=t5*((5*x5)+4)*n*(n-1.0)*(n-2.0)*(rho**(n-3))*(B**(2.0/3.0))
    term6=a/24.0
    return(term6)
    
def d3E_sym_d3rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a1=d3E_sym_d3rho_T_k(rho,I)
    a3=d3E_sym_d3rho_H31(rho,I,x31,t31,sigma1)
    a4=d3E_sym_d3rho_H32(rho,I,x32,t32,sigma2)
    a5=d3E_sym_d3rho_H33(rho,I,x33,t33,sigma3)
    a6=d3E_sym_d3rho_Heff(rho,I,t1,t2,x1,x2)
    a7=d3E_sym_d3rho_Heff2(rho,I,t4,x4,delta)
    a8=d3E_sym_d3rho_Heff3(rho,I,t5,x5,gamma)
    add=a1-a3-a4-a5+a6-a7-a8
    return(add)
    
def Qsym(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma):
    a=27*(rho**3)*d3E_sym_d3rho(rho,I,t0,t1,t2,t31,t32,t33,t4,t5,x0,x1,x2,x31,x32,x33,x4,x5,sigma1,sigma2,sigma3,delta,gamma)
    return(a)
    
    

