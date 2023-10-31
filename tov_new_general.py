import numpy as np

import matplotlib.pyplot as plt

import h5py

def TOV_SOLVER(n0):

    from scipy.interpolate import CubicSpline

    from skyrm2 import beta_eos

    from skyrm2 import e

    pi=3.141592653589793

    n=np.linspace(0.0001,1.92,1000)

    I=[]

    P_beta=[]


    for i in range(len(n)):
         I0_old,rhop0,rhon0,rhoe0,rhomu0,Pe0,Pmu0,epsilone0,epsilonmu0,mue0,mumu0,Pb0,Pbeta0=beta_eos(n[i],1.0)
         I.append(I0_old)
         P_beta.append(Pbeta0)
    
    f=CubicSpline(n,I)

    g=CubicSpline(P_beta,n)

    G=6.67e-08                                                                                                      # Gravitational wave constant value 
  
    C=3.e10                                                                                                         # velocity of light
  
    M_sun=1.989e33                                                                                                  # solar mass

    I0_old,rhop0,rhon0,rhoe0,rhomu0,Pe0,Pmu0,epsilone0,epsilonmu0,mue0,mumu0,Pb0,Pbeta0=beta_eos(n0,1.0)

    e0_old=e(n0,I0_old)

    e0=(e0_old*(10**39)*G)/(625000*(C**4))

    r0=100
   
    rf=2500000
   
    n=100001                                                                                              
  
    h=(rf-r0)/(n-1) 
   
    r=np.linspace(r0,rf,n)
   
    m=np.zeros([n])            
   
#rho=np.zeros([n])
   
    p=np.zeros([n])

    m[0]=(4*pi*(r0**3)*e0)/3.0

    p[0]=(Pbeta0*(10**39)*G)/(625000*(C**4))
   
    e1=np.zeros([n])

    e1[0]=e0

    def dpdr(r,P,m):
         P1=(P*625000*(C**4))/((10**39)*G)
         n=g(P1)
         I=f(n)
         e0=e(n,I)
         e0_new=(e0*(10**39)*G)/(625000*(C**4))
         T=P+e0_new
         W=m+(4*pi*(r**3)*P)
         V=(2*m*r)-(r**2)
         Z=(T*W)/V
         return(Z)
  
  
    def dmdr(r,P,m):
         P1=(P*625000*(C**4))/((10**39)*G)
         n=g(P1)
         I=f(n)
         e0=e(n,I)
         e0_new=(e0*(10**39)*G)/(625000*(C**4))
         H=e0_new*(4*pi*(r**2))
         return(H)
  
    for i in range(1,n):
 
  
         if(p[i-1]>=1.e-20):
           k1=h*dpdr(r[i-1],p[i-1],m[i-1])
           l1=h*dmdr(r[i-1],p[i-1],m[i-1])
           k2=h*dpdr((r[i-1]+(h/2)),(p[i-1]+(k1/2)),(m[i-1]+(l1/2)))
           l2=h*dmdr((r[i-1]+(h/2)),(p[i-1]+(k1/2)),(m[i-1]+(l1/2)))
           k3=h*dpdr((r[i-1]+(h/2)),(p[i-1]+(k2/2)),(m[i-1]+(l2/2)))
           l3=h*dmdr((r[i-1]+(h/2)),(p[i-1]+(k2/2)),(m[i-1]+(l2/2)))
           k4=h*dpdr((r[i-1]+(h)),(p[i-1]+(k3)),(m[i-1]+(l3)))
           l4=h*dmdr((r[i-1]+(h)),(p[i-1]+(k3)),(m[i-1]+(l3)))
           p[i]=p[i-1]+(((k1+(2*k2)+(2*k3)+k4)/6))
           m[i]=m[i-1]+(((l1+(2*l2)+(2*l3)+l4)/6))
           p_I=(p[i]*625000*(C**4))/((10**39)*G)
           n_I=g(p_I)
           I_I=f(n_I)
           e_I=e(n_I,I_I)
           e1[i]=(e_I*(10**39)*G)/(625000*(C**4))
         else:
           M=(m[i-2]*(C**2))/((G)*(M_sun))
           #print(M)
           R=(r[i-2]*1.e-05)
           P_max=(p[0]*(C**4))/G
           E_max=(e1[0]*(C**4))/G
           rho_max=(e1[0]*(C**4))/(G*(C**2))
           break
    return(M,R,P_max,E_max,rho_max)
