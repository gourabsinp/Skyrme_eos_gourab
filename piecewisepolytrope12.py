import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import sys


def Tov_solver(rho0):
  
  
  
  r0=10                                                                                                          # initial value of radial point 
  
  h=(2500000-10)/100000                                                                                          # step by which increase the radial point
  r=np.linspace(10,2500000,100001)
  m=np.zeros([100001])            
  rho=np.zeros([100001])
  p=np.zeros([100001])
  e=np.zeros([100001])
  
  G=6.67e-08                                                                                                      # Gravitational wave constant value 
  
  C=3.e10                                                                                                         # velocity of light
  
  M_sun=1.989e33                                                                                                  # solar mass
  
  #rho_lower1=((10**3)*G)/(C**2)
  
  rho_lower1=0.0
  
  rho_lower_2=((2.44034e+07)*G)/(C**2)
  
  rho_lower_3=((3.78358e+11)*G)/(C**2)
  
  rho_lower_4=((2.62780e+12)*G)/(C**2)
  
  rho_lower_5=((2.7e+14)*G)/(C**2)
  
  p0=((10**34.2)*G)/(C**4)
  
  p1=((10**34.384)*G)/(C**4)                                                                                      # the value of pressure at rest mass density rho1
  
  rho1=((10**14.7)*G)/(C**2)                                                                                      # First transition rest mass density.
  
  rho2=((10**15)*G)/(C**2)                                                                                        # second transition rest mass density.
  
  gama1_lower=1.58425
  
  gama2_lower=1.28733
  
  gama3_lower=0.62223
  
  gama4_lower=1.35692
  
  gama1=3.005                                                                                                     # value of polytropic index(rho<=rho1)
  
  gama2=2.988                                                                                                     # value of polytropic index(rho1<rho<=rho2)
  
  gama3=2.851                                                                                                  # value of polytropic index(rho>rho2)
  
  k1_lower=(6.80110e-09)/((G**(gama1_lower-1))*(C**(4-2*gama1_lower)))
  
  k2_lower=(1.06186e-06)/((G**(gama2_lower-1))*(C**(4-2*gama2_lower)))
  
  k3_lower=(5.32697e+01)/((G**(gama3_lower-1))*(C**(4-2*gama3_lower)))
  
  k4_lower=(3.99874e-08)/((G**(gama4_lower-1))*(C**(4-2*gama4_lower)))
  
  print("k1=",k1_lower)
  
  print("k2=",k2_lower)
  
  print("k3=",k3_lower)
  
  print("k4=",k4_lower)
  
  #k_x=p0/(rho_lower_5**gama1)   
  
  #k_x=(k4_lower*(rho_lower_5**gama4_lower))/(rho_lower_5**gama1)                                                                                        # value of polytropic constant(rho<=rho1)
  
  
  
  #k_x=k4_lower*((rho_lower_5)**(gama4_lower-gama1))
  
  k_x=p1/(rho1**gama1)
  
  p1_compare=k_x*(rho1**gama1)
  
  k_y=p1/(rho1**gama2)                                                                                           # value of polytropic constant(rho1<rho<=rho2)
  
  k_y_compare=k_x*(rho1**(gama1-gama2))
  
  #k_y=k_x*(rho1**(gama1-gama2))
  
  k_z=k_y*(rho2**(gama2-gama3))                                                                               # value of polytropic constant(rho>rho2)
  #print("p_0",k3_lower*(rho_lower_4**gama3_lower))
  #print("p0=",k4_lower*(rho_lower_5**gama4_lower))
  #print("p1=",p1)
  #print("p1_compare=",p1_compare)
  #print("p2=",k_y*(rho2**gama2))
  #print(k_y)
  #print(k_y_compare)
  
  
  
  a1_lower=0
  
  a2_lower=0
  
  a3_lower=0
  
  a4_lower=0
  
  
  #a2_lower=((k1_lower/(gama1_lower-1))*(rho_lower_2**(gama1_lower-1)))-((k2_lower/(gama2_lower-1))*(rho_lower_2**(gama2_lower-1)))
  
  #a3_lower=a2_lower+((k2_lower/(gama2_lower-1))*(rho_lower_3**(gama2_lower-1)))-((k3_lower/(gama3_lower-1))*(rho_lower_3**(gama3_lower-1)))
  
  #a4_lower=a3_lower+((k3_lower/(gama3_lower-1))*(rho_lower_4**(gama3_lower-1)))-((k4_lower/(gama4_lower-1))*(rho_lower_4**(gama4_lower-1)))
  
  #a_1=a4_lower+((k4_lower/(gama4_lower-1))*(rho_lower_5**(gama4_lower-1)))-((k_x/(gama1-1))*(rho_lower_5**(gama1-1)))
  
  a_1=0
  
  a_2=a_1+((k_x/(gama1-1))*(rho1**(gama1-1)))-((k_y/(gama2-1))*(rho1**(gama2-1)))
  
  a_3=a_2+((k_y/(gama2-1))*(rho2**(gama2-1)))-((k_z/(gama3-1))*(rho2**(gama3-1)))
  
  
  # The condition use to know the value of k & gama for initial/central rest mass density & energy density.
  
  if(rho0<rho_lower_2):
   k=k1_lower
   gama=gama1_lower
   a=a1_lower
   e0=rho0+(a1_lower*rho0)+((k1_lower*(rho0**gama1_lower))/(gama1_lower-1))
     
  if(rho_lower_2<=rho0<rho_lower_3):
   k=k2_lower
   gama=gama2_lower
   a=a2_lower
   e0=rho0+(a2_lower*rho0)+((k2_lower*(rho0**gama2_lower))/(gama2_lower-1))
     
  if(rho_lower_3<=rho0<rho_lower_4):
   k=k3_lower
   gama=gama3_lower
   a=a3_lower
   e0=rho0+(a3_lower*rho0)+((k3_lower*(rho0**gama3_lower))/(gama3_lower-1))
     
  if(rho_lower_4<=rho0<rho_lower_5):
   k=k4_lower
   gama=gama4_lower
   a=a4_lower
   e0=rho0+(a4_lower*rho0)+((k4_lower*(rho0**gama4_lower))/(gama4_lower-1))
  
  
  if(rho_lower_5<=rho0<rho1):
   k=k_x
   gama=gama1
   a=a_1
   e0=rho0+(a_1*rho0)+((k_x*(rho0**gama1))/(gama1-1))
  
  if(rho1<=rho0<rho2):
   k=k_y
   gama=gama2
   a=a_2
   e0=rho0+(a_2*rho0)+((k_y*(rho0**gama2))/(gama2-1))
  
  if(rho0>=rho2):
   k=k_z
   gama=gama3
   a=a_3
   e0=rho0+(a_3*rho0)+((k_z*(rho0**gama3))/(gama3-1))
  
  
  #e0=rho0+((k*(rho0**gama))/(gama-1))                                                                            # The value of  mass density/energy density for central point
   
  m0=(4*np.pi*(r0**3)*e0)/3                                                                                      # The mass of star at the central point
  
  rho[0]=rho0
  
  m[0]=m0
  
  e[0]=e0
  
  # Here I define piecewise polytropic equation of state.
  
  def Eos(rho):
   """ 
   Input: 
       rho: rest-mass density
   Output:
       P: pressure
   """
   if(rho<rho_lower_2):
     k=k1_lower
     gama=gama1_lower
     p=k1_lower*(rho**gama1_lower)
     return p
     
   if(rho_lower_2<=rho<rho_lower_3):
     k=k2_lower
     gama=gama2_lower
     p=k2_lower*(rho**gama2_lower)
     return p
     
   if(rho_lower_3<=rho<rho_lower_4):
     k=k3_lower
     gama=gama3_lower
     p=k3_lower*(rho**gama3_lower)
     return p
     
   if(rho_lower_4<=rho<rho_lower_5):
     k=k4_lower
     gama=gama4_lower
     p=k4_lower*(rho**gama4_lower)
     return p
   
   if(rho_lower_5<=rho<rho1):
     k=k_x
     gama=gama1
     p=k_x*(rho**gama1)
     return p
   
   if(rho1<=rho<rho2):
     k=k_y
     gama=gama2
     p=k_y*(rho**gama2)
     return p
   
   if(rho>=rho2):
     k=k_z
     gama=gama3
     p=k_z*(rho**gama3)
     return(p)
    
   
     
  p[0]=Eos(rho[0])
  
  
  # Here I start to solve TOV equation.
  
  # Using function f I define the derivative of rest mass density w.r.t the radial point
  
  
  
  def f(r,rho,m):
   
   T=rho+(a*rho)+((k*(rho**gama))/(gama-1))+(k*(rho**gama))
   U=m+(4*np.pi*(r**3)*(k*(rho**gama)))
   V=k*gama*(rho**(gama-1))
   W=(2*m*r)-(r**2)
   f=(T*U)/(V*W)
   return(f)
   
   
   
  # Using the function g I define the derivative of mass w.r.t the radial point
   
  
  
  def g(r,rho,m):
   
   L=rho+(a*rho)+((k*(rho**gama))/(gama-1))
   g=((L)*4*np.pi*(r**2))
   return(g)
   
   
   
 # To solve the rest mass density and mass at each radial point I use 4th order runge kutta method by this for loop.
 
 
   
  for i in range(1,1000001):
 
  
   if(rho[i-1]>=0.0):#((1.e+3)*G/(C**2))):#(((10**3)*G)/(C**2))):
   
    if(rho[i-1]<=rho_lower_2):
     k=k1_lower
     gama=gama1_lower
     a=a1_lower
     
    if(rho_lower_2<rho[i-1]<=rho_lower_3):
     k=k2_lower
     gama=gama2_lower
     a=a2_lower
     
    if(rho_lower_3<rho[i-1]<=rho_lower_4):
     k=k3_lower
     gama=gama3_lower
     a=a3_lower
     
    if(rho_lower_4<rho[i-1]<=rho_lower_5):
     k=k4_lower
     gama=gama4_lower
     a=a4_lower
   
    if(rho_lower_5<rho[i-1]<=rho1):
     k=k_x
     gama=gama1
     a=a_1
   
    if(rho1<rho[i-1]<=rho2):
     k=k_y
     gama=gama2
     a=a_2
   
    if(rho[i-1]>rho2):
     k=k_z
     gama=gama3
     a=a_3
    k1=h*f(r[i-1],rho[i-1],m[i-1])
    l1=h*g(r[i-1],rho[i-1],m[i-1])
    k2=h*f((r[i-1]+(h/2)),(rho[i-1]+(k1/2)),(m[i-1]+(l1/2)))
    l2=h*g((r[i-1]+(h/2)),(rho[i-1]+(k1/2)),(m[i-1]+(l1/2)))
    k3=h*f((r[i-1]+(h/2)),(rho[i-1]+(k2/2)),(m[i-1]+(l2/2)))
    l3=h*g((r[i-1]+(h/2)),(rho[i-1]+(k2/2)),(m[i-1]+(l2/2)))
    k4=h*f((r[i-1]+(h)),(rho[i-1]+(k3)),(m[i-1]+(l3)))
    l4=h*g((r[i-1]+(h)),(rho[i-1]+(k3)),(m[i-1]+(l3)))
    rho[i]=rho[i-1]+(((k1+(2*k2)+(2*k3)+k4)/6))
    m[i]=m[i-1]+(((l1+(2*l2)+(2*l3)+l4)/6))
    p[i]=Eos(rho[i])
    e[i]=(rho[i]*(1+a))+(p[i]/(gama-1))
    
   
   else:
    M=(m[i-2]*(C**2))/((G)*(M_sun))
    R=(r[i-2]*1.e-05)
    RHO=(rho[0]*(C**2))/G
    E=(e[0]*(C**4))/G
    P=(p[0]*(C**4))/G
    break
    
    
    
  return M, R, RHO, E, P
  
  

G=6.67e-08                                                                                  # Gravitational wave constant value  

C=3.e10                                                                                     # velocity of light

rho01=(((10**14.6))*G)/(C**2)                                                                    # the starting value of central rest mass density

rho02=((10**(15.6))*G)/(C**2)                                                                    # the end value of central rest mass density

rho0=np.logspace(np.log10(rho01),np.log10(rho02),100)

M1=[]                                                                                       # this array I store the maximum mass in each central density

R1=[]                                                                                       # this array I store the radius at maximum mass in each central density

RHO1=[]                                                                                     # this array I store the central rest mass density

E1=[]                                                                                       # this array I store the central energy density

P1=[]




# Here I plot mass - radius curve for different central densities from the above T0V solution.

                  

for i in range(100):
 M,R,RHO,E,P=Tov_solver(rho0[i])
 #print(Tov_solver(rho0[i]))
 M1.append(M)
 R1.append(R)
 RHO1.append(RHO)
 E1.append(E)
 P1.append(P)




plt.scatter(R1,M1)
plt.xlabel("R [km]")
plt.ylabel("M [M_sun]")
plt.title("AP3 EOS")
plt.show()


# Here I calculate the maximum mass of star, corresponding radius, corresponding rest mass density & corresponding mass density/energy density



M_max=max(M1)
maxpos1=M1.index(M_max)
R_max=R1[maxpos1]
RHO_max=RHO1[maxpos1]
E_max=E1[maxpos1]
print(" The values from direct Tov solution:")
print("maximum mass of star= ",M_max)
print("Radius of maximum mass star= ",R_max)
print("central rest mass density of maximum mass star= ",RHO_max)
print("central energy density of maximum mass star= ",E_max)








# Here I use Interpolation formula to get mass, radius, central density of any arbitary star



f=CubicSpline(E1,M1)                                                            # Where f is the interpolation function between energy density(x axis) & mass(y axis).
E2=np.logspace(np.log10(E1[0]),np.log10(E1[-1]),10000)                            # In E2 store the array of E1 for interpolation
M2=f(E2)
M3=[]                                                                        # In M2 store the masses for each E2 value
for i in range(10000):
  M2[i]=f(E2[i])
  M3.append(M2[i])
M1_max=max(M3)
maxpos=M3.index(M1_max)
E1_max=E2[maxpos]   
print("The values from Tov solution after interpolation:")
print("The maximum mass of the star from the interpolation is: ",M1_max)
print("The central energy density of the corresponding stable mass: " ,E1_max)
plt.plot(np.log10(E2),M3)
plt.scatter(np.log10(E1),M1)
plt.xlabel("log(E) [erg/cm^3]")
plt.ylabel("M [M_sun]")
plt.title("AP3 EOS ")
plt.show()







g=CubicSpline(E1,R1)                                                            # Where g is the interpolation function between energy density(x axis) & radius(y axis).
E2=np.logspace(np.log10(E1[0]),np.log10(E1[-1]),10000)                            # In E2 store the array of E1 for interpolation
R2=g(E2)                                                                        # In R2 store the radius of star for each E2 value corresponding to the M2 value.
R1_max=g(E1_max)
print("The radius of the stable maximum mass star: ",R1_max)
plt.plot(np.log10(E2),R2)
plt.scatter(np.log10(E1),R1)
plt.xlabel("log(E) [erg/cm^3]")
plt.ylabel("R [km]")
plt.title("AP3 EOS ")
plt.show()







U=CubicSpline(E1,RHO1)                                                            # Where U is the interpolation function between energy density(x axis) & rest mass density(y axis).
E2=np.logspace(np.log10(E1[0]),np.log10(E1[-1]),10000)                              # In E2 store the array of E1 for interpolation
RHO2=U(E2)

plt.plot(np.log10(E2),np.log10(RHO2))
plt.scatter(np.log10(E1),np.log10(RHO1))
plt.xlabel("log(E) [erg/cm^3]")
plt.ylabel("log(RHO) [gm/cm^3]")
plt.title("AP3 EOS")
plt.show()






t=CubicSpline(E1,P1)
E2=np.logspace(np.log10(E1[0]),np.log10(E1[-1]),10000)
P2=t(E2)
plt.plot(np.log10(E2),np.log10(P2))
plt.scatter(np.log10(E1),np.log10(P1))
plt.xlabel("log(E) [erg/cm^3]")
plt.ylabel("log(P) [dyne/cm^2]")
plt.title("AP3 EOS ")
plt.show()





W=CubicSpline(RHO1,P1)
RHO2=np.logspace(np.log10(RHO1[0]),np.log10(RHO1[-1]),10000)
P2=W(RHO2)
plt.plot(np.log10(RHO2),np.log10(P2))
plt.scatter(np.log10(RHO1),np.log10(P1))
plt.xlabel("log(RHO) [gm/cm^3]")
plt.ylabel("log(P) [dyne/cm^2]")
plt.title("AP3 EOS")
plt.show()



M4=[]                                   # M3 is array where I store only the stable value of M1 in M1-R1 curve i.e. starting to peak value & neglect the decay part because it is unstable part        
 
E4=[]

R4=[]                                  #  R3 is the array where I store the radius value of masses M1 in same way.
                                                                                
for i in range(maxpos1+1):
   M4.append(M1[i])



for i in range(maxpos1+1):
   E4.append(E1[i])
   
for i in range(maxpos1+1):
   R4.append(R1[i])




A=CubicSpline(M4,E4)                                                            # Where A is the interpolation function between stable part of mass(x axis) & radius(y axis) of star.
M5=np.linspace(M4[0],M4[-1],10000)                                                # In M4 store the array of M3 for interpolation
E5=A(M5)                                                                        # In R4 store the radius of star for each M4 value.
E6=A(1.4)
F=CubicSpline(E4,R4)                                                            # Where g is the interpolation function between energy density(x axis) & radius(y axis).
E7=np.logspace(np.log10(E4[0]),np.log10(E4[-1]),10000)                            # In E2 store the array of E1 for interpolation
R7=F(E7)
R8=F(E6)
print("Radius of star-----",R8)                                                                        # In R2 store the radius of star for each E2 value corresponding to the M2 value.
                           

 
 
 
 
 
 
 
