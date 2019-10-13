
import numpy as np
from numpy import genfromtxt
from numpy  import array

from scipy import signal
import matplotlib.pyplot as plt


t=np.arange(start = 0, stop = .2, step = 10e-7)
v1=2000000*(np.sin(2*np.pi*60*t))+100000*(np.sin(2*np.pi*1200*t))
t1=1000000*t #Time vector in micro-seconds



print(len(t))

Ts=t[100]-t[99] # Sample Time
Fs=1/Ts # sampling Frequency
FNy=Fs/2 # Nyquist Frequncy

print('Sampling Time (us):',round(float(1e6*Ts),4))
print('Sampling Frequency(MHz) :',round(float(Fs/1e6),2))
print('Nyquist Frequency(MHz) :',round(float(FNy/1e6),2))

def runningMean(x, N):# function to calculate moving average of a signal
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
         y[ctr] = np.sum(x[ctr:(ctr+N)])
    return y/N
  
v1m=runningMean(v1,1)

samples=len(v1m)


Ar=1.5e-3 #Cross section area of the coil
Ns=240 #Secondary winding turns 
Np=1 #primary winding turns
l=1.0076 #Mean flux path length

ks1=1154 # JA theory parameter of the core (k)
Beta=0.9 # shoulder adjustment of the BH curve proposed by Prof Annakkage
Ms=1.57e6 #Saturation magnetizationm
C=0.0198  # Domain flexing parameter of JA model
A1=7.09e-4#Inter domain coupling parameter of JA model
M0=4*3.14159*1e-7 #Permeability in vaccum
A=499 #JA model parameter A for the Langevin function
Beta2=1e-9



def simrun(Ip,sam,Rd1,Ts1,b2):
  H=np.zeros(sam)
  He=np.zeros(sam)
  V=np.zeros(sam)
  He=np.zeros(sam)
  B=np.zeros(sam)
  M=np.zeros(sam)
  MN=np.zeros(sam)
  Dir=np.zeros(sam)
  Mirr=np.ones(sam)
  k1=np.zeros(sam)
  V=np.zeros(sam)
  Const5=np.ones(sam)  

  
  
  for i in range(1,len(Ip)-1):
    
    n=.5*1e-6
    
    dIp=(Ip[i]-Ip[i-1])
    
    dH=(1/(l*Ns))*dIp
      
    H[i]=H[i-1]+dH+n*dH/Ts
  
    if dH>=0.0:
      Dir[i]=1
    else:
      Dir[i]=-1
      
      
    if H[i]*Dir[i]>0.0:
      k1u=ks1*(1+Beta*1*(M[i-1]/Ms)**2)
    else:
      k1u=ks1*(1-Beta*0*(M[i-1]/Ms)**2)
      
      
    if dH>01.0:
      k1m=k1u*(1+b2*1*(dH/Ts)**1)
    elif dH<-01.0:
      k1m=k1u*(1-b2*1*(dH/Ts)**1)
    else:
      k1m=k1u
      
      
    k1[i]=k1m
   
      
    He[i]=H[i]+A1*MN[i-1]   
    
    dHe=He[i]-He[i-1]
    
    coth=np.cosh(He[i]/A)/np.sinh(He[i]/A)
    x=A/He[i]
    
    
    MN[i]=Ms*(coth-x)
    
    
    coth=np.cosh((H[i]+A1*MN[i])/A)/np.sinh((H[i]+A1*MN[i])/A)
    
    x=A/(H[i]+A1*MN[i])
    
   
    
    if dH>0.0:
      Const5[i]=1+Rd1*dH/(Ts)   
    else:
      Const5[i]=1-Rd1*dH/(Ts)    
      
    
    dMN=dHe*(Ms*(((1-coth**2)*(1/A))+A*(1/(H[i]+A1*MN[i])**2)))
    
    dM=(dH*(1-C)*(MN[i]-M[i-1])+C*Dir[i]*k1[i]*dMN)/(Dir[i]*k1[i]-A1*(1-C)*(MN[i]-M[i-1]))/Const5[i]
    
    dMirr=(dM-C*dMN)/(1-C)
    
    Mirr[i]=Mirr[i-1]+ dMirr
    
    M[i]=(C*MN[i]+(1-C)*Mirr[i])
      
    B[i]=M0*(H[i]+M[i])
   
    V[i]=Ns*Ar*M0*(1/Ts1)*(B[i]-B[i-1])
    
  return B,H,Const5,k1,V


ks1=1154 # JA theory parameter of the core (k)
Beta=0.09 # shoulder adjustment of the BH curve proposed by Prof Annakkage
Ms=1.57e6 #Saturation magnetizationm
C=0.0198  # Domain flexing parameter of JA model
A1=7.09e-4#Inter domain coupling parameter of JA model
M0=4*3.14159*1e-7 #Permeability in vaccum
A=499 #JA model parameter A for the Langevin function
Rd=1e-12
Beta2=0.0

Beta=.96
B,H,C5,k,Vol=simrun(v1m,samples,Rd,Ts,Beta2)
plt.figure(1)
plt.plot(H[50000:len(H)-1],B[50000:len(B)-1])
plt.xlabel('Magnetic feild strength H (A/m)')
plt.ylabel('Flux density B (T) ')
plt.show()
