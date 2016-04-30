# -*- coding: utf-8 -*-
"""
Created on Sun Jan 03 20:52:38 2016

@author: Shubham Deshmane
"""
import scipy
import matplotlib.pyplot as plt
import scipy.optimize as opti
mA=1#kg/s
mB=2#kg/s
TAin=400#k
TBin=300#k
U=100#w/m^2k
p=1#m
L=10.0#m
n=10#grid pts
TAguess=n*[TAin]
TBguess=n*[TBin]
Tguess=scipy.array(TAguess+TBguess)
def CpA(T):
    CpA=4000+10.0*T+0.01*T**2
    return CpA
    
    
def CpB(T):
    CpB=3000+5.0*T+0.02*T**2
    return CpB
def residuals(T,U,L,TAin,TBin,mA,mB):
    n=len(T)
    TA=T[:n/2]
    TB=T[n/2:]
    dL=L/((n-1))
    errorALHS=((U*(TAin-TB[0]))/(mA*CpA(TAin)))+((TA[1]-TAin)/dL)
    errorBLHS=((U*(TAin-TB[0]))/(mB*CpB(TB[0])))+((TB[1]-TB[0])/dL)
    errorARHS=((U*(TA[-1]-TBin))/(mA*CpA(TA[-1])))+((TA[-1]-TA[-2])/dL)
    errorBRHS=((U*(TA[-1]-TBin))/(mB*CpB(TBin)))+((TBin-TB[-2])/dL)
    errorA=scipy.zeros(n/2)
    errorB=scipy.zeros(n/2)
    errorA[0]=errorALHS;errorA[-1]=errorARHS
    errorB[0]=errorBLHS;errorB[-1]=errorBRHS
    ##central difference
    errorA[1:-1]=(U*(TA[1:-1]-TB[1:-1])/(mA*CpA(TA[1:-1])))+(((TA[2:])-(TA[1:-1]))/dL)
    errorB[1:-1]=(U*(TA[1:-1]-TB[1:-1])/(mB*CpB(TB[1:-1])))+(((TB[2:])-(TB[1:-1]))/dL)
    return scipy.concatenate((errorA,errorB))
    
n=len(Tguess) 
ans=opti.leastsq(residuals,Tguess,args=(U,L,TAin,TBin,mA,mB))
print (ans)
Tans=ans[0]
TAans=Tans[:n/2]
TAans[0]=TAin
TBans=Tans[n/2:]
TBans[-1]=TBin
print(TAans)
print(TBans)
b=scipy.linspace(0,100,10)
plt.plot(b,TAans,'r')
plt.show()
plt.plot(b,TBans,'b')
plt.show()

#Energy of A nd B
EA=4000*(TAin-TAans[9])+5*(TAin**2-TAans[9]**2)+(0.01/3)*(TAin**3-TAans[9]**3)
EB=6000*(TBans[0]-TBin)+5*(TBans[0]**2-TBin**2)+(0.04/3)*(TBans[0]**3-TBin**3)
print TAans[9]
print TAin
print EA
print TBans[0]
print TBin
print EB
ERROR=((EA-EB)/(EA))*100
print ERROR

