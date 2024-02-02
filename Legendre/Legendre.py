import matplotlib.pyplot as plt
import numpy as np
from scipy.special import eval_legendre,legendre
from scipy.integrate import quad
import math
n = float(input("enter the positive integer n : "))
x = float(input("enter the value of x : "))   
def gm(n):  
    if n == 1:
        return 1    
    elif n == 0.5:
        return np.sqrt(np.pi)   
    else :
        return (n-1)*gm(n-1)
def leg(n,x,m=0,p=0):    
    if (n %2) == 0 :
        m = int(n/2)
    else : 
        m = int((n-1)/2)
    for i in range(m+1):
        p+= (((-1)**i) * (gm(2*n-2*i+1))* (x**(n-2*i)))/((2**n) * gm(i+1) * gm(n-i+1) * gm(n-2*i+1))        
    return p

def leg_dif(n,x,m=0,p=0):    
    m = 0
    if (n %2) == 0 :
        m = int(n/2)
    else : 
        m = int((n-1)/2)
    p = 0
    for i in range(m+1):
        p+= ((n-2*i)*((-1)**i) * (gm(2*n-2*i+1))* (x**(n-2*i-1)))/((2**n) * gm(i+1) * gm(n-i+1) * gm(n-2*i+1))       
    return p

print("legendre polynomial: ",leg(n,x,m=0,p=0))
print("Legendre derivative: ",leg_dif(n,x,m=0,p=0))
print('inbuilt leg :',eval_legendre(n,x))  
print('inbuilt leg diff',np.polyval(legendre(n).deriv(),x))
xdata = np.linspace(-0.999999,0.999999,100)
p0,p1,p2,p3,p0_diff,p1_diff,p2_diff,p3_diff,=[],[],[],[],[],[],[],[]


for i in xdata:
    p0.append(leg(0 ,i))
    p1.append(leg(1 ,i))
    p2.append(leg(2, i))
    p3.append(leg(3, i))
    p0_diff.append(leg_dif(0, i))
    p1_diff.append(leg_dif(1, i))
    p2_diff.append(leg_dif(2, i))
    p3_diff.append(leg_dif(3, i))        
def writeintofile(file1,xdata,p0,p1,p2):
    with open(file1,'w') as file :    
        for i in range(len(xdata)):
            file.write(str(xdata[i])+' , '+str(p0[i])+' , '+str(p1[i])+' , '+str(p2[i])+'\n')                   
writeintofile('leg00.dat', xdata, p0, p1, p2) 
writeintofile('leg01.dat', xdata, p1_diff, p2_diff, p3_diff)    
#plotting graph
plt.title("PLOT OF xdata VS Pn ",c= 'b')                                   
plt.plot(xdata,p0)
plt.plot(xdata,p1)
plt.plot(xdata,p2)
plt.legend(['p0','p1','p2'],loc = 'best')
plt.xlabel('xdata')
plt.ylabel('Pn')
plt.grid()
plt.show()
plt.title("PLOT OF xdata VS DIFF Pn")
plt.plot(xdata,p1)
plt.plot(xdata,p0_diff)
plt.plot(xdata,p2_diff)
plt.legend(['p1','p0 diff','p2 diff'],loc = 'best')
plt.xlabel('xdata')
plt.ylabel('Pn diff')
plt.grid()
plt.show()


n,LHS,RHS=2,[],[]
for i in range(len(xdata)):
    LHS.append(n*p2[i])
    RHS.append(xdata[i]*p2_diff[i] - p1_diff[i])
print("LHS \n",LHS[0:10])
print("RHS \n", RHS[0:10])
if np.allclose(LHS,RHS):
        print("Relation1 satisfied")
else:
    print("Relation1 not satisfied")
writeintofile('leg02.dat', xdata, p2, p2_diff, p1_diff) 

# relation 2
n,LHS,RHS=2,[],[]
for i in range(len(xdata)):
    LHS.append((2*n+1)*xdata[i]*p2[i])
    RHS.append((n+1)*p3[i]+n*p1[i])
print("LHS \n",LHS[0:20])
print("RHS \n", RHS[0:20])
if np.allclose(LHS,RHS):
        print("Relation2 satisfied")
else:
    print("Relation2 not satisfied")
writeintofile('leg03.dat', xdata, p2, p1, p3)    
#Relation 3
n,LHS,RHS=3,[],[]
for i in range(len(xdata)):
    LHS.append(n*p3[i])
    RHS.append((2*n-1)*xdata[i]*p2[i] - (n-1)*p1[i])
print("LHS \n",LHS[0:20])
print("RHS \n", RHS[0:20])
if np.allclose(LHS,RHS):
        print("Relation3 satisfied")
else:
    print("Relation3 not satisfied")
writeintofile('leg04.dat', xdata, p3, p2, p1) 

#orthogonality
A=[] ; B =[]
for n in range(3):
    for m in range(3):
        if n == m:
           A.append(2/(2*n+1))
        else:
           A.append(0)
        f=legendre(n)*legendre(m)
        inte , err = quad(f, -1, 1)
        B.append(inte)
RHS = np.array(B).reshape(3,3)
LHS = np.array(A).reshape(3,3)
print("RHS = ",RHS)
print("LHS = ",LHS)
if np.allclose(LHS,RHS):
   print("Orthogonality verified")
else:
    print("Orthogonality not verified")
