import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as integrate
def trapeziodal(func,a,b,n):
    y=[]
    h=(b-a)/n
    for i in range(n+1):    
        y.append(func(a+i*h))  #y at limit points
    trp=h*(func(a)+func(b))/2
    for j in range(1,len(y)-1):
        trp=trp+h*(y[j])
    return(trp)
def simpson(func,a,b,n):
    h=(b-a)/(2*n)
    simp=h*(func(a)+func(b))/3
    for i in range(1,2*n): 
        if(i%2==0):
            simp=simp+2*h*func(a+i*h)/3  
        elif(i%2==1):
            simp=simp+4*h*func(a+i*h)/3 
    return(simp)


def graph(func,a,b,n):
    N=np.arange(1,n+1,1)    #Number of Subintervals
    H2=(b-a)/(2*N)
    It, Is, Iq =[], [], []
    for i in N:
        z=trapeziodal(func,a,b,i)
        It.append(z)
        z1=simpson(func,a,b,i)
        Is.append(z1)
        analytic = integrate.quad(func, a, b)
        Iq.append(analytic)
    plt.scatter(H2,It,label="Trapezoidal",marker="*")
    plt.scatter(H2,Is,label="Simpson",marker=".")
    plt.scatter(H2,Is,label="Scipy Quad",marker=".")
    plt.yscale("log")
    plt.xscale("log")      
    plt.legend()
    plt.xlabel("h")
    plt.ylabel("I(h)")
    plt.title("Convergence Test")
    plt.grid(True)
    plt.show()


def Q2a():
    trp=trapeziodal(func,a,b,n)
    print("Integration Trapezoidal method",trp)

def Q2b():
    simp=simpson(func,a,b,n)
    print("Result by Simpson Method",float(simp))

def Q2c():
    
    analytic = integrate.quad(func, a, b)
    print("Analytic Solution = ", analytic[0])
    trap=trapeziodal(func,a,b,n)
    print("Solution by trapezoidal method: ",trap)
    err=abs(analytic[0]-trap)
    print("Truncation Error = ",float(err))

def Q2d():
    graph(func,a,b,n)


def Q3a():
    v=np.array([0.0, 0.5, 2.0, 4.05, 8.0, 12.5, 18.0, 24.5, 32.0, 40.5, 50.0])
    c=np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    h=(c[-1]-c[0])/(len(c)-1)
    pwr_trap=h*(v[0]+v[-1])/2
    for j in range(1,len(v)-1):
        pwr_trap=pwr_trap+h*(v[j])
    print("Power using trapezoidal method: ", float(pwr_trap), "J")
    pwr_simp=h*(v[0]+v[-1])/3
    for k in range(1,(len(v)-1)): 
        if(k%2==0):
            pwr_simp=pwr_simp+2*h*(v[k])/3  
        elif(k%2==1):
            pwr_simp=pwr_simp+4*h*(v[k])/3
    print("Power using simpson method : ", float(pwr_simp), "J")
    
def Q3b():
    h=(b-a)/n
    trp=trapeziodal(func,a,b,n)
    print("Integral is {:.8} using Trapezoidal Method".format(float(trp)))
    #Simpson
    simp=simpson(func,a,b,n)
    print("Integral is {:.8} using Simpson Method".format(float(simp)))
    graph(func,a,b,n)

if __name__ == "__main__":
    func=eval("lambda x:"+input("F: "))
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    Q2a()
    Q2b()
    Q2c()
    Q2d()
    Q3a()
    Q3b()