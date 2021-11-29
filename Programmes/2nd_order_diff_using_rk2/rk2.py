import pandas as pd
#%matplotlib inline
import matplotlib.pyplot as plt
from scipy import integrate
import numpy as np
def simpe_harmonic(X, t, cons):  #simple harmonic oscillator
    k,m = cons
    x, y = X
    dotx = y
    dot2x = -k*x/m
    return np.array([dotx, dot2x])
def damped_harmonic(X, t, cons):  #damped oscillator
    k,m,b = cons
    x, y = X
    dotx = y
    doty = -k*x-b*dotx
    return np.array([dotx, doty])
def simple_pendulum(X, t, cons):  #Simple Pendulum
    g,L = cons
    x, y = X
    dotx = y
    doty = -g*x/L
    return np.array([dotx, doty])
def RK2(func, X0, t,cons):
    dt = t[1] - t[0]
    nt = len(t)
    X  = np.zeros([nt, len(X0)])
    X[0] = X0
    for i in range(nt-1):
        k1 =dt* func(X[i], t[i], cons)
        k2 = dt*func(X[i] +  k1, t[i] + dt, cons)
        X[i+1] = X[i] + (k1 +  k2 )/2
    return X, t


if __name__ == "__main__":
    k = 1; m = 1; cons=(k,m); t_p = 2*np.pi*np.sqrt(m/k); tmin=0; tmax=7*t_p; Nt = 1000; x0 = [2,0] 
    t = np.linspace(tmin,tmax,Nt)
    Xrk2 = RK2(simpe_harmonic, x0, t, cons)
    x,y= Xrk2[0].T
    t=Xrk2[1]
    tdimless=t/t_p
    #Plotting
    
    plt.title('Simple Harmonic Oscillator', fontsize=20)    
    plt.plot(tdimless,x,color='black',label="Displacement");plt.plot(tdimless,y,color='brown',label="Velocity")
    plt.grid(); plt.legend()
    plt.show()
#Damped Harmonic Oscillator
    k = 1; m = 1; t_p = 2*np.pi*np.sqrt(m/k); tmin=0; tmax=30*t_p; Nt = 1000; x0 = [1,0]
    t = np.linspace(tmin,tmax,Nt)
    tdimless=t/t_p; dis=[];vel=[];time=[]
    ba=[0.15,2,5]
    for b in ba:
        cons=(m,k,b)
        Xrk2=RK2(damped_harmonic,x0,tdimless, cons)
        x,y= Xrk2[0].T
        tdimless=Xrk2[1]
        dis.append(x)
        vel.append(y)
    fig, axs = plt.subplots(3,figsize=(11,15))
    fig.suptitle('Damped Harmonic Oscillator', fontsize=20)
    axs[0].plot(tdimless,dis[0],label = "displacement")
    axs[0].set(xlabel="time ",title="Underdamped, b = 0.15")
    axs[0].plot(tdimless,vel[0],label = "velocity")
    axs[0].plot(tdimless,1/2*k*dis[0]**2)
    axs[0].grid(); axs[0].legend()
    axs[1].plot(tdimless,dis[1],label = "displacement")
    axs[1].set(xlabel="time ",title="Critically Damped, b = 2")
    axs[1].plot(tdimless,vel[1],label = "velocity")
    axs[1].plot(tdimless,1/2*k*dis[1]**2)
    axs[1].grid(); axs[1].legend()
    axs[2].plot(tdimless,dis[2],label = "displacement")
    axs[2].plot(tdimless,1/2*k*dis[2]**2)
    axs[2].set(xlabel="time ",title="Overdamped, b =5")
    axs[2].plot(tdimless,vel[2],label = "velocity")
    axs[2].grid(); axs[2].legend()
    plt.show()
#Simple pendulum
    g = 9.8; L = 1; cons=(g,L); t_p = 2*np.pi*np.sqrt(L/g); tmin=0; tmax=7*t_p; Nt = 1000; x0 = [2,0] 
    t = np.linspace(tmin,tmax,Nt)
    Xrk2 = RK2(simple_pendulum, x0, t, cons)
    x,y= Xrk2[0].T
    t=Xrk2[1]
    tdimless=t/t_p
#Plotting    
    plt.title('Simple Pendulum', fontsize=20)    
    plt.plot(tdimless,x,color='black',label="Angular Displacement");plt.plot(tdimless,y,color='brown',label="Angular Velocity")
    plt.grid(); plt.legend()
    plt.show()
