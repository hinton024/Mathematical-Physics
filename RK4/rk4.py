import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def eqn(X, t):
    x, y = X
    dotx = y + x - x**3
    doty = -x
    return np.array([dotx, doty])
def RK4(func, X0, t):
    dt = t[1] - t[0]
    nt = len(t)
    X  = np.zeros([nt, len(X0)])
    X[0] = X0
    for i in range(nt-1):
        k1 = func(X[i], t[i])
        k2 = func(X[i] + dt/2. * k1, t[i] + dt/2.)
        k3 = func(X[i] + dt/2. * k2, t[i] + dt/2.)
        k4 = func(X[i] + dt    * k3, t[i] + dt)
        X[i+1] = X[i] + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X

def graph(t,x,y,x1,y1,x2,y2,x3,y3,title):
    fig,axs=plt.subplots(2,2,figsize=(7,7))
    fig.suptitle(title, fontsize=30)
    ax11,ax12,ax21,ax22=axs[0][0],axs[0][1],axs[1][0],axs[1][1]
    
    ax11.scatter(t,x,color='black',label="x"),ax11.scatter(t,y,color='brown',label="y")
    ax11.set_title("x(0)=0 & y(0)=-1")      
    ax12.scatter(t,x1,color='black',label="x"),ax12.scatter(t,y1,color='brown',label="y")
    ax12.set_title("x(0)=0 & y(0)=-2")      
    ax21.scatter(t,x2,color='black',label="x"),ax21.scatter(t,y2,color='brown',label="y")
    ax21.set_title("x(0)=0 & y(0)=-3"),ax21.set_xlabel("Time")      
    ax22.scatter(t,x3,color='black',label="x"),ax22.scatter(t,y3,color='brown',label="y")
    ax22.set_title("x(0)=0 & y(0)=-4"),ax22.set_xlabel("Time")    

    ax11.legend(),ax11.grid(True),ax12.legend(),ax12.grid(True),ax21.legend(),ax21.grid(True),ax22.legend(),ax22.grid(True)
    plt.show()
def graph1(x,y,x1,y1,x2,y2,x3,y3,title):
    fig,axs=plt.subplots(2,2,figsize=(7,7))
    fig.suptitle(title, fontsize=30)
    ax11,ax12,ax21,ax22=axs[0][0],axs[0][1],axs[1][0],axs[1][1]   
    ax11.scatter(x,y,color='black')
    ax11.set_title("1st condition");ax11.set_ylabel("y")      
    ax12.scatter(x1,y1,color='black')
    ax12.set_title("2nd condition");      
    ax21.scatter(x2,y2,color='black')
    ax21.set_title("3rd condition");ax21.set_xlabel("x");ax21.set_ylabel("y")      
    ax22.scatter(x3,y3,color='black')
    ax22.set_title("4th condition");ax22.set_xlabel("x")   
    ax11.legend(),ax11.grid(True),ax12.legend(),ax12.grid(True),ax21.legend(),ax21.grid(True),ax22.legend(),ax22.grid(True)
    plt.show()

if __name__ == "__main__":
    x0 = [0,0,0,0];y0 = [-1,-2,-3,-4];Nt = 100;tmax = 15
    t = np.linspace(0.,tmax, Nt)
    X0 = [x0[0], y0[0]] 
    X1 = [x0[1], y0[1]]
    X2 = [x0[2], y0[2]]
    X3 = [x0[3], y0[3]]
    res = RK4(eqn, X0, t)
    res1 = RK4(eqn, X1, t)
    res2 = RK4(eqn, X2, t)
    res3 = RK4(eqn, X3, t)
    x, y = res.T; 
    x1, y1 = res1.T
    x2, y2 = res2.T
    x3, y3 = res3.T
    graph(t,x,y,x1,y1,x2,y2,x3,y3,"RK4 Method")
    graph1(x,y,x1,y1,x2,y2,x3,y3,"x vs y")
    data = {"Time":t,"x_rk4":x2,"y_rk4":y2}
    print(pd.DataFrame(data))