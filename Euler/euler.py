#!/usr/bin/env python
# coding: utf-8

# In[106]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import texttable as tt
tab = tt.Texttable()
def euler(f,a,b,n,yinit):
    h=(b-a)/(n)
    xs = a+np.arange(n)*h
    ys=np.zeros(n)
    y = yinit
    for j,x in enumerate(xs):
        ys[j]=y
        y+=h*f(x,y)
    return xs, ys
def rk2(f,a,b,n,yinit):
    h=(b-a)/(n)
    xs = a+np.arange(n)*h
    ys=np.zeros(n)
    y = yinit
    for j,x in enumerate(xs):
        ys[j]=y
        k0 = h*f(x,y)
        y+=h*f(x+h/2,y+k0/2)
    return xs, ys
def rk4(f,a,b,n,yinit):
    h=(b-a)/(n)
    xs = a+np.arange(n)*h
    ys=np.zeros(n)
    y = yinit
    for j,x in enumerate(xs):
        ys[j]=y
        k0 = h*f(x,y)
        k1 = h*f(x+h/2,y+k0/2)
        k2 = h*f(x+h/2,y+k1/2)
        k3 = h*f(x+h,y+k2)
        y+=(k0+2*k1+2*k2+k3)/6
    return xs, ys
def Analytic(yinit,a,b,h,tau):
    xs,ys=[],[]
    p=np.arange(a,b,h)
    for t in p:        
        D=yinit*np.exp(-1*t/tau)
        ys.append(D)
        xs.append(t)        
    return xs,ys
def graph(xs_1,ys_1,xs_2,ys_2,xs_3,ys_3,xs_4,ys_4,xs_5,ys_5,xs_6,ys_6,xs_7,ys_7,xs_8,ys_8,xs_9,ys_9,xs_10,ys_10,title):
    fig,axs=plt.subplots(3,2,figsize=(15,15))
    fig.suptitle(title, fontsize=30)
    ax11,ax12,ax21,ax22,ax31,ax32=axs[0][0],axs[0][1],axs[1][0],axs[1][1],axs[2][0],axs[2][1]
    ax11.plot(xs_1,ys_1,'^', color='green',label="euler"),ax11.plot(xs_4,ys_4,'-', color='black',label="rk2")
    ax11.plot(xs_7,ys_7,'>',color='black',label="rk4"),ax11.plot(xs_7,ys_7,'*',color='brown',label="Analytic")
    ax11.set_title("Analytic v/s Euler v/s rk2 v/s rk4"),ax11.set_ylabel("N"),ax11.set_xlabel("Time")      
    ax12.plot(xs_1,ys_1,'^', color='green',label="h ={}".format(xs_1[1]-xs_1[0])),ax12.plot(xs_2,ys_2,'-', color='black',label="h ={}".format(xs_2[1]-xs_2[0]))
    ax12.plot(xs_3,ys_3,'>',color='black',label="h ={}".format(xs_3[1]-xs_3[0])),ax12.plot(xs_7,ys_7,'*',color='brown',label="Analytic")
    ax12.set_title("Euler for different stepsize(h)"),ax12.set_ylabel("N"),ax12.set_xlabel("Time")         
    ax21.plot(xs_4,ys_4,'^', color='green',label="h ={}".format(xs_4[1]-xs_4[0])),ax21.plot(xs_5,ys_5,'-', color='black',label="h ={}".format(xs_5[1]-xs_5[0]))
    ax21.plot(xs_6,ys_6,'>',color='black',label="h ={}".format(xs_6[1]-xs_6[0])),ax21.plot(xs_7,ys_7,'*',color='brown',label="Analytic")
    ax21.set_title("rk2 for different stepsize(h)"),ax21.set_ylabel("N"),ax21.set_xlabel("Time")                             
    ax22.plot(xs_7,ys_7,'^', color='green',label="h ={}".format(xs_7[1]-xs_7[0])),ax22.plot(xs_8,ys_8,'-', color='black',label="h ={}".format(xs_8[1]-xs_8[0]))
    ax22.plot(xs_9,ys_9,'>',color='black',label="h ={}".format(xs_9[1]-xs_9[0])),ax22.plot(xs_10,ys_10,'*',color='brown',label="Analytic")
    ax22.set_title("rk4 for different stepsize(h)"),ax22.set_ylabel("N"),ax22.set_xlabel("Time")           
    ax31.plot(xs_1,(ys_10-ys_1)/ys_10,'.', color='green',label="euler"),ax31.plot(xs_1,(ys_10-ys_4)/ys_10,'.', color='black',label="rk2")
    ax31.plot(xs_1,(ys_10-ys_7)/ys_10,'.', color='red',label="rk4")
    ax31.set_title("Error Plot at h = 0.4"),ax31.set_ylabel("Absolute Error"),ax31.set_xlabel("Time")               
    ax32.plot(xs_1,(ys_1-ys_4)/ys_1,'.', color='green',label="euler-rk2"),ax32.plot(xs_1,(ys_7-ys_4)/ys_7,'.', color='black',label="rk4-rk2")
    ax32.plot(xs_1,(ys_1-ys_7)/ys_1,'.', color='red',label="euler-rk4")
    ax32.set_title("Comparative Error Plot at h = 0.4"),ax32.set_ylabel("Absolute Error"),ax32.set_xlabel("Time")                 
    ax11.legend(),ax11.grid(True),ax12.legend(),ax12.grid(True),ax21.legend(),ax21.grid(True),ax22.legend(),ax22.grid(True)
    ax31.legend(),ax31.grid(True),ax32.legend(),ax32.grid(True)
    plt.show()
def q3_a(a,yinit,t_half):
    b = 5*t_half 
    tau=t_half/np.log(2)
    h = t_half/10
    n = int((b-a)/h)
    decay = lambda x, y: -1*y/tau
    xs_1, ys_1 = euler(decay,a,b,n,yinit)
    xs_2, ys_2 = euler(decay,a,b,2*n,yinit)
    xs_3, ys_3 = euler(decay,a,b,4*n,yinit)    
    xs_4, ys_4 = rk2(decay,a,b,n,yinit)
    xs_5, ys_5 = rk2(decay,a,b,2*n,yinit)
    xs_6, ys_6 = rk2(decay,a,b,4*n,yinit)
    xs_7, ys_7 = rk4(decay,a,b,n,yinit)
    xs_8, ys_8 = rk4(decay,a,b,2*n,yinit)
    xs_9, ys_9 = rk4(decay,a,b,4*n,yinit)
    xs_10, ys_10 = Analytic(yinit,a,b,h,tau)
    print("Radioactive Decay", "h =", h)
    headings_1 = ["t" ,"Analytic","euler","rk2","rk4","Ab_error euler","Ab_error rk2","Ab_error rk4"]
    tab.header(headings_1)
    for row in zip(xs_1,ys_10,ys_1,ys_4,ys_7,(ys_10-ys_1)/ys_10,(ys_10-ys_4)/ys_10,(ys_10-ys_7)/ys_10):
        tab.add_row(row)
        tab.set_max_width(0)
        tab.set_precision(6)
    s = tab.draw()
    print(s)
    tab.reset()
   
    graph(xs_1,ys_1,xs_2,ys_2,xs_3,ys_3,xs_4,ys_4,xs_5,ys_5,xs_6,ys_6,xs_7,ys_7,xs_8,ys_8,xs_9,ys_9,xs_10,ys_10,"Radioactive Decay")



def q3_b(a,yinit,R,C):
    b = 5*R*C 
    tau=R*C
    h = tau/10
    n = int((b-a)/h)
    rc = lambda x, y: -1*y/tau    
    xs_1, ys_1 = euler(rc,a,b,n,yinit)
    xs_2, ys_2 = euler(rc,a,b,2*n,yinit)
    xs_3, ys_3 = euler(rc,a,b,4*n,yinit)
    
    xs_4, ys_4 = rk2(rc,a,b,n,yinit)
    xs_5, ys_5 = rk2(rc,a,b,2*n,yinit)
    xs_6, ys_6 = rk2(rc,a,b,4*n,yinit)
    
    xs_7, ys_7 = rk4(rc,a,b,n,yinit)
    xs_8, ys_8 = rk4(rc,a,b,2*n,yinit)
    xs_9, ys_9 = rk4(rc,a,b,4*n,yinit)
    xs_10, ys_10 = Analytic(yinit,a,b,h,tau)
    print("RC Circuit", "h =", h)
    headings_1 = ["t" ,"Analytic","euler","rk2","rk4","$\delta$","Ab_error rk2","Ab_error rk4"]
    tab.header(headings_1)
    for row in zip(xs_1,ys_10,ys_1,ys_4,ys_7,(ys_10-ys_1)/ys_10,(ys_10-ys_4)/ys_10,(ys_10-ys_7)/ys_10):
        tab.add_row(row)
        tab.set_max_width(0)
        tab.set_precision(6)
    s = tab.draw()
    print(s)
    tab.reset()
   
    graph(xs_1,ys_1,xs_2,ys_2,xs_3,ys_3,xs_4,ys_4,xs_5,ys_5,xs_6,ys_6,xs_7,ys_7,xs_8,ys_8,xs_9,ys_9,xs_10,ys_10,"RC Circuit")

def q3_c(a,yinit,eta,rad,m):
    tau=m/((np.pi)*6*rad*eta)
    b = 5*tau
    h = tau/10
    n = int((b-a)/h)
    stokes = lambda x, y: -1*y/tau
    xs_1, ys_1 = euler(stokes,a,b,n,yinit)
    xs_2, ys_2 = euler(stokes,a,b,2*n,yinit)
    xs_3, ys_3 = euler(stokes,a,b,4*n,yinit)    
    xs_4, ys_4 = rk2(stokes,a,b,n,yinit)
    xs_5, ys_5 = rk2(stokes,a,b,2*n,yinit)
    xs_6, ys_6 = rk2(stokes,a,b,4*n,yinit)    
    xs_7, ys_7 = rk4(stokes,a,b,n,yinit)
    xs_8, ys_8 = rk4(stokes,a,b,2*n,yinit)
    xs_9, ys_9 = rk4(stokes,a,b,4*n,yinit)
    xs_10, ys_10 = Analytic(yinit,a,b,h,tau)   
    xs_10.pop(int(xs_10[-1]))
    ys_10.pop(int(ys_10[-1]))
    print("Stokes Law", "h =", h)
    headings_1 = ["t" ,"Analytic","euler","rk2","rk4","Ab_error euler","Ab_error rk2","Ab_error rk4"]
    tab.header(headings_1)
    for row in zip(xs_1,ys_10,ys_1,ys_4,ys_7,(ys_10-ys_1)/ys_10,(ys_10-ys_4)/ys_10,(ys_10-ys_7)/ys_10):
        tab.add_row(row)
        tab.set_max_width(0)
        tab.set_precision(6)
    s = tab.draw()
    print(s)
    tab.reset()    
    graph(xs_1,ys_1,xs_2,ys_2,xs_3,ys_3,xs_4,ys_4,xs_5,ys_5,xs_6,ys_6,xs_7,ys_7,xs_8,ys_8,xs_9,ys_9,xs_10,ys_10,"Stokes Law")
if __name__ == "__main__":
    q3_a(0,20000,4)
    q3_b(0,10,1e3,1e-6)
    q3_c(0,10,10,0.2,200)
   
   


# In[103]:






# In[ ]:





# In[104]:



   


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# \delta

# In[ ]:





# In[ ]:





# In[ ]:




