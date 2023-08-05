#this is a program to solve simple pendulum, damped and coupled pendulum and then plotting the graphs using RK2 method

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def damped_oscillator(t,Var,cons,b): 
    k,m= cons
    x,v = Var
    dx_dt = v    
    dv_dt = -(b/m)*v -(k/m)*x                         
    return np.array([dx_dt,dv_dt])
def simple_pendulum(t,Var,cons,θ):
    g,L=cons
    θ,v = Var
    dθ_dt = v    
    dv_dt = -(g/L)*θ
    return np.array([dθ_dt,dv_dt])
def coupled_pendulum(t,dep,cons,wo):
    k,m=cons
    dxa_dt=dep[1]
    dva_dt=-(wo**2+k/m*dep[0])
    dxb_dt=dep[3]
    dvb_dt=-(wo**2+k/m*dep[2])
    return np.array([dxa_dt,dva_dt,dxb_dt,dvb_dt])
def RK_2(tn,IC,func,h,cons,b):        #IC= [0,0,-1] or [0,0] where t=0,x=0 and y =-1 & t=0,i=0, tn is the final value
    time_vect = np.array([IC[0]])   #[0,0,-1] where t=0
    a = []
    for i in range(1,len(IC)):
        a.append([IC[i]])
    y_vect = np.array(a)
    N = int((tn - IC[0])/h)  
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i],cons,b)   #slope 1: m1 = hf(tn, xn);l1 = hg(tn,xn,yn)
        m2_vect = h*func(time_vect[i] + h ,y_vect[:,i] + m1_vect,cons,b)  #slope 2 : m2 = hf(tn + h, xn + m1);l2 = hg(tn + (h/2), xn + l1) 
        t = time_vect[i] + h        #ti = a + ih
        t_vect = np.append(time_vect,t)
        time_vect = t_vect 
        y_next = y_vect[:,i] + ((m1_vect+m2_vect)/2)   #x(t + h) = x(t) + ((m1+m2)/2 , (l1+l2)/2)
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j]);Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]

def graph():
    fig,ax = plt.subplots(2) 
    ax[0].scatter(RK2_1[0]/T1,RK2_1[1][0],c='indianred',marker="*",label='Displacement')
    ax[0].set(xlabel = "Time/Time Period",ylabel ='Displacement',title="Displacement vs Time")
    ax[0].legend(loc='best');ax[0].grid()
    ax[1].scatter(RK2_1[0]/T1,RK2_1[1][1],marker="*",label = 'Velocity',c='lightseagreen')
    ax[1].set(xlabel = "Time/Time Period",ylabel ='Velocity',title="Velocity vs Time")
    ax[1].legend(loc='best');ax[1].grid()
    fig.suptitle("Simple Harmonic osicallator (b=0)")
    plt.show()
def graph1():
     fig,ax = plt.subplots(3,2) 
     ax[0,0].scatter(RK2_2[0],RK2_2[1][0],c='deeppink',marker="*",label='Displacement')
     ax[0,0].set(xlabel = "Time",ylabel ='Displacement',title="Underdamped Harmonic osicallator(b^2-4km<0)")
     ax[0,0].legend(loc='best');ax[0,0].grid()
     ax[0,1].scatter(RK2_2[0],RK2_2[1][1],marker="*",label = 'Velocity',c='green')
     ax[0,1].set(xlabel = "Time",ylabel ='Velocity',title="Underdamped Harmonic osicallator(b^2-4km<0)")
     ax[0,1].legend(loc='best');ax[0,1].grid()
     ax[1,0].scatter(RK2_3[0],RK2_3[1][0],c='violet',marker="*",label='Displacement')
     ax[1,0].set(xlabel = "Time",ylabel ='Displacement',title="Critically Damped Harmonic osicallator(b^2-4km=0)")
     ax[1,0].legend(loc='best');ax[1,0].grid()
     ax[1,1].scatter(RK2_3[0],RK2_3[1][1],marker="*",label = 'Velocity',c='crimson')
     ax[1,1].set(xlabel = "Time",ylabel ='Velocity',title="Critically Damped Harmonic osicallator(b^2-4km=0)")
     ax[1,1].legend(loc='best');ax[1,1].grid()
     ax[2,0].scatter(RK2_4[0],RK2_4[1][0],c='purple',marker="*",label='Displacement')
     ax[2,0].set(xlabel = "Time",ylabel ='Displacement',title="Over Damped Harmonic osicallator(b^2-4km>0)")
     ax[2,0].legend(loc='best');ax[2,0].grid()
     ax[2,1].scatter(RK2_4[0],RK2_4[1][1],marker="*",label = 'Velocity',c='chocolate')
     ax[2,1].set(xlabel = "Time",ylabel ='Velocity',title="Over Damped Harmonic osicallator(b^2-4km>0)")
     ax[2,1].legend(loc='best');ax[2,1].grid()
     fig.suptitle("Damped Harmonic osicallator")
     plt.show()
    
def graph2():
    fig,ax = plt.subplots(2) 
    ax[0].scatter(RK2[0]/T,RK2[1][0],c='orange',marker="*",label='Displacement')
    ax[0].set(xlabel = "Time/Time Period",ylabel ='Angular Displacement',title="Displacement vs Time")
    ax[0].legend(loc='best');ax[0].grid()
    ax[1].scatter(RK2[0]/T,RK2[1][1],marker="*",label = 'Velocity',c='blue')
    ax[1].set(xlabel = "Time/Time Period",ylabel ='Angular Velocity',title="Velocity vs Time")
    ax[1].legend(loc='best');ax[1].grid()
    fig.suptitle("Simple Pendulum")
    plt.show()
    
def graph3():
    fig,ax = plt.subplots(2,2) 
    ax[0,0].scatter(RK_A[0],RK_A[1][0],c='limegreen',marker="*",label='Displacement')
    ax[0,0].set(xlabel = "Time",ylabel ='Angular Displacement',title="Angular Displacement ($\dfrac {d^2xa}{dt}$)")
    ax[0,0].legend(loc='best');ax[0,0].grid()
    ax[0,1].scatter(RK_A[0],RK_A[1][1],marker="*",label = 'Velocity',c='magenta')
    ax[0,1].set(xlabel = "Time",ylabel ='Angular Velocity',title="Angular Velocity vs Time ($\dfrac {d^2xa}{dt}$)")
    ax[0,1].legend(loc='best');ax[0,1].grid()
    ax[1,0].scatter(RK_A[0],RK_A[1][2],c='dodgerblue',marker="*",label='Displacement')
    ax[1,0].set(xlabel = "Time",ylabel ='Angular Displacement',title="Angular Displacement vs Time ($\dfrac {d^2xb}{dt}$)")
    ax[1,0].legend(loc='best');ax[1,0].grid()
    ax[1,1].scatter(RK_A[0],RK_A[1][3],marker="*",label = 'Velocity',c='peru')
    ax[1,1].set(xlabel = "Time",ylabel ='Angular Velocity',title="Angular Velocity vs Time ($\dfrac {d^2xb}{dt}$)")
    ax[1,1].legend(loc='best');ax[1,1].grid()
    fig.suptitle("Coupled Pendulum")
    plt.show()


if __name__ == "__main__" :
    y1=[0,1,0];k=4;m=0.5;cons=k,m
    T1=(2*np.pi*np.sqrt(m/k))
    b1=[0,0.5,np.sqrt(4*k*m),9]
    RK2_1 = RK_2(20,y1,damped_oscillator,0.1,cons,b1[0])
    RK2_2 = RK_2(20,y1,damped_oscillator,0.1,cons,b1[1])
    RK2_3 = RK_2(20,y1,damped_oscillator,0.1,cons,b1[2])
    RK2_4 = RK_2(20,y1,damped_oscillator,0.1,cons,b1[3])
    T1=(2*np.pi*np.sqrt(m/k))
    g=9.8;L=2;cons1=g,L
    T=(2*np.pi*np.sqrt(L/g))
    RK2=RK_2(20,y1,simple_pendulum,0.1,cons1,np.pi/2)
    y2=[0,1,2,2,1];k=4;m=0.5;cons=k,m
    RK_A=RK_2(20,y2,coupled_pendulum,0.1,cons,2)
    print("Table for Simple Pendulum")
    data={"Time":RK2_1[0]/T1,"Displacement":RK2_1[1][0],"Velocity":RK2_1[1][1]}
    print(pd.DataFrame(data))
    print("Table for Damped Harmonic Oscillator(Under Damped)")
    data1={"Time":RK2_2[0],"Displacement":RK2_2[1][0],"Velocity":RK2_2[1][1]}
    print(pd.DataFrame(data1))
    print("Table for Simple Pendulum(Underdamped)")
    data2={"Time":RK2[0]/T,"Angular Displacement":RK2[1][0],"Angular Velocity":RK2[1][1]}
    print(pd.DataFrame(data2))
    print("Table for Coupled Pendulum(Under Damped)")
    data2={"Time":RK_A[0],"Angular Displacement(for xa)":RK_A[1][0],"Angular Velocity(for xa)":RK_A[1][1],"Angular Displacement(for xb)":RK_A[1][2],"Angular Velocity(for xb)":RK_A[1][3]}
    print(pd.DataFrame(data2))
    pd.set_option('display.max_rows', None)
    graph()
    graph1()
    graph2()
    graph3()