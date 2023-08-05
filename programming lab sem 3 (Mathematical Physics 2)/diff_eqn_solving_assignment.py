# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 11:48:54 2021

@author: chand
"""
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tab 

def rk2(ini_cond,inde_f,func,h,val):        
    time_vect = np.array([ini_cond[0]])
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    N = int((inde_f - ini_cond[0])/h)  
    s_no = [1,]
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i],val)
        m2_vect = h*func(time_vect[i] + h,y_vect[:,i]+m1_vect,val)
        mrk2 = (1/2)*(m1_vect + m2_vect)
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + mrk2
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)
        y_vect = np.array(Y)
        s_no.append(i+2)
    return [s_no,time_vect,y_vect]
def harmonic_osc(t,X,val):
    x,v = X
    k,m,b = val
    dx_dt = v
    dv_dt = -(b*v + k*x)/m
    return np.array([dx_dt,dv_dt])

def simple_pen(t,X,val):
    θ,v = X
    g,L = val
    dθ_dt = v
    dv_dt = -(g/L)*θ
    return np.array([dθ_dt,dv_dt])

def coupled_pen(t,X,val):
    k,m,w0 = val
    xa,va,xb,vb = X
    dxa_dt = va
    dva_dt = -(w0*w0)*xa + (k/m)*(xa - xb)
    dxb_dt = vb
    dvb_dt = -(w0*w0)*xb - (k/m)*(xa - xb)
    return np.array([dxa_dt,dva_dt,dxb_dt,dvb_dt])

def graph_sketch(indep,dep,title,axis,labels,savefigs,queno):  
    if queno == 1 : 
        fig,axs = plt.subplots(2, 2)
        fig.tight_layout()
        ax = [[0,0],[0,1],[1,0],[1,1]]
        for i in range(len(dep)):
            a = ax[i]
            for j in range(len(dep[i])):
                axs[a[0],a[1]].scatter(indep[i],dep[i][j],label=labels[j],s=6) 
            axs[a[0],a[1]].set_xlabel(axis[0])
            axs[a[0],a[1]].set_ylabel(axis[1])
            axs[a[0],a[1]].set_title(title[i])
            axs[a[0],a[1]].grid()
            axs[a[0],a[1]].legend(loc = "upper right")
        plt.savefig(savefigs[0])
        
        plt.show()
      
    else : 
        fig,axs = plt.subplots(2)
        fig.tight_layout()
        ax = [0,1]
        for i in range(len(dep)):
            for j in range(len(dep[i])):
                axs[ax[i]].scatter(indep,dep[i][j],label=labels[j],s=6)
            axs[ax[i]].set(xlabel=axis[0],ylabel=axis[1][i],title = title[i])
            axs[ax[i]].grid()
            axs[ax[i]].legend(loc = "upper right")
        plt.savefig(savefigs[0])
        plt.show()
   
def q1():
    cond = [0,2,1]
    k = 4;m = 0.5;a = np.sqrt(4*k*m)
    b = [0,a-2,a,a+2]
    Tlist=[2*np.pi*(np.sqrt(m/k)),2,2,2]
    list1 = [0,0,0,0]
    for i in range(len(b)):
        list1[i] = rk2(cond, 4*Tlist[i], harmonic_osc, Tlist[i]/50, [k,m,b[i]])
    indep = np.array([list1[0][1]/Tlist[0],list1[1][1]/Tlist[1],list1[2][1]/Tlist[2],list1[3][1]/Tlist[3]])
    dep = [list1[0][2],list1[1][2],list1[2][2],list1[3][2]]
    print("FOR SIMPLE HARMONIC OSCILLATOR (b = 0)","\n")
    print(tab({'S.NO':list1[0][0],'TIME':list1[0][1],'DISPLACEMENT':list1[0][2][0],'VELOCITY':list1[0][2][1]},headers = 'keys'))
    print("\n","FOR UNDERDAMPED HARMONIC OSCILLATOR (b ="+str(round(b[1],1))+" )","\n")
    print(tab({'S.NO':list1[1][0],'TIME':list1[1][1],'DISPLACEMENT':list1[1][2][0],'VELOCITY':list1[1][2][1]},headers = 'keys'))
    print("\n","FOR CRITICALLY DAMPED HARMONIC OSCILLATOR (b ="+str(round(b[2],1))+" )","\n")
    print(tab({'S.NO':list1[2][0],'TIME':list1[2][1],'DISPLACEMENT':list1[2][2][0],'VELOCITY':list1[2][2][1]},headers = 'keys'))
    print("\n","FOR OVERDAMPED HARMONIC OSCILLATOR (b ="+str(round(b[3],1))+" )","\n")
    print(tab({'S.NO':list1[3][0],'TIME':list1[3][1],'DISPLACEMENT':list1[3][2][0],'VELOCITY':list1[3][2][1]},headers = 'keys'))
    title = ["SIMPLE HARMONIC (b = 0)","UNDERDAMPED OSCILLATOR","CRITICALLY DAMPED OSCILLATOR","OVERDAMPED OSCILLATOR"]
    axis = ["TIME/TIME PERIOD","MAGNITUDE"]
    labels = ["Displacement (x)","Velocity (v)"]
    savefigs = ["fig_1.png"]
    graph_sketch(indep, dep, title, axis, labels, savefigs, 1)
def q2():
    cond2 = [0,np.radians(15),1]
    L,g = 0.5,9.8 ; val2 = [g,L]
    T2=2*np.pi*(np.sqrt(L/g))
    sp = rk2(cond2, 5*T2,simple_pen, T2/50, val2)
    print("FOR SIMPLE PENDULUM","\n")
    print(tab({'S.NO':sp[0],'TIME':sp[1],'ANGULAR DISPLACEMENT':sp[2][0],'ANGULAR VELOCITY':sp[2][1]},headers = 'keys'))
    plt.scatter(sp[1]/T2,sp[2][0],marker = "*",label = 'ANGULAR DISPLACEMENT')
    plt.scatter(sp[1]/T2,sp[2][1],marker = ".",label = 'ANGULAR VELOCITY')
    plt.xlabel("TIME/TIME PERIOD")
    plt.ylabel("MAGNITUDE")
    plt.title("SIMPLE PENDULUM")
    plt.grid();plt.legend(loc = "upper right")
    plt.savefig("fig_2.png")
    plt.show()
def q3():
    cond3 = [0,np.radians(15),1,np.radians(15),1]
    k2=9.8;m2=10;w0=np.sqrt(k2/m2);T3 = 2*np.pi/w0
    val3 = [k2,m2,w0]
    cp = rk2(cond3,5*T3,coupled_pen, T3/50, val3)
    print("FOR COUPLED PENDULUM ","\n")
    print(tab({'S.NO':cp[0],'TIME':cp[1],'DISPLACEMENT (a)':cp[2][0],'VELOCITY (a)':cp[2][1],'DISPLACEMENT (b)':cp[2][2],'VELOCITY (b)':cp[2][3]},headers = 'keys'))
    dep1 = [cp[2][0],cp[2][1]]
    dep2 = [cp[2][2],cp[2][3]]
    axis = ["TIME/TIME PERIOD","MAGNITUDE"]
    labels = ["Displacement (x)","Velocity (v)"]
    title = ["COUPLED PENDULUM (a)","COUPLED PENDULUM (b)"]
    savefigs = ["fig_3.png"]
    graph_sketch(cp[1]/T3, np.array([dep1,dep2]), title, axis, labels, savefigs, 3)

q1(),q2(),q3()