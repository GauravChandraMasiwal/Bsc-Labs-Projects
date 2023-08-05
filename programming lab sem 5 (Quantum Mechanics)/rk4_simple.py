# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 11:04:07 2022

@author: chand
"""

def Rk4(ini_cond,inde_f,func,N,e):  #RK4 METHOD for solving simulatenous differential equation
    h = (inde_f - ini_cond[0])/(N)       
    indep =[ini_cond[0]]
    dep1 = [ini_cond[1]]
    dep2 = [ini_cond[2]]
    for i in range(N):
        m1 = h*func(indep[i],[dep1[i],dep2[i]],e)
        
        m2 = h*func(indep[i] + (h/2),[dep1[i] + (m1[0]/2),dep2[i]+(m1[1]/2)],e)
        
        m3 = h*func(indep[i] + (h/2),[dep1[i] + (m2[0]/2),dep2[i]+(m2[1]/2)],e)
        
        m4 = h*func(indep[i] + h,[dep1[i]+m3[0],dep2[i]+m3[1]],e)
        
        mrk4 = (1/6)*(m1 + 2*m2 + 2*m3 + m4)
        
        indep.append(indep[i]+h)
        dep1.append(dep1[i]+mrk4[0])
        dep2.append(dep2[i]+mrk4[1])

    return [indep,[dep1,dep2]]