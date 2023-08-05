#University Roll No : 20068567021
#Name : GAURAV CHANDRA
#College Roll No : 2020PHY1122 

#PART B : COMPUTATION

import numpy as np
import matplotlib.pyplot as plt

#EULER AND RK4 CODE FOR A SYSTEM OF N NUMBERS OF SIMULTANEOUS EQUATIONS,SO IT CAN ALSO BE USED TO SOLVE FOR A SECOND ORDER DIFFERENTIAL EQUATION

def Euler(ini_cond,inde_f,func,h,values):
    time_vect = np.array([ini_cond[0]])  
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)  
    N = int((inde_f - ini_cond[0])/h)  
    for i in range(N):
        m_vect = h*func(time_vect[i],y_vect[:,i],values)  
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + m_vect
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)       
        y_vect = np.array(Y)
    return [time_vect,y_vect]

def RK4(ini_cond,inde_f,func,h,values):       
    time_vect = np.array([ini_cond[0]])
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    N = int((inde_f - ini_cond[0])/h)  
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i],values)
        m2_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m1_vect/2),values)
        m3_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m2_vect/2),values)
        m4_vect = h*func(time_vect[i] + h,y_vect[:,i]+m3_vect,values)
        mrk4 = (1/6)*(m1_vect + 2*m2_vect + 2*m3_vect + m4_vect)
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + mrk4
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]

def simple_harmonic(t,X,values):
    
    x,v = X
    k,m = values
    dx_dt = v
    dv_dt = -(k*x)/m
    return np.array([dx_dt,dv_dt])


#function for plotting graph between the independent and dependent variables
def graph_sketch(x1,y1,x2,y2,title,savefigs):
    plt.scatter(x1,y1[0],label = "x(euler)")
    plt.scatter(x1,y1[1],label = "v(euler)")
    plt.scatter(x2,y2[0],label = "x(rk4)")
    plt.scatter(x2,y2[1],label = "v(rk4)")
    plt.title(title[0])
    plt.xlabel("TIME/TIME-PERIOD")
    plt.ylabel("MAGNITUDE")
    plt.grid()
    plt.legend()
    plt.savefig(savefigs[0])
    plt.show()
    
    # phasor plot between dependent variables
    plt.scatter(y1[0],y1[1],label = "EULER")
    plt.scatter(y2[0],y2[1],label = "RK4")
    plt.title(title[1])
    plt.xlabel("DISPLACEMENT")
    plt.ylabel("VELOCITY")
    plt.grid()
    plt.legend()
    plt.savefig(savefigs[1])
    plt.show()
    
#this function contains solutions of question 2,3 and 4
def q2_q3_q4():
    initial_condition = [0,1,1]  #list contains initial time,displacement and velocity
    k = 441 #N/m
    m = 49 #gm
    
    final_domain = 2*np.pi*(np.sqrt(m/k))  #sec  #this is the time peroid
    result_euler = Euler(initial_condition,2* final_domain, simple_harmonic, final_domain/100,[k,m])
    result_rk4 = RK4(initial_condition, 2*final_domain, simple_harmonic, final_domain/100,[k,m])
    title = ["SIMPLE HARMONIC OSCILLATOR (x(0)=1 , v(0) = 1)","PHASE PLOT FOR X VS V (x(0)=1 , v(0) = 1)"]
    savefigure = ["2020phy1122_1_A.png","2020phy1122_1_B.png"]
    
    graph_sketch(result_euler[0], result_euler[1], result_rk4[0],result_rk4[1], title, savefigure)
    
    print(final_domain)
    initial_condition2 = [0,0,5]
    
    result_euler_2 = Euler(initial_condition2,2* final_domain, simple_harmonic, final_domain/100,[k,m])
    result_rk4_2 = RK4(initial_condition2,2* final_domain, simple_harmonic, final_domain/100,[k,m])
    title2 = ["SIMPLE HARMONIC OSCILLATOR (x(0)=0 , v(0) = 5)","PHASE PLOT FOR X VS V (x(0)=0, v(0) = 5)"]
    savefigure2 = ["2020phy1122_2_A.png","2020phy1122_2_B.png"]
    
    graph_sketch(result_euler_2[0], result_euler_2[1], result_rk4_2[0],result_rk4_2[1], title2, savefigure2)
    
q2_q3_q4()


