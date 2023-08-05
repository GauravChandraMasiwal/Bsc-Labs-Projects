
# PARACHUTE PROBLEM SOLVING

import csv
import numpy as np

from team14_program_E import euler
from team14_program_R import rk4
def func_vector(indep,dep): #dep is an array of dependent variables [y,v]   
    if n == 0:         #SATGE 1: FREE FALL
        Cd_A_net = Cd_flat * b0
    
    elif n == 1:        #STAGE 2: PULLING OF RIPCORD
        Cd_A_net = Cd_flat * b0 + (Cd_cyl * b1 * l * (indep - T0)/(T1 - T0))
        
    elif n == 2:        #STAGE 3: OPENING OF PARACHUTE UPTO a1 AREA
        Cd_A_net = Cd_cyl * b1 * h + (Cd_hem * alpha0 * np.exp(beta0 * (indep - T1)/(T2 - T1)))
        
    elif n == 3:        #STAGE 4 : OVER INFLATION OF PARACHUTE
        Cd_A_net = Cd_cyl * b1 * h + Cd_hem * alpha1 * (1 + beta1 *np.round( np.sin(np.pi * (indep - T2)/(T3 - T2)),10))

    else :              #STAGE 5 : REACHED STEADY AREA a1
        Cd_A_net = Cd_cyl * b1 * h + Cd_hem * a1
        
    k = 1/2 * p *Cd_A_net       #Cd_A_net = Cd(body)*A(body) + Cd(equip)*A(equip) 
    kappa = k/k1
    dy_new_dT = dep[1] 

    dvelo_new_dT = kappa*dep[1]**2- 1

    f_vector = np.array([dy_new_dT,dvelo_new_dT])
    return f_vector
    
def cal_k(t_arr,k_list):        #k_list IS THE LIST IN WHICH K VALUES WILL BE APPENDED FOR TIME VALUES OF t_arr           
    #TO CALCULATE THE DRAG CONSTANT FOR DIFFERENT STAGES
    for t in t_arr:
        if t <= T0:
            Cd_A_net = Cd_flat * b0
    
        elif t > T0 and t <= T1:
            Cd_A_net = Cd_flat * b0 + (Cd_cyl * b1 * l * (t - T0)/(T1 - T0))
        
        elif t >T1 and t <= T2:
            Cd_A_net = Cd_cyl * b1 * h + (Cd_hem * alpha0 * np.exp(beta0 * (t - T1)/(T2 - T1)))
        
        elif t > T2 and t <= T3:
            Cd_A_net = Cd_cyl * b1 * h + Cd_hem * alpha1 * (1 + beta1 *np.round( np.sin(np.pi * (t - T2)/(T3 - T2)),10))

        else :
            Cd_A_net = Cd_cyl * b1 * h + Cd_hem * a1
        
        k = 1/2 * p *Cd_A_net
        kappa = k / k1
        k_list.append(kappa)

if __name__ == "__main__" :
    g = 9.81        #ms^-2
    m = 97.2        #kg
    H = 5480.9136   #metres
    p = 1           #kg/m^3
    Cd_flat = 1.95
    Cd_cyl = 0.35
    Cd_hem = 1.33
    b0 = 0.5        #m^2
    b1 = 0.1        #m^2
    a1 = 43.8       #m^2
    h = 1.78        #metre
    l = 8.96        #metre
    
    t1 = np.sqrt(H/g)       #TIME FACTOR TO MAKE TIME VARIABLE DIMENTIONLESS
    v1 = np.sqrt(g*H)       #VELOCITY FACTOR TO MAKE VELOCITY VARIABLE DIMENTIONLESS
    k1 = m / H              #DRAG FACTOR TO MAKE DRAG CONSTANT DIMENTIONLESS
    T0 = 20/t1      #sec    # t0 IS THE DIMENTIONLESS TIME FOR FREE FALL
    T1 = 20.5/t1      #sec  
    T2 = 21.5/t1     #sec
    T3 = 23.2/t1       #sec
    T4 = 35/t1
    alpha0 = (Cd_flat * b0 + Cd_cyl * b1 * (l - h))/Cd_hem
    beta0  = np.log(a1/alpha0)
    beta1 = 0.15
    alpha1 = a1
    
    n = 0       #STAGE 1
    initial0 = [0,1,0]   # INITIAL CONDITIONS
    
    result1 = euler(initial0,T0,func_vector)
    result2 = rk4(initial0,T0,func_vector)
    
    n = 1       #STAGE 2
    
    initial_1_euler = [result1[0][-1],result1[1][0][-1],result1[1][1][-1]]      #INITIAL CONDITIONS FOR STAGE 2 WILL BE THE FINAL VALUES OBTAINED IN STAGE 1
    initial_1_rk4 = [result2[0][-1],result2[1][0][-1],result2[1][1][-1]]
    
    result3 = euler(initial_1_euler,T1,func_vector)
    result4 = rk4(initial_1_rk4,T1,func_vector)
    
    n = 2       #STAGE 3
    
    initial_2_euler = [result3[0][-1],result3[1][0][-1],result3[1][1][-1]]
    initial_2_rk4 = [result4[0][-1],result4[1][0][-1],result4[1][1][-1]]
    
    result5 = euler(initial_2_euler,T2,func_vector)
    result6 = rk4(initial_2_rk4,T2,func_vector)
    
    n = 3       #STAGE 4
    
    initial_3_euler = [result5[0][-1],result5[1][0][-1],result5[1][1][-1]]
    initial_3_rk4 = [result6[0][-1],result6[1][0][-1],result6[1][1][-1]]
    
    result7 = euler(initial_3_euler,T3,func_vector)
    result8 = rk4(initial_3_rk4,T3,func_vector)
    
    n = 4       #STAGE 5
    
    initial_4_euler = [result7[0][-1],result7[1][0][-1],result7[1][1][-1]]
    initial_4_rk4 = [result8[0][-1],result8[1][0][-1],result8[1][1][-1]]
    
    result9 = euler(initial_4_euler,T4,func_vector)
    result10 = rk4(initial_4_rk4,T4,func_vector)
    
    # CONCATENATION OF ALL STAGES ARRAYS
    
    res11 = np.concatenate([result1[0],result3[0],result5[0],result7[0],result9[0]])
    res12 = np.concatenate([result1[1][0],result3[1][0],result5[1][0],result7[1][0],result9[1][0]])
    res13 = np.concatenate([result1[1][1],result3[1][1],result5[1][1],result7[1][1],result9[1][1]])
    
    res21 = np.concatenate([result2[0],result4[0],result6[0],result8[0],result10[0]])
    res22 = np.concatenate([result2[1][0],result4[1][0],result6[1][0],result8[1][0],result10[1][0]])
    res23 = np.concatenate([result2[1][1],result4[1][1],result6[1][1],result8[1][1],result10[1][1]])
    
    res1 = np.array([res11,res12,res13])
    res2 = np.array([res21,res22,res23])
    
    drag_cons_1 = []
    drag_cons_2 = []
    
    # FOR DRAG CONSTANT ARRAYS 
    
    cal_k(res1[0], drag_cons_1)
    cal_k(res2[0], drag_cons_2)
    
    # WRITING THE VALUES OBTAINED INTO CSV FILES FOR EULER AND RK4
    
    header = [['(T) (EULER)','(Y) (EULER)','(V) (EULER)','(K) (EULER)'],['(T) (RK4)','(Y) (RK4)','(V) (RK4)','(K) (RK4)']]    
    
    with open('parachute_euler.csv','w',newline = "") as file :      # newline = " " , FOR SKIPPING THE EMPTY ROWS CREATED USING writerow
        writer = csv.writer(file)
        writer.writerow(header[0])
        
        for i in range(len(res1[0])):
        
            writer.writerow([res1[0][i],res1[1][i],res1[2][i],drag_cons_1[i]])   #TO WRITE THE VALUES OF DIFFRENT ARRAYS AS COLUMNS
    
    with open('parachute_rk4.csv','w',newline = "") as file :
        writer = csv.writer(file)
        writer.writerow(header[1])
        
        for i in range(len(res2[0])):
        
            writer.writerow([res2[0][i],res2[1][i],res2[2][i],drag_cons_2[i]])
        
    
       
