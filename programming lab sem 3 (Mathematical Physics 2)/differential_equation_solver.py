#Here i solved a set of differential equation using different computation methods and then compared the results.
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def func_vector(t,dep):   #where t is the independent variable and x&y are the dependent variables
    x,y = dep
    dx_dt = x + y -x**3        #dx/dt = x+y-x^3
    dy_dt = -x                             #dy/dt = -x
    fnct_vector = np.array([dx_dt,dy_dt])
    return fnct_vector
def funct_vector1(t,dep):
    I1,I2 = dep
    dx_dt = -4*I1 + 4*I2 +12;dy_dt = -1.6*I1+1.2*I2+4.8                             
    fnct_vector1 = np.array([dx_dt,dy_dt])
    return fnct_vector1
def graph_sketch1(X,Y1,X1,Y2,X2,Y11):
    fig,ax = plt.subplots()
    plt.plot(X,Y1,c='deeppink',marker = '^',label='EULER')
    plt.plot(X1,Y2,c='green',marker = '^',label = 'RK2')
    plt.plot(X2,Y11,c='orange',marker = '^',label='RK4')
    ax.set_xlabel('I1(t)');ax.set_ylabel('I2(t)')
    plt.title("Graph for I2(t) vs I1(t)");plt.legend(loc = 'best')
    plt.grid()
    plt.show()
def graph_sketch2(X,Y1,Y2,X1,Y11,Y12,X2,Y21,Y22):
    fig,ax = plt.subplots()
    plt.plot(X,Y1,c='deeppink',marker = '^',label='I1(t)(EULER)')
    plt.plot(X,Y2,c='green',marker = '^',label = 'I2(t)(EULER)')
    plt.plot(X1,Y11,c='orange',marker = '^',label='I1(t)(RK2)')
    plt.plot(X1,Y12,c='chocolate',marker = '^',label = 'I2(t)(RK2)')
    plt.plot(X2,Y21,marker = '^',label='I1(t)(RK4)')
    plt.plot(X2,Y22,marker = '^',label = 'I2(t)(RK4)')
    ax.set_xlabel('time');ax.set_ylabel('I(t)')
    plt.title("Graph for t vs I(t)");plt.legend(loc = 'best')
    plt.grid()
    plt.show()
def graph(X,Y,Z,X_1,Y_1,Z_1,X_2,Y_2,Z_2,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X11,Y11,Z11,X12,Y12,Z12,X13,Y13,Z13,X22,Y22,Z22,X21,Y21,Z21,X31,Y31,Z31): 
     fig,ax = plt.subplots(2,2)  
     fig.tight_layout()
     ax[0,0].scatter(X, Y,marker = "*", label="x (EULER)", s = 6, c='darkorange')
     ax[0,0].scatter(X, Z,marker = "*", label="y (EULER)", s = 6 , c='blue')
     ax[0,0].scatter(X_1, Y_1,marker = "s", label="x (RK2)", s = 6, c='purple')
     ax[0,0].scatter(X_1, Z_1,marker = "s", label="y (RK2)", s = 6 , c='crimson')
     ax[0,0].scatter(X_2, Y_2,marker = "^", label="x (RK4)", s = 6, c='deeppink')
     ax[0,0].scatter(X_2, Z_2,marker = "^", label="y (RK4)", s = 6 , c='green')
     ax[0,0].set(xlabel = "time",ylabel ='Magnitude',title="x(0) = 0 and y(0) = -1")
     ax[0,0].legend(loc = "upper right"); ax[0,0].grid()    
     ax[0,1].scatter(X1, Y1,marker = "*", label="x (EULER)", s = 6, c='limegreen')
     ax[0,1].scatter(X1, Z1,marker = "*", label="y (EULER)", s = 6 , c='magenta')
     ax[0,1].scatter(X2, Y2,marker = "s", label="x (RK2)", s = 6, c='chocolate')
     ax[0,1].scatter(X2, Z2,marker = "s", label="y (RK2)", s = 6 , c='crimson')
     ax[0,1].scatter(X3, Y3,marker = "s", label="x (RK2)", s = 6, c='chocolate')
     ax[0,1].scatter(X3, Z3,marker = "s", label="y (RK2)", s = 6 , c='crimson')
     ax[0,1].set(xlabel = "time",ylabel ='Magnitude',title="x(0) = 0 and y(0) = -2")
     ax[0,1].legend(loc = "upper right"); ax[0,1].grid()     
     ax[1,0].scatter(X11, Y11,marker = "*", label="x (EULER)", s = 6, c="chocolate")
     ax[1,0].scatter(X11, Z11,marker = "*", label="y (EULER)", s = 6 , c = "palevioletred")
     ax[1,0].scatter(X12, Y12,marker = "s", label="x (RK2)", s = 6, c="green")
     ax[1,0].scatter(X12, Z12,marker = "s", label="y (RK2)", s = 6 , c = "pink")
     ax[1,0].scatter(X13, Y13,marker = "^", label="x (RK4)", s = 6, c="red")
     ax[1,0].scatter(X13, Z13,marker = "^", label="y (RK4)", s = 6 , c = "blue")
     ax[1,0].set(xlabel = "time",ylabel ='Magnitude',title="x(0) = 0 and y(0) = -3")
     ax[1,0].legend(loc = "upper right");ax[1,0].grid()   
     ax[1,1].scatter(X22, Y22,marker = "*", label="x (EULER)", s = 6)
     ax[1,1].scatter(X22, Z22,marker = "*", label="Y (EULER)", s = 6)
     ax[1,1].scatter(X21, Y21,marker = "s", label="x (RK2)", s = 6)
     ax[1,1].scatter(X21, Z21,marker = "s", label="Y (RK2)", s = 6)
     ax[1,1].scatter(X31, Y31,marker = "^", label="x (RK4)", s = 6)
     ax[1,1].scatter(X31, Z31,marker = "^", label="Y (RK4)", s = 6)
     ax[1,1].set(xlabel = "time",ylabel ='Magnitude',title="x(0) = 0 and y(0) = -4")
     ax[1,1].legend(loc = "best");ax[1,1].grid()
     plt.show()
def graph1(X,Y,X1,Y1,X2,Y2,X11,Y11,X12,Y12,X13,Y13,X21,Y21,X22,Y22,X23,Y23,X14,Y14,X24,Y24,X34,Y34):
     fig1,ax1 = plt.subplots(2,2)
     ax1[0,0].scatter(X, Y,marker = "*",label="EULER",s = 6,c='darkgreen')
     ax1[0,0].scatter(X1, Y1,marker = "s",label = 'RK2',s = 6,c='red')
     ax1[0,0].scatter(X2, Y2,marker = "^",label = "RK4",s = 6,c='pink')
     ax1[0,0].set(xlabel = "x",ylabel = 'y',title="1st Condition")
     ax1[0,0].legend(loc = "upper right");ax1[0,0].grid()
     ax1[0,1].scatter(X11,Y11,marker = "*",label = "EULER",s = 6,c='aqua')
     ax1[0,1].scatter(X12,Y12,marker = "s",label = "RK2",s = 6,c="pink")
     ax1[0,1].scatter(X13,Y13,marker = "^",label = "RK4",s = 6,c="khaki")
     ax1[0,1].set(xlabel = "x",ylabel = 'y',title="2nd Condition")
     ax1[0,1].legend(loc = "upper right");ax1[0,1].grid()   
     ax1[1,0].scatter(X21,Y21,marker = "*",label = "EULER",s = 6,c='darkblue')
     ax1[1,0].scatter(X22,Y22,marker = "s",label = "RK2",s = 6,c='deeppink')
     ax1[1,0].scatter(X23,Y23,marker = "s",label = "RK2",s = 6,c='green')
     ax1[1,0].set(xlabel = "x",ylabel = 'y',title="3rd Condition")
     ax1[1,0].legend(loc = "upper right");ax1[1,0].grid()    
     ax1[1,1].scatter(X14,Y14,marker = "*",label = "EULER",s = 6,c='deeppink')
     ax1[1,1].scatter(X24,Y24,marker = "s",label = "RK2",s = 6,c='orange')
     ax1[1,1].scatter(X34,Y34,marker = "^",label = "RK4",s = 6,c='blue')
     ax1[1,1].set(xlabel = "x",ylabel = 'y',title="4th Condition")
     ax1[1,1].legend(loc = "upper right");ax1[1,1].grid()   
     fig1.suptitle("y vs x");fig1.tight_layout()
     plt.show()
def euler(tn,IC,func,h):  #IC= [0,0,-1] or [0,0] where t=0,x=0 and y =-1 & t=0,i=0, tn is the final value
    time_vect = np.array([IC[0]])   #[0,0,-1] where t=0
    a = []
    for i in range(1,len(IC)):
        a.append([IC[i]])     #appended all the values of initial conditions
    y_vect = np.array(a)
    N = int((tn - IC[0])/h)     #N are the nodal points
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i])  #slope:m = hf(tn,xn,yn); l = hg(tn,xn,yn)
        t = time_vect[i] + h          #ti = a + ih
        t_vect = np.append(time_vect,t)
        time_vect = t_vect 
        y_next = y_vect[:,i] + m1_vect      #x(t + h) = x(t) + (m,l)
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j]);Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]
def RK_2(tn,IC,func,h):        #IC= [0,0,-1] or [0,0] where t=0,x=0 and y =-1 & t=0,i=0, tn is the final value
    time_vect = np.array([IC[0]])   #[0,0,-1] where t=0
    a = []
    for i in range(1,len(IC)):
        a.append([IC[i]])
    y_vect = np.array(a)
    N = int((tn - IC[0])/h)  
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i])   #slope 1: m1 = hf(tn, xn);l1 = hg(tn,xn,yn)
        m2_vect = h*func(time_vect[i] + h ,y_vect[:,i] + m1_vect)  #slope 2 : m2 = hf(tn + h, xn + m1);l2 = hg(tn + (h/2), xn + l1) 
        t = time_vect[i] + h        #ti = a + ih
        t_vect = np.append(time_vect,t)
        time_vect = t_vect 
        y_next = y_vect[:,i] + ((m1_vect+m2_vect)/2)   #x(t + h) = x(t) + ((m1+m2)/2 , (l1+l2)/2)
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j]);Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]   
def RK_4(tn,IC,func,h):        #IC= [0,0,-1] where t=0,x=0 and y =-1
    time_vect = np.array([IC[0]])  #[0,0,-1] where t=0
    a = []
    for i in range(1,len(IC)):
        a.append([IC[i]])
    y_vect = np.array(a)
    N = int((tn - IC[0])/h)   #N is the nodal points
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i])     #slope 1: m1 = hf(tn,xn,yn);l1 = hg(tn,xn,yn)
        m2_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m1_vect/2))   #slope 2: m2 = hf(tn + (h/2), xn + m1);l2 = hg(tn + (h/2), xn + l1) 
        m3_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m2_vect/2))   #slope 3: m3 = hf(tn + (h/2), xn + m2);l3 = gf(tn + (h/2), xn + l2)
        m4_vect = h*func(time_vect[i] + h,y_vect[:,i]+m3_vect)  #slope 4: m4 = hf(tn + h, xn + m3);l4 = gf(tn + h, xn + l3)
        mrk4 = (1/6)*(m1_vect + 2*m2_vect + 2*m3_vect + m4_vect) #mrk4 = (m1+2m2+2m3+m4)/6
        t = time_vect[i] + h     #ti = a + ih
        t_vect = np.append(time_vect,t)
        time_vect = t_vect 
        y_next = y_vect[:,i] + mrk4
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j]);Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]  
if __name__ == "__main__" :   
    y1 = [0,0,-1];y11 = [0,0,0];y2 = [0,0,-2];y3 = [0,0,-3];y4 = [0,0,-4]
    euler_1= euler(15,y1,func_vector,0.1);euler_2 = euler(15,y2,func_vector,0.1)
    euler_3 = euler(15,y3,func_vector,0.1);euler_4 = euler(15,y4,func_vector,0.1)    
    RK2_1 = RK_2(15,y1,func_vector,0.1);RK2_2 = RK_2(15,y2,func_vector,0.1)
    RK2_3 = RK_2(15,y3,func_vector,0.1);RK2_4 = RK_2(15,y4,func_vector,0.1)    
    RK4_1 = RK_4(15,y1,func_vector,0.1);RK4_2 = RK_4(15,y2,func_vector,0.1)
    RK4_3 = RK_4(15,y3,func_vector,0.1);RK4_4 = RK_4(15,y4,func_vector,0.1)
    euler_11= euler(6,y11,funct_vector1,0.1);RK2_11= RK_2(6,y11,funct_vector1,0.1)
    RK4_11= RK_2(6,y11,funct_vector1,0.1) 
    data1 = {"t":euler_1[0],"x EULER":euler_1[1][0],"y (EULER)":euler_1[1][1],"x (RK2)":RK2_1[1][0],"y(RK2)":RK2_1[1][1],"x (RK4)":RK4_1[1][0],"y (RK4)":RK4_1[1][1]}
    print('FOR 1st CONDITION [t,x,y] = [0,0,-1]')
    print(pd.DataFrame(data1))
    print('FOR 2nd CONDITION [t,x,y] = [0,0,-2]')
    data2 = {"t":euler_2[0],"x EULER":euler_2[1][0],"y (EULER)":euler_2[1][1],"x (RK2)":RK2_2[1][0],"y(RK2)":RK2_2[1][1],"x (RK4)":RK4_2[1][0],"y (RK4)":RK4_2[1][1]}
    print(pd.DataFrame(data2))
    print('FOR 3rd CONDITION [t,x,y] = [0,0,-3]')
    data3 = {"t":euler_3[0],"x EULER":euler_3[1][0],"y (EULER)":euler_3[1][1],"x (RK2)":RK2_3[1][0],"y(RK2)":RK2_3[1][1],"x (RK4)":RK4_3[1][0],"y (RK4)":RK4_3[1][1]}
    print(pd.DataFrame(data3))
    print('FOR 4th CONDITION [t,x,y] = [0,0,-4]')
    data4 = {"t":euler_4[0],"x EULER":euler_4[1][0],"y (EULER)":euler_4[1][1],"x (RK2)":RK2_4[1][0],"y(RK2)":RK2_4[1][1],"x (RK4)":RK4_4[1][0],"y (RK4)":RK4_4[1][1]}
    print(pd.DataFrame(data4))
    print("For Linear Electric Network")
    data5={"t":euler_11[0],"I1(t)(EULER)":euler_11[1][0],"I2(t)(EULER)":euler_11[1][1],"I1(t)(RK2)":RK2_11[1][0],"I2(t)(RK2)":RK2_11[1][1],"I1(t)(RK4)":RK4_11[1][0],"I2(t)(RK4)":RK4_11[1][1]}
    print(pd.DataFrame(data5))
    pd.set_option('display.expand_frame_repr',False)
    pd.set_option('display.max_rows', None)  
    graph(euler_1[0],euler_1[1][0],euler_1[1][1],RK2_1[0],RK2_1[1][0],RK2_1[1][1],RK4_1[0],RK4_1[1][0],RK4_1[1][1],euler_2[0],euler_2[1][0],euler_2[1][1],RK2_2[0],RK2_2[1][0],RK2_2[1][1],RK4_2[0],RK4_2[1][0],RK4_2[1][1],euler_3[0],euler_3[1][0],euler_3[1][1],RK2_3[0],RK2_3[1][0],RK2_3[1][1],RK4_3[0],RK4_3[1][0],RK4_3[1][1],euler_4[0],euler_4[1][0],euler_4[1][1],RK2_4[0],RK2_4[1][0],RK2_4[1][1],RK4_4[0],RK4_4[1][0],RK4_4[1][1])
    graph1(euler_1[1][0],euler_1[1][1],RK2_1[1][0],RK2_1[1][1],RK4_1[1][0],RK4_1[1][1],euler_2[1][0],euler_2[1][1],RK2_2[1][0],RK2_2[1][1],RK4_2[1][0],RK4_2[1][1],euler_3[1][0],euler_3[1][1],RK2_3[1][0],RK2_3[1][1],RK4_3[1][0],RK4_3[1][1],euler_4[1][0],euler_4[1][1],RK2_4[1][0],RK2_4[1][1],RK4_4[1][0],RK4_4[1][1])
    graph_sketch2(euler_11[0],euler_11[1][0],euler_11[1][1],RK2_11[0],RK2_11[1][0],RK2_11[1][1],RK4_11[0],RK4_11[1][0],RK4_11[1][1])
    graph_sketch1(euler_11[1][0],euler_11[1][1],RK2_11[1][0],RK2_11[1][1],RK4_11[1][0],RK4_11[1][1])