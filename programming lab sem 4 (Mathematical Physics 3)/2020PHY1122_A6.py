#NAME : GAURAV CHANDRA
#ROLLNO : 2020PHY1122
#PARTNER NAME : KHUSHI
#PARTNER ROLLNO : 2020PHY1155

#curvefiiting graph
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.stats import linregress 

#fnct for leastsquarefitting and calculations of other constants
    
def Mylsf(x,y):            
    Y_CAL = []
    E = []
    sum_of_residual = 0
    sum_of_residual_squares = 0
    a=len(x)
    b=len(y)
    a0 = 0
    a1 = 0
    if a==b and a>2:
        sum_xi = 0
        sum_yi = 0
        sum_xisq = 0
        sum_xiyi = 0
        sum_yisq = 0

        for i in range(a):
            sum_xi = sum_xi + x[i]
            sum_yi = sum_yi + y[i]
            sum_xisq = sum_xisq + x[i]**2
            sum_xiyi = sum_xiyi + x[i]*y[i]
            sum_yisq = sum_yisq + y[i]**2

        a0=(sum_xi * sum_yi - a * sum_xiyi)/(sum_xi **2 - a*sum_xisq)              #slope

        a1=(sum_xisq * sum_yi - sum_xi * sum_xiyi)/(a * sum_xisq - sum_xi**2)      #intercept
        
        bxy = sum_xiyi/sum_yisq   #regression coefficient of TENSION on LAMBDA SQUARE
        
        byx = sum_xiyi/sum_xisq   #regression coefficient of LAMBDA SQUARE on TENSION
        
        r = np.sqrt(bxy*byx)   #correlation coeffcient
        
        for i in range(a):
            Ycal = x[i]*a0 + a1
            Y_CAL.append(Ycal)
            error = Ycal - y[i]
            E.append(error)
            sum_of_residual += error   #sum of residuals
            sum_of_residual_squares += error**2  #sum of residual squares
        
        SSxx = sum_xisq - ((sum_xi)**2/a)   
        std_slope = np.sqrt((sum_of_residual_squares)/(SSxx * (a-2)))  #standard deviation of slope
        std_intercept = np.sqrt(((std_slope)**2)*(sum_xisq/a))   #standard deviation of intercept             
        #delta = (a * sum_xisq - sum_xi**2)  
        #sigma_y = np.sqrt(sum_of_residual_squares/(a - 2))
        #sigma_c = sigma_y*np.sqrt(sum_xisq/delta)
        #sigma_m = sigma_y*np.sqrt(a/delta)
    else:
        print(" error ! check observations again")

    return a0,a1,sum_of_residual,sum_of_residual_squares,std_slope,std_intercept,r,Y_CAL


def Mywlsf(x,y,w):
    Y_CAL = []
    E = []
    sum_of_residual = 0
    sum_of_residual_squares = 0
    a=len(x)
    b=len(y)
    m = 0
    c = 0
    i = 0
    if a==b and a>2:
        sum_wi = 0
        sum_wixi = 0
        sum_wiyi = 0
        sum_wixiyi = 0
        sum_wixisq = 0
        delta = 0
        #sum_wi = sum(w)
        x_mean,y_mean = 0,0
        S_xx,S_yy = 0,0
        for i in range(a):
            sum_wi = sum_wi + w[i]
        while i < a: 
                   x_mean += x[i]*w[i]/sum_wi
                   y_mean += y[i]*w[i]/sum_wi
                   i = i+1
        for i in range(a):
            
            sum_wixi = sum_wixi + w[i]*x[i]
            sum_wiyi = sum_wiyi + w[i]*y[i]
            sum_wixiyi = sum_wixiyi + w[i]*x[i]*y[i]
            sum_wixisq = sum_wixisq + w[i]*(x[i]**2)
            S_xx = S_xx + w[i]*(x[i]- x_mean)*2
            S_yy =  S_yy +w[i]*(y[i] - y_mean)*2  
        
        delta = delta + sum_wi*sum_wixisq - (sum_wixi)**2
        
        m = ((sum_wi * sum_wixiyi) - (sum_wixi * sum_wiyi))/delta
        c = ((sum_wixisq * sum_wiyi) - (sum_wixi * sum_wixiyi))/delta
        corr_coeff = (m*np.sqrt(S_xx))/(np.sqrt(S_yy))
        
        for i in range(a):
            Ycal = x[i]*m + c
            Y_CAL.append(Ycal)
            error = Ycal - y[i]
            E.append(error)
            sum_of_residual += error   #sum of residuals
            sum_of_residual_squares += error**2  #sum of residual squares
        
        std_slope = np.sqrt(sum_wi/delta)
        std_intercept = np.sqrt(sum_wixisq/delta)
                
    else :
        print("check observations again !!!")
    return m,c,sum_of_residual,sum_of_residual_squares,std_intercept,std_slope,corr_coeff,Y_CAL


ldata1=np.loadtxt('D:\\python work\\progg class\\1122.dat',delimiter=',',dtype='double').T
xx = ldata1[0]
yy = ldata1[1]

M=[];T1=[];T2=[];T3=[];T4=[];T5=[];T6=[];T7=[];T8=[];T9=[];T10=[]   
for i in range(len(ldata1[0])):             #to convert elements of list from str into float
    M.append(float(ldata1[0][i]))
    T1.append(float(ldata1[1][i]**2))
    T2.append(float(ldata1[2][i]**2))
    T3.append(float(ldata1[3][i]**2))
    T4.append(float(ldata1[4][i]**2))
    T5.append(float(ldata1[5][i]**2))
    T6.append(float(ldata1[6][i]**2))
    T7.append(float(ldata1[7][i]**2))
    T8.append(float(ldata1[8][i]**2))
    T9.append(float(ldata1[9][i]**2))
    T10.append(float(ldata1[10][i]**2))
    
T=[]
std_err = []
for i in range(len(ldata1[0])):
    lis = [T1[i],T2[i],T3[i],T4[i],T5[i],T6[i],T7[i],T8[i],T9[i],T10[i]]
    s1=np.mean(lis)
    s2=np.std(lis)  
    T.append(s1)
    std_err.append(s2)
std_err = np.array(std_err)
w_yy = 1/(std_err)**2

print("putput using ordinary least square method:")
m,c,residual_sum,residual_sum_sq,sigma_m,sigma_c,correlation_coef,y_cal= Mylsf(xx, T)
print("\n","Slope = ",m ,"\n","Intercept =",c,"\n","Error in Slope=",sigma_m,"\n","Error in Intercept =",sigma_c,"\n","Sum of residuals = ",residual_sum,"\n","Sum of residuals sqaure = ",residual_sum_sq,"\n","Coefficient of Correlation =",correlation_coef)

print()
print("putput using weighted least square method:")

M,C,Residual_sum,Residual_sum_sq,Sigma_m,Sigma_c,Correlation_coef,Y_CAL= Mywlsf(xx,T,w_yy)
print("\n","Slope = ",M ,"\n","Intercept =",C,"\n","Error in Slope=",Sigma_m,"\n","Error in Intercept =",Sigma_c,"\n","Sum of residuals = ",Residual_sum,"\n","Sum of residuals sqaure = ",Residual_sum_sq,"\n","Coefficient of Correlation =",Correlation_coef)

plt.errorbar(xx,y_cal,yerr=std_err,xerr=None,fmt='o',ecolor = 'red',color='black')
plt.scatter(xx,T,label = 'data')
plt.plot(xx,y_cal,label="OLS",c='yellow')
plt.plot(xx,Y_CAL,label="WLS",linestyle="dashed",c='g')
plt.xlabel(r"Mass (g)")
plt.ylabel(r"$T^2$ ($s^2$)")
plt.title("OLS vs WLS")
plt.legend()
plt.grid()
plt.savefig("1122graph.pdf")
plt.show()

#part e

print()
print("calculation of k and m :")
print()
print("slope and intercept from wlsqf are :",M,"and ",C)
print("k = 4*pi*pi/slope = ",4*(np.pi)**2/M,'N/m')
print("m = intercept/slope = ",C/M,'grams')
print()

#part f 
print("OUTPUT FROM LINREGRESS")

print(linregress(xx,T))
print()
print("comparison between lsqf and linregress : ")
print()
print("the absoulte error in slope between both methods is : ",abs(m-linregress(xx,T)[0]))
print("the absoulte error in intercept between both methods is : ",abs(c-linregress(xx,T)[1]))
print("the absoulte error in correlation coefficient between both methods is : ",abs(correlation_coef-linregress(xx,T)[2]))
