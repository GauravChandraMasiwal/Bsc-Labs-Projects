#Here i have created functions for linear and weighted curve fitting methods and then i compared the results with the inbuilt linear regression model

import numpy as np
import matplotlib.pyplot as plt
#(a) ’Mylsf’ and ’Mywlsf’ Functions
print("(a) Least Square Fitting and Weighted Least Sqaure Fitting")
def Mylsf(x, y):
    n1 = len(x)
    n2 = len(y)
    slope = 0
    intercept = 0
    if n1 == n2 and n1 > 3:
        sigma_xi = 0
        sigma_yi = 0
        sigma_xi_yi = 0
        sigma_xisq = 0
        res_sum = 0
        res_sum_sq = 0
        
        count = 0
        while count < n1:
            sigma_xi = sigma_xi + x[count]
            sigma_yi = sigma_yi + y[count]
            sigma_xi_yi = sigma_xi_yi + x[count] * y[count]
            sigma_xisq = sigma_xisq + x[count]**2
            count = count + 1
            residual_sum = 0
            residual_sum_sq = 0
            corr_num = 0
            corr_deno1 = 0
            corr_deno2 = 0
    
    delta = (n1 * sigma_xisq - sigma_xi ** 2)  
    m = (n1 * sigma_xi_yi - sigma_xi * sigma_yi) / delta
    c = (sigma_xisq * sigma_yi - sigma_xi * sigma_xi_yi) /delta
    res_sum = res_sum + (n2-intercept*n1-slope)
    res_sum_sq = res_sum_sq + (n2-intercept*n1-slope)**2
    sigma_y = np.sqrt(res_sum_sq)/(n1-2)
    sigma_c = sigma_y*np.sqrt(sigma_xisq / delta)
    for i in range(0,len(x)):
        residual_sum = residual_sum + (y[i] - m*x[i]-c)
        residual_sum_sq = residual_sum_sq + (y[i] - m*x[i] - c)**2
        corr_num = corr_num + (y[i] - np.mean(x))*(y[i] - np.mean(y))
        corr_deno1 =corr_deno1 + (x[i] - np.mean(y))**2
        corr_deno2 = corr_deno2 + (y[i] - np.mean(x))**2

    sigma_y = np.sqrt(residual_sum_sq/(len(x) - 2))
    sigma_c = sigma_y*np.sqrt(sigma_xisq/delta)
    sigma_m = sigma_y * np.sqrt(len(x)/delta)
    correlation_coef = corr_num/(np.sqrt(corr_deno1 * corr_deno2))
    return  sigma_y,m,c,sigma_m ,sigma_c ,residual_sum,residual_sum_sq,correlation_coef

ldata1=np.loadtxt('D:\\python work\\progg class\\1122.dat',delimiter=',',dtype='double').T
xx = ldata1[0]
yy = ldata1[1]
sigma_y,m,c,sigma_m ,sigma_c ,residual_sum,residual_sum_sq,correlation_coef= Mylsf(xx, yy)
print("\n","Slope = ",m ,"\n","Intercept =",c,"\n","Error in Slope=",sigma_m,"\n","Error in Intercept =",sigma_c,"\n","Sum of residuals = ",residual_sum,"\n","Sum of residuals sqaure = ",residual_sum_sq,"\n","Coefficient of Correlation =",correlation_coef)

'''
import matplotlib.pyplot as plt
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.grid()
plt.scatter(xx,yy, c="red")
plt.title("LEAST SQUARE FITTING")

for i_x1,i_y1 in zip(xx,yy):
    plt.text(i_x1,i_y1,'({},{})'.format(i_x1,i_y1))
plt.show()
print()
print()'''

def wlsf(x,y,e):
    n1=len(x)
    n2=len(y)
    slope=0                        
    intercept=0                                        
    w=1/e**2                   #WEIGHT
    i=0
    χ_2 =0                    # χ**2
                     
    if n1==n2 and n1>3: 
        S_wi_xi=0
        S_wi_yi=0
        S_wi_xi_yi=0
        S_wi_xisq=0
        S_wi=0
        x_mean = 0
        y_mean = 0
        S_xx = 0
        S_yy=0
        Sw=sum(w)
        while i < n1: 
                   x_mean += x[i]*w[i]/Sw
                   y_mean += y[i]*w[i]/Sw
                   i = i+1
        count=0
        while count<n1:
            S_wi_xi=S_wi_xi+ x[count]*w[count]
            S_wi_yi= S_wi_yi+y[count]*w[count]
            S_wi_xi_yi=S_wi_xi_yi+ x[count]*y[count]*w[count]
            S_wi_xisq=S_wi_xisq+(x[count]**2)*w[count]
            S_wi=S_wi+w[count]
            S_xx = S_xx + w[count]*(x[count]- x_mean)*2
            S_yy =  S_yy +w[count]*(y[count] - y_mean)*2            
            count+=1

            
        # SLOPE , INTERCEPT AND CORRELATION COEFFICIENT
        intercept=(S_wi_xisq*S_wi_yi-S_wi_xi*S_wi_xi_yi)/(S_wi*S_wi_xisq-S_wi_xi**2)
        slope=(S_wi*S_wi_xi_yi-S_wi_xi*S_wi_yi)/(S_wi*S_wi_xisq-S_wi_xi**2)
        corr_coeff = (slope*np.sqrt(S_xx))/(np.sqrt(S_yy))
        
        
        #DEFINING CORRESPONDING BEST FITTED Y VALUE FOR X
        c=np.array([intercept]*n1)
        xm=np.dot(slope,x)
        y_calc=xm+c                     #Best Fitted y
       
        # ERROR IN SLOPE , INTERCEPT , STANDARD DEVIATION OF PREDICTED Y , CHI**2
        j=0
        while j < n1:
                   χ_2 += w[j]*(y[j] - y_calc[j])*2
                   j = j+1
        e_s=((S_wi)/(S_wi*S_wi_xisq - S_wi_xi*2))*(1/2)
        e_i=((S_wi_xisq)/(S_wi*S_wi_xisq-S_wi_xi*2))*(1/2)
        stdev = np.sqrt((n1*χ_2)/((n1-2)*Sw))

        # DETERMINING NO. OF DATA POINTS MORE THAN 1 SIGMA DEVIATION FROM FITTED LINE
        def my_condition(x,y):
               return abs(x-y)>stdev

            
        print('X (1/\u03BB\u00b2 ) = ', x)
        print('Y (\u03BC) = ',y)
        print('Y_calculated (\u03BC) =  ', y_calc)
        print('SLOPE =  ',slope)
        print('INTERCEPT=  ',intercept)
        print('CORRELATION COEFFICIENT =  ' , corr_coeff)
        print('ERROR IN SLOPE =  ',e_s)
        print('ERROR IN INTERCEPT =  ',e_i)
        print('χ SQAURED =  ',χ_2)
        print('χ =  ',np.sqrt( χ_2) )
        
'''
        # PLOTTING FITTED LINE AND ERROR BARS 
        plt.ylim(1.51511934,1.53909402)
        plt.scatter(x,y, c = 'NAVY',s= 70 ,edgecolors='black',alpha=0.95,label = 'observed data')
        plt.plot(x,y_calc, c='PURPLE', linewidth = 1.5, label = 'Fitted line')
        plt.errorbar(x,y,yerr=e,fmt = ' ',c = 'DARKGREEN',label = 'Errors asssociated with y')
        plt.legend()
        plt.grid(True)
        plt.title("WEIGHTED LEAST SQUARE FIT FOR CAUCHY'S CONSTANTS",c='DEEPPINK',fontsize=20,fontweight='bold')
        plt.xlabel(' WAVELENGTH INVERSE SQAURED 1/\u03BB\u00b2 (1/m\u00b2) ' ,c='DODGERBLUE',fontsize=15)
        plt.ylabel(' REFRACTIVE INDEX (\u03BC) ',c='DODGERBLUE',fontsize=15)
        plt.show()'''



y = np.mean(ldata1[1])
a = []
my_l= 1/(np.std(yy))**2
a.append(my_l)
e=np.array(my_l)
wlsf(xx,np.mean(yy),Mylsf(xx, yy)[0])


#(f) Comparison of Results from the inbuilt function of Lingress
from scipy.stats import linregress
print("(f) Comparison with Inbuilt Function")
print(linregress(xx, yy))