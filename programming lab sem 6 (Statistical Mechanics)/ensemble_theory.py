#NAME : GAURAV CHANDRA
#ROLLNO : 2020PHY1122

import numpy as np
import matplotlib.pyplot as plt
import random


def microstates(n_t,n_c):
    Arr = []
    i = 0
    while i<n_c:
        r = random.randint(0,1)
        Arr.append(r)
        i = i+1
    Arr = np.array(Arr)
    j = 1
    while j<n_t:
        states = np.zeros(n_c)
        i = 0
        while i<n_c:
            r = random.randint(0,1)
            states[i] = r
            states = np.array(states)
            i = i+1
        data = np.vstack([Arr,states])
        arr = data
        j = j+1    
     
    return data    


def add(arr):
    Arr2 = []
    for i in Arr2:
        sum1 = 0
        for j in i:
            sum1 = sum1 + j
        Arr2.append(sum1)            
    return Arr2        
    
    
def count(arr,n_t,n_c):
    count = []
    for i in range(n_c+1):
        c = arr.count(i)
        count.append(c)
    freq = []
    for i in count:
        freq.append(i/n_t)    
    return freq


#For fixed n_c and variable n_t
n_t = [10,100,1000,2000,3000]
n_c = 3
n_h = np.arange(0,n_c+1,1)
Freq = []
for i in n_t:
    data_a = microstates(i,n_c)
    data_b = add(data_a)
    freq = count(data_b,i,n_c)
    Freq.append(freq)

for i in range(len(Freq)):
    plt.plot(n_h,Freq[i],label = "No of trials "+str(n_t[i]))
    plt.scatter(n_h,Freq[i])
plt.xlabel("NO. OF HEADS")
plt.ylabel("FREQUENCY")    
plt.title("NO. OF COINS = "+str(n_c))
plt.grid()
plt.legend()    
plt.show()


#For fixed n_t and variable n_c
n_t = 100
n_c = np.arange(1,11,1)
n_h = [0]

for i in n_c:
    n_h.append(i)
Freq = []
for i in n_c:
      data_a = microstates(n_t,i)
      data_b = add(data_a)
      freq = count(data_b,n_t,10)
      Freq.append(freq)

for i in range(len(Freq)):
    j = i+1
    plt.plot(n_h,Freq[i],label = "NO. OF COINS "+str(j))
    plt.scatter(n_h,Freq[i])
    
    
plt.xlabel("NO. OF HEADS")
plt.ylabel("FREQUENCY")
plt.title("NO. OF TRIALS = "+str(n_t))
plt.legend()
plt.grid()
plt.show()


def cumulative(n_t,n_c):
    mic_states = microstates(n_t,n_c)
    heads = []             #1
    tails = []             #0
    for i in mic_states:
        k = []
        for items in i:
            k.append(items)
        h = k.count(1)
        t = k.count(0)
        heads.append(h)
        tails.append(t)
    p = []   #heads
    q = []   #tails
    k = 0   
    h = 0
    t = 0
    for i in range(len(heads)):
        p.append((heads[i]+h)/(heads[i]+tails[i]+k))
        q.append((tails[i]+t)/(heads[i]+tails[i]+k))
        k = k + heads[i]+tails[i]
        h = heads[i] + h
        t = tails[i] + t   
    return p,q


n_t = 30000
N_t = np.arange(1,n_t+1,1)
n_c = 3
p,q = cumulative(n_t,n_c)
plt.plot(N_t,p,label = "p (HEADS)")
plt.plot(N_t,q,label = "q (TAILS)")
plt.scatter(N_t,p,marker=".")
plt.scatter(N_t,q,marker=".")
plt.title("FOR NO. OF COINS = "+str(n_c))
plt.xlabel("N0. OF TRIALS")
plt.ylabel("CUMULATIVE FREQUENCY")
plt.legend()
plt.grid()
plt.show()





