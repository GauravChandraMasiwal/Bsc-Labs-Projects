'''
Roll No.: 2020PHY1122
Name: Gaurav Chandra
'''

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import random


def Generator(N_c, N_t):

    out_arr = []    #outcomes_array

    for i in range(N_t):
        out = []  #outcomes

        for j in range(N_c):
            
            out.append(random.randint(0,1))

        out_arr.append(out)

    n_heads = []  # counting the number of heads and tails

    for n in out_arr:
        n_heads.append(sum(n))

    n_heads = np.array(n_heads)
    n_tails = N_c - n_heads

    freq_count = []  # frequency count of macro-states

    for i in range(N_c + 1):
        freq_count.append(list(n_heads).count(i))

    freq_count = np.array(freq_count)

    probability = freq_count/N_t

    binomial_distri_prob = [] # data from bionomial distribtion
    Nc_arr = np.arange(0, N_c+1)

    for i in range(len(Nc_arr)):
        binomial_distri_prob.append(math.comb(N_c, i)/(2**N_c))

    p = []   # calculating p and q
    for i in range(len(n_heads)): 
        p.append(np.sum(n_heads[:i+1])/((i+1)*N_c))
        
    p = np.array(p)
    q = 1-p

    table = pd.DataFrame({'Trials': np.arange(
        1, N_t + 1), 'Outcomes': out_arr, 'No. of Heads': n_heads, 'No. of Tails': n_tails, 'p': p, 'q': q})

    table.set_index('Trials', inplace=True)

    return table, probability, binomial_distri_prob


# plot1

N_c = 7
N_t = 20
while N_t <= 20000:

    y = Generator(N_c, N_t)[1]
    x = np.arange(0, N_c+1)

    plt.plot(x, y, label=f'N_t = {N_t}', marker='o')
    N_t *= 10

y_bd = Generator(N_c, N_t)[2]
plt.plot(x, y_bd, label='bionomial distribution', marker='*')
plt.legend()
plt.grid()
plt.xlabel('NUMBER OF HEADS')
plt.ylabel('PROBABILITY')
plt.savefig('PLOT1_1122')
plt.title('TRIALS VARIATION PLOT')
plt.show()

# plot2

N_t = 20000
for coins in range(2, 10, 2):

    x = np.arange(0, coins+1)
    y = Generator(coins, N_t)[1]

    plt.plot(x, y, label=f'N_c = {coins}', marker='o')

plt.legend()
plt.grid()
plt.xlabel('NUMBER OF COINS')
plt.ylabel('PROBABILITY')
plt.title('COIN VARIATION PLOT')
plt.savefig('PLOT2_1122')

plt.show()

# plot3

N_c = 3
N_t = 10000

data = Generator(N_c, N_t)[0]
y1 = data['p'].to_numpy()
y2 = data['q'].to_numpy()

x = np.arange(0, N_t)

plt.plot(x, y1, label='p')
plt.plot(x, y2, label='q')
plt.legend()
plt.grid()
plt.xlabel('NUMBER OF TRIALS')
plt.ylabel('p, q')
plt.title('CUMULATIVE PLOT')
plt.savefig('PLOT3_1122')

plt.show()


    

        
    
    
    
        
        
  

    