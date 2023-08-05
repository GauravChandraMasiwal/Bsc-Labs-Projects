
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 2(a)----------------------------------------------------------------------------
def f(x):
    return 2*(e - V(x))

def V(x):
    return (x**2)/2

def Numerov(a, b, N, f, n):
    h = (b-a)/N
    x = np.linspace(a, b, N+1)

    u = np.zeros([N+1])

    c = 1 + ((h**2)/12)*f(x)

    if n % 2 == 0:  # even states
        u[0] += 1
        u[1] += (6 - 5*c[0])*u[0]/c[1]

    else:  # odd states
        u[0] += 0
        u[1] += h

    for i in range(1, N):
        u[i+1] += ((12 - 10*c[i])*u[i] - c[i-1]*u[i-1])/c[i+1]

    return x, u

e = 0+0.5+10**(-6)
a,b = Numerov(0, 4, 120, f, 0)
plt.plot(a,b)
'''

def simps(x, y):
    h = (x[-1]-x[0])/len(x)
    integral = (h/3)*(2*np.sum(y[2:-2:2]) + 4*np.sum(y[1:-1:2]) + y[0] + y[-1])

    return integral


def V(x):
    return (x**2)/2


def wfn_plot(n, de=None, p=None):

    if de != None:
        del_e = np.array([de])

    else:
        del_e = np.array([10**(-2), 10**(-4), 10**(-6), 10**(-8)])

    e = n + del_e + 0.5

    for i in range(len(e)):

        def f(x):
            return 2*(e[i] - V(x))

        x = Numerov(0, 4, 100, f, n)[0]
        u = Numerov(0, 4, 100, f, n)[1]

        u = u/np.sqrt(simps(x, u**2))  # normalizing the wave fnc

        def parity(x, u):
            x_ = -x[1:]
            X = np.append(x_[-1::-1], x)

            if n % 2 != 0:  # odd states
                u_ = -u[1:]
                U = np.append(u_[-1::-1], u)

            else:  # even states
                u_ = u[1:]
                U = np.append(u[-1::-1], u_)

            return X, U

        if p != None:
            wfn = (parity(x, u)[1])**p
        else:
            wfn = (parity(x, u)[1])

        x_pts = parity(x, u)[0]

        plt.plot(x_pts, wfn, label=f'e = {e[i]}')
        plt.title(f'for n = {n}')

    plt.legend()
    plt.grid()
    plt.xlabel('x')
    if p == None:
        plt.ylabel('Ψ(x)')
    else:
        plt.ylabel(r'$Ψ^{p}(x)$'.format(p=p))

    plt.show()


wfn_plot(0)  # ground state

wfn_plot(1, 10**(-6))  # first excited state
wfn_plot(2, 10**(-6))  # second excited state
wfn_plot(3, 10**(-6))  # third excited state


# 2(b)----------------------------------------------------------------------------

wfn_plot(1, 10**(-6), 2)  # probability density for the first excited state
wfn_plot(2, 10**(-6), 2)  # probability density for the second excited state
wfn_plot(3, 10**(-6), 2)  # probability density for the third excited state

# 2(c)----------------------------------------------------------------------------

h_cut = 1.054571817*(10**(-34))  # in Js
omega = 5.5*(10**(14))

n = np.array([0, 1, 2, 3])
E_ana = (n + 0.5)*h_cut*omega

E_cal = []

for i in n:
    E_cal.append((0.5 + i + 10**(-6))*h_cut*omega)

E_cal = np.array(E_cal)

table = pd.DataFrame({'e_ana': E_ana, 'e_cal': E_cal})
print(table)
'''