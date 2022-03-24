import numpy as np
from math import *
import matplotlib.pyplot as plt


#Dans la suite r est la première coordonnée polaire

def fonction(r):
    x=(-887384455.8905718+(-1.9852204086241444e+31/(r**2))+(2.6549172e+20/r))**(1/2) #x=r'
    return x


def rungekutta4(f, r0, t0, tf, nb): #(fonction à intégrer; r initial; t initial; t final; nombre de points)
    T=np.linspace(t0,tf,nb)
    n = len(T)
    R = np.zeros(n)
    R[0] = (r0)
    h=(T[len(T)-1]-T[0])/n
    for i in range(n-2):
        k1 = f(R[i])
        k2 = f(R[i] + k1 * h / 2)
        k3 = f(R[i] + k2 * h / 2)
        k4 = f(R[i] + k3 * h)
        R[i+1] = R[i]+(h / 6) * (k1 + 2*k2 + 2*k3 + k4)
    R[len(R)-1]=152093407000
    return R,T


R = list(rungekutta4(fonction, 147091144000.001053, 0, 15778800, 10000))[0]
T = list(np.linspace(0,1,10000))


def m_en_km(T):
    for i in range(len(T)):
        T[i] = T[i]/1000
    return T
    



P = 10000*[147091144000]
A = 10000*[152093407000]


COS=[]
for i in range(len(T)):
    COS.append(149592275500-2501131500*cos(pi * T[i]))


plt.figure()
plt.xlim(0,1)
plt.title("RK4 (rouge)- Première approximation (bleu)")
plt.xlabel('temps(demi-année)')
plt.ylabel('Distance T- S (m)')
plt.plot(T,COS, color='blue')
plt.plot(T,R, color='red')
plt.plot(T, P, '--', color='black')
plt.plot(T,A, '--', color='black')
plt.show()


#le reste du programme n'a pas d'interêt, il reprend le même principe