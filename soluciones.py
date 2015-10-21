#!/usr/bin/env python
# -*- coding: utf-8 -*-

from planeta import Planeta
import numpy as np
import matplotlib.pyplot as plt

'''
CI = condiciones iniciales, tal que:
x_0 = 10
y_0 = 0
vx_0 = 0
vy_0 libre
'''

# - - - Solucion usando Euler - - - #
CI = np.array([10., 0, 0, 0.35])
p = Planeta(CI)

t_total = 1000.
N = 1000
t_paso= t_total/N

x = np.zeros(N); y = np.zeros(N)
vx = np.zeros(N); vy = np.zeros(N)

E = np.zeros(N)

[x[0],y[0],vx[0],vy[0]] = CI

E[0] = p.energia_total()
t= np.linspace(1,t_total,N)

for i in range (1,N):

    p.avanza_euler(t_paso)
    x[i] = p.y_actual[0]
    y[i] = p.y_actual[1]
    vx[i] = p.y_actual[2]
    vy[i] = p.y_actual[3]
    E[i] = p.energia_total()

fig=plt.figure(1)
fig.subplots_adjust(hspace=0.4)

plt.subplot(2, 1, 1)
plt.plot(x, y, 'c-')
plt.title("Trayectoria potencial central (Metodo Euler)")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(2, 1, 2)
plt.plot(t,E, 'c-')
plt.title("Energia orbita")
plt.xlabel("Tiempo")
plt.ylabel("Energia")

plt.show()

# - - - solucion usando RK4 - - - #
CI = np.array([10., 0, 0, 0.35])
p = Planeta(CI)

t_total =  1000.
N = 1000
t_paso= t_total/N

x = np.zeros(N); y = np.zeros(N)
vx = np.zeros(N); vy = np.zeros(N)

E = np.zeros(N)

[x[0],y[0],vx[0],vy[0]] = CI

E[0] = p.energia_total()

for i in range (1,N):

    p.avanza_rk4(t_paso)
    x[i] = p.y_actual[0]
    y[i] = p.y_actual[1]
    vx[i] = p.y_actual[2]
    vy[i] = p.y_actual[3]
    E[i] = p.energia_total()


fig=plt.figure(1)
fig.subplots_adjust(hspace=0.4)

plt.subplot(2, 1, 1)
plt.plot(x, y, 'r-')
plt.title("Trayectoria potencial central (Metodo RK4)")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(2, 1, 2)
plt.plot(t,E, 'r-')
plt.title("Energia orbita")
plt.xlabel("Tiempo")
plt.ylabel("Energia")

plt.show()

# - - - Solucion usando Verlet - - - #
CI = np.array([10., 0, 0, 0.35])
p = Planeta(CI)

t_total = 5000.
N = 1000
t_paso= t_total/N

x = np.zeros(N); y = np.zeros(N)
vx = np.zeros(N); vy = np.zeros(N)

E = np.zeros(N)

[x[0],y[0],vx[0],vy[0]] = CI

E[0] = p.energia_total()

for i in range (1,N):

    p.avanza_verlet(t_paso)
    x[i] = p.y_actual[0]
    y[i] = p.y_actual[1]
    vx[i] = p.y_actual[2]
    vy[i] = p.y_actual[3]
    E[i] = p.energia_total()


fig=plt.figure(1)
fig.subplots_adjust(hspace=0.4)

plt.subplot(2, 1, 1)
plt.plot(x, y, 'm-')
plt.title("Trayectoria potencial central (Metodo Verlet)")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(2, 1, 2)
plt.plot(t,E, 'm-')
plt.title("Energia Orbita")
plt.xlabel("Tiempo")
plt.ylabel("Energia")

plt.show()