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


plt.plot(x, y)
plt.title("Trayectoria potencial central (Metodo Euler)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

plt.plot(t,E)
plt.title("Energia orbita")
plt.xlabel("Tiempo")
plt.ylabel("Energia")

plt.show()
