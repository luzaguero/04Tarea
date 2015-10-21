#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
G=1
M=1
m=1

class Planeta(object):
    '''
    Simula, con opcion de 3 distintos metodos de integracion (Euler, Runge-Kutta4 y Verlet),
    el movimiento de un planeta que orbita al rededor del Sol, usando como base su pontecial.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un metodo especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual = condicion_inicial
        self.t_actual = 0.
        self.alpha = alpha

    def ecuacion_de_movimiento(self):
        '''
        Implementa la ecuacion de movimiento, como sistema de ecuaciones de
        primer orden.
        '''
        x, y, vx, vy = self.y_actual
        fx=lambda x,y,t: G*M*x*((2*self.alpha)/((x**2 + y**2)**2) - 1/((np.sqrt(x**2 + y**2))**3))
        fy=lambda x,y,t: G*M*y*((2*self.alpha)/((x**2 + y**2)**2) - 1/((np.sqrt(x**2 + y**2))**3))

        return np.array([vx, vy, fx, fy])

    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
        t = self.t_actual
        x, y, vx, vy = self.y_actual

        fx = self.ecuacion_de_movimiento()[2]
        fy = self.ecuacion_de_movimiento()[3]

        fxn = fx(x,y,t)
        fyn = fy(x,y,t)

        self.y_actual = x + dt*(vx + fxn), y + dt*(vy + fyn), vx + dt*fxn, vy + dt*fyn
        pass

    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        '''

        t=self.t_actual
        x,y,vx,vy=self.y_actual
        fx=self.ecuacion_de_movimiento()[2]
        fy=self.ecuacion_de_movimiento()[3]

        k1x=dt*vx; k1y=dt*vy
        fx1=fx(x,y,t)
        fy1=fy(x,y,t)

        k2x=dt*(vx+dt*fx1/2.); k2y=dt*(vy+dt*fy1/2.)
        fx2=fx(x+k1x/2.,y+k1y/2.,t+dt/2.)
        fy2=fy(x+k1x/2.,y+k1y/2.,t+dt/2.)

        k3x=dt*(vx+dt*fx2/2.); k3y=dt*(vy+dt*fy2/2.)
        fx3=fx(x+k2x/2.,y+k2y/2.,t+dt/2.)
        fy3=fy(x+k2x/2.,y+k2y/2.,t+dt/2.)

        k4x=dt*(vx+dt*fx3); k4y=dt*(vy+dt*fy3)
        fx4=fx(x+k3x,y+k3y,t+dt)
        fy4=fy(x+k3x,y+k3y,t+dt)

        xn=x+(k1x+2*k2x+2*k3x+k4x)/6.
        vxn=vx+(fx1+2*fx2+2*fx3+fx4)/6.
        yn=y+(k1y+2*k2y+2*k3y+k4y)/6.
        vyn=vy+(fy1+2*fy2+2*fy3+fy4)/6.

        self.y_actual=xn,yn,vxn,vyn
        pass

    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        t=self.t_actual
        x,y,vx,vy=self.y_actual

        fx=self.ecuacion_de_movimiento()[2]
        fy=self.ecuacion_de_movimiento()[3]

        fxn = fx(x,y,t)
        fyn = fy(x,y,t)

        xn=x+vx*dt+(fxn*(dt**2))/2.
        yn=y+vy*dt+(fyn*(dt**2))/2.

        vxn=vx+((fxn+fx(xn,yn,t+dt))*dt)/2.
        vyn=vy+((fyn+fy(xn,yn,t+dt))*dt)/2.

        self.y_actual=xn,yn,vxn,vyn
        pass

    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        '''
        x, y, vx, vy = self.y_actual
        r = np.sqrt(x**2 + y**2)
        U = - G*M*m/r  + self.alpha*G*M*m /r**2
        E_total = m*(vx**2 +vy**2)/2. + U

        return E_total
        pass


