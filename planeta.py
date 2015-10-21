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
        pass

    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        '''
        pass

    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        pass

    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        '''
        pass


