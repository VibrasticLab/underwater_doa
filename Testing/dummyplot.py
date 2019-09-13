#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 21:15:13 2019

@author: viblab
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-1.25,0.01,1.25)
m = -0.66
b = 1
y1 = m*x + b
p = np.array([0.4,0.9,1,1.5,2])

fi=0

for i in range(len(p)):
      val = np.power(np.abs(y1),p[i])+np.power(np.power(np.abs(x),p[i]),(1/p[i]))
      Np = np.min(val)
      ind = np.argmin(val)
      xp = np.linspace(-Np,Np,1000)
      yp = (Np**p[i] - np.power(np.power(np.abs(xp),p[i]),1/p[i]))
      
      plt.figure(fi)
      fi=fi+1
      plt.plot(xp,yp,'-k',xp,-yp,'-k',x,y1,'-r')
      
      