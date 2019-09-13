#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:01:39 2019

@author: viblab
"""
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import matplotlib.lines as mline

radius = 3
std_angle = np.array([0,30,60,90,120,150,180])
std_point = np.zeros((std_angle.size,2))
sen_point = np.array([-0.45,-0.15,0.15,0.45])

## estimasi angle (edit disini)
est_angle = np.array([
0,
36.19,
67.899,
90.67,
117.8437,
170.985,
167.72,
])

## estimasi point
est_point = np.zeros((est_angle.size,2))

for i in range(std_angle.size):
      std_point[i,0] = radius*mt.cos(mt.radians(std_angle[i]))
      std_point[i,1] = radius*mt.sin(mt.radians(std_angle[i])) + 1

## calculate estimasi     
for i in range(std_angle.size):
      est_point[i,0] = radius*mt.cos(mt.radians(est_angle[i]))
      est_point[i,1] = radius*mt.sin(mt.radians(est_angle[i])) + 1
      
fig, ax = plt.subplots()

for i in range(std_angle.size):
       txtnum = "%s" % std_angle[i]
       ax.add_patch(plt.Circle((std_point[i,0],std_point[i,1]), 0.1, color='r'))
       plt.text(std_point[i,0],std_point[i,1],txtnum)
       
## plot estinasi
for i in range(std_angle.size):
       txtnum = "%s" % est_angle[i]
       ax.add_patch(plt.Circle((est_point[i,0],est_point[i,1]), 0.1, color='b'))
       plt.text(est_point[i,0],est_point[i,1],txtnum)
       
for i in range(sen_point.size):
      ax.add_patch(plt.Circle((sen_point[i],1), 0.1, color='g'))

#### contoh kode goblok nan males #######       
sqr0=[(-4.5,5.5),(4.5,5.5)]
(line_x0, line_y0) = zip(*sqr0)
ax.add_line(mline.Line2D(line_x0, line_y0, linewidth=2*radius, color='black'))

sqr1=[(-4.5,0),(4.5,0)]
(line_x1, line_y1) = zip(*sqr1)
ax.add_line(mline.Line2D(line_x1, line_y1, linewidth=2*radius, color='black'))

sqr2=[(-4.5,5.5),(-4.5,0)]
(line_x2, line_y2) = zip(*sqr2)
ax.add_line(mline.Line2D(line_x2, line_y2, linewidth=2*radius, color='black'))

sqr3=[(4.5,0),(4.5,5.5)]
(line_x3, line_y3) = zip(*sqr3)
ax.add_line(mline.Line2D(line_x3, line_y3, linewidth=2*radius, color='black'))

################################       

#ax.set_aspect('equal', adjustable='datalim')
ax.plot()
plt.show()