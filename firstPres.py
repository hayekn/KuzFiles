from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scint

T = np.linspace(0, 2, 200)
S = np.linspace(0, 2, 15)
def sig(t):
  return 1/(1+e**t)

fig, ax = plt.subplots(figsize=(5, 5))
T = np.linspace(-10, 10, 100)
ax.set_xlabel("t")
ax.plot(T, [sig(t) for t in T])
plt.axvline(0, -1, 1, color='k')
plt.show()

def y1(y2, W12, D1, S1):
  return sig(S1*(W12*y2-D1))
def y2(y1, W21, D2, S2):
  return sig(S2*(W21*y1-D2))

def vecFieldY1(W12, W21, D1, D2, S1, S2):
  return [[-y1+sig(S1*(W12*y2-D1)) for y1 in S] for y2 in S]
def vecFieldY2(W12, W21, D1, D2, S1, S2):
  return [[-y2+sig(S2*(W21*y1-D2)) for y1 in S] for y2 in S]

def odeTraj(t, Y, W12, W21, D1, D2, S1, S2):
  y1 = Y[0]
  y2 = Y[1]

  dy1dt = -y1 + sig(S1*(W12*y2-D1))
  dy2dt = -y2+sig(S2*(W21*y1-D2))

  return [10*dy1dt, 10*dy2dt]

def findInt(orange, blue):
  A = []
  for i in orange:
    for j in blue:
      if abs(i[0]-j[0]) < .005 and abs(i[1] - j[1]) < .005:
        #print("("+str(i[0])+","+str(i[1])+") SAME AS ("+str(j[0])+","+str(j[1])+")")
        A.append([i[0], i[1]])
  return A

def executeStuff(W12, W21, D1, D2, S1, S2,n):
  fig, ax = plt.subplots(figsize=(5, 5))
  T = np.linspace(0, 1.2, 200)
  S = np.linspace(0, 1.2 , 15)
  

  ax.plot(T, [y2(t, W12, D1, S1) for t in T], color='orange', label = 'd$y_2$/dt = 0')
  ax.plot([y1(t, W21, D2, S2) for t in T], T, color='teal',label = 'd$y_1$/dt = 0')
  plt.quiver(S, S, vecFieldY1(W12, W21, D1, D2, S1, S2), vecFieldY2(W12, W21, D1, D2, S1, S2), width=.002, headwidth=5)


  Y1 = [.2, 1.7]
  Y2 = [1, 2]
  Y3 = [.75, 1.1]

  Y4 = [1.7, .2, ]
  Y5 = [2, 1]
  Y6 = [1.1, .75]

  Y7 = [1,1]

  sol1=scint.solve_ivp(odeTraj, [0, 2], Y1, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol1.y[0,:], sol1.y[1,:], 'k-', linewidth=2)
  sol2=scint.solve_ivp(odeTraj, [0, 2], Y2, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol2.y[0,:], sol2.y[1,:], 'k-', linewidth=2)
  sol3=scint.solve_ivp(odeTraj, [0, 2], Y3, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol3.y[0,:], sol3.y[1,:], 'k-', linewidth=2)
  plt.plot(sol3.y[0,-1],sol3.y[1,-1],'ro', markersize = 10)

  sol4=scint.solve_ivp(odeTraj, [0, 2], Y4, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol4.y[0,:], sol4.y[1,:], 'k-', linewidth=2)
  sol5=scint.solve_ivp(odeTraj, [0, 2], Y5, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol5.y[0,:], sol5.y[1,:], 'k-', linewidth=2)
  sol6=scint.solve_ivp(odeTraj, [0, 2], Y6, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol6.y[0,:], sol6.y[1,:], 'k-', linewidth=2)
  plt.plot(sol6.y[0,-1],sol6.y[1,-1],'ro', markersize = 10)


  sol7=scint.solve_ivp(odeTraj, [0, 2], Y7, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
  plt.plot(sol7.y[0,:], sol7.y[1,:], 'k-', linewidth=2)
  plt.plot(sol7.y[0,-1],sol7.y[1,-1],'ro', markersize = 10)
  plt.xlim(0,1.2)
  plt.ylim(0,1.2)
  ax.set_xlabel("$y_1$")
  ax.set_ylabel("$y_2$")
  plt.legend(loc ='lower left')
  slope = str(n)
  plt.savefig(''+slope+'')

  # orange = [[t, y2(t, W12, D1, S1)] for t in T]
  # blue = [[y1(t, W21, D2, S2), t] for t in T]
  # plt.scatter([findInt(orange, blue)[n][0] for n in range(len(findInt(orange, blue)))],
  # [findInt(orange, blue)[n][1] for n in range(len(findInt(orange, blue)))],
  #           color='r', zorder=2)
  plt.show()

executeStuff(.7, .7, .6, .6, 20, 20,1) #(W12, W21, D1, D2, S1, S2)

def slideshow():
  count=-1
  for i in np.linspace(20, 1, 150):
    count += 1
    executeStuff(1, 1, 1, 1, i, i, count)

slideshow()
