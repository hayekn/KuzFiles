from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scint

T = np.linspace(0, 2, 200)
S = np.linspace(0, 2, 15)

def sig(t):
  return 1/(1+e**t)
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

  return [3*dy1dt, 3*dy2dt]

def findInt(orange, blue):
  A = []
  for i in orange:
    for j in blue:
      if abs(i[0]-j[0]) < .0055 and abs(i[1] - j[1]) < .0055:
        #print("("+str(i[0])+","+str(i[1])+") SAME AS ("+str(j[0])+","+str(j[1])+")")
        A.append([i[0], i[1]])
  return [A[0], A[len(A)-1]]

def executePlot(W12, W21, D1, D2, S1, S2, count):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.set_xlabel("y1")
    ax.set_ylabel("y2")

    ax.plot(T, [y2(t, W12, D1, S1) for t in T], color='orange', label='dy2/dt = 0')
    ax.plot([y1(t, W21, D2, S2) for t in T], T, color='teal', label='dy1/dt = 0')
    ax.legend()
    plt.quiver(S, S, vecFieldY1(W12, W21, D1, D2, S1, S2), vecFieldY2(W12, W21, D1, D2, S1, S2), width=.002, headwidth=5)


    Y1 = [.2, 1.7]
    Y2 = [1, 2]
    Y3 = [.75, 1.1]

    Y4 = [1.7, .2, ]
    Y5 = [2, 1]
    Y6 = [1.1, .75]


    sol1=scint.solve_ivp(odeTraj, [0, 2], Y1, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    sol2=scint.solve_ivp(odeTraj, [0, 2], Y2, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    sol3=scint.solve_ivp(odeTraj, [0, 2], Y3, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    sol4=scint.solve_ivp(odeTraj, [0, 2], Y4, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    sol5=scint.solve_ivp(odeTraj, [0, 2], Y5, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    sol6=scint.solve_ivp(odeTraj, [0, 2], Y6, t_eval=T, args=(W12, W21, D1, D2, S1, S2))
    SOLS = [sol1, sol2, sol3, sol4, sol5, sol6]

    for sol in SOLS:
        plt.plot(sol.y[0,:], sol.y[1,:], 'k-', linewidth=1.7)
        plt.arrow(sol.y[0, 25], sol.y[1, 25], sol.y[0, 25]-sol.y[0, 24], sol.y[1, 25]-sol.y[1, 24], shape='full', lw=0, length_includes_head=True, head_width=.05, color='black')


    orange = [[t, y2(t, W12, D1, S1)] for t in T]
    blue = [[y1(t, W21, D2, S2), t] for t in T]
    '''plt.scatter([findInt(orange, blue)[n][0] for n in range(len(findInt(orange, blue)))],
    [findInt(orange, blue)[n][1] for n in range(len(findInt(orange, blue)))],
            color='r', zorder=2)'''
    plt.savefig("Downloads/dynSlides/"+str(count)+"noDotRev.png", dpi=750)
    #plt.show()

def slideshow():
  count=-1
  for i in np.linspace(20, 1, 150):
    count += 1
    executePlot(1, 1, 1, 1, i, i, count)

slideshow()
