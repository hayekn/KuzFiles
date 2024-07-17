from final import *
import ruptures as rpt
from scipy import integrate

def test():
    algo = rpt.Dynp(model="l2", min_size=1, jump=1).fit(d)
    result = algo.predict(n_bkps=1)
    print(result)
    rpt.display(d, result)
def loadChangeArray():
    A = []
    C = np.linspace(1, top, grad)
    for c in C:
        y = xppaut_model(t, False, csTOvta=c)['Int'][6]
        d = np.gradient(y)
        algo = rpt.Dynp(model="l1", min_size=1, jump=.5).fit(d)
        result = algo.predict(n_bkps=1)[0]
        A.append(result)
        print(result)
    with open("Code/Models/FINAL/changePoint.npy", "wb") as f:
        np.save(f, A)

def loadSlopeArray():
    S = []
    with open("Code/Models/FINAL/changePoint.npy", "rb") as f:
            A = np.load(f)

    C = np.linspace(1, top, grad)
    for i in range(np.size(A)):
        y = xppaut_model(t, False, csTOvta=C[i])['Int'][6]
        sec1 = (y[A[i]]-y[0])/A[i]
        sec2 = (y[-1] - y[A[i]])/(np.size(A) - A[i])
        S.append(sec1/sec2)
        
    return S

def makeGraphs():
    with open("Code/Models/FINAL/changePoint.npy", "rb") as f:
            A = np.load(f)
    A = [a*(120/grad) for a in A]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))

    C = np.linspace(1, top, grad)
    M = []
    I = []
    for c in C:
        vta = xppaut_model(t, False, csTOvta=c)['Int'][5]
        M.append(np.max(vta))
        I.append(integrate.simpson(vta, x=t))

    axs[0].plot(I, A, color='purple')
    axs[0].set_ylabel("Change Point (min)")
    axs[0].set_xlabel("Total DA Release")
    axs[1].plot(I,  loadSlopeArray())
    axs[1].set_ylabel("Maintenance Coefficient")
    axs[1].set_xlabel("Total DA Release")

    axs[0].set_ylim(0, 45)
    axs[1].set_ylim(0, 10)
    plt.show()
    return
    axs[0,0].set_ylabel("Change Point (min)")
    axs[0,0].set_xlabel("Sensitivity of VTA")
    axs[0,0].plot(C, A)

    axs[0,1].plot(M, A, color='red')
    axs[0,1].set_ylabel("Change Point (min)")
    axs[0,1].set_xlabel("Peak DA Level")

    axs[0,2].plot(I, A, color='purple')
    axs[0,2].set_ylabel("Change Point (min)")
    axs[0,2].set_xlabel("Total DA Release")

    S = loadSlopeArray()
    axs[1, 0].plot(C, S)
    axs[1,0].set_ylabel("Maintenance Coefficient")
    axs[1,0].set_xlabel("Sensitivity of VTA")

    axs[1, 1].plot(M, S)
    axs[1,1].set_ylabel("Maintenance Coefficient")
    axs[1,1].set_xlabel("Peak DA Level")

    axs[1,2].plot(I, S)
    axs[1,2].set_ylabel("Maintenance Coefficient")
    axs[1,2].set_xlabel("Total DA Release")

    for i in range(3):
         axs[1,i].set_ylim(0, 10)
         axs[0,i].set_ylim(0, 45)

    '''y = xppaut_model(t, False, csTOvta=1)['Int'][6]
    alc = axs[1,1].plot(t, y, label='Alcohol Vol.', color = 'green')[0]
    dot = axs[1,1].plot((A[0], A[0]), (0, 30), '--', color='gray')[0]

    frames=grad
    def update(frame):
        y = xppaut_model(t, False, csTOvta=2*(frame/frames)+1)
        alc.set_data(t, y['Int'][6])
        dot.set_data((A[frame], A[frame]), (0, 30))

    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=1)'''
    plt.tight_layout()
    plt.show()


grad = 1500
top = 3
t = np.linspace(0, 120, grad)

#loadChangeArray()
makeGraphs()

#loadArray()
#makeGraphs()