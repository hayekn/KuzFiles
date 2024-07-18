from final import *
from matplotlib import colors

def generate(reso):
    B = []
    W = np.linspace(0, 6, reso)
    t = np.linspace(0, 2000, 800)
    for SB in W:
        for BS in W:
            y = xppaut_model(t, seekTObin=SB, binTOseek=BS)['Int']
            if y[1][-1] > .95 and y[2][-1] > .95:
                B.append(1)
            elif y[1][round((12/2000)*800)] < .1:
                B.append(-1)
            else:
                B.append(0)
            print(str(B[-1]) + " seekTObin=" + str(SB) + " and binTOseek="+ str(BS))
            print(np.size(B)/(reso**2))

    with open("Code/Models/FINAL/cortical.npy", "wb") as f:
                B = np.array(B).reshape((reso, reso))
                np.save(f, B)

def colorGraphs():
    with open("Code/Models/FINAL/cortical.npy", "rb") as f:
        B = np.load(f)

    cmap = colors.ListedColormap(['white', 'gray', 'black'])

    fig, ax = plt.subplots(figsize=(8, 10))

    cax = ax.imshow(B, extent=[0, 6, 0, 6], origin='lower', cmap=cmap, aspect='auto')
    cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal', shrink=.5)
    cbar.ax.set_xticklabels(['FL Unengaged', 'Two-Phased', 'FL Always Engaged'])

    plt.xlabel('binTOseek')
    plt.ylabel('seekTObin')
    plt.title('Binge Characteristic')
    plt.show()

def linGraphs(reso):
    B = []
    W = np.linspace(0, 6, reso)
    t = np.linspace(0, 2000, 800)
    for w in W:
        y = xppaut_model(t, seekTObin=w, binTOseek=w)['Int']
        if y[1][-1] > .95 and y[2][-1] > .95:
            B.append(1)
        elif y[1][round((12/2000)*800)] < .1:
            B.append(-1)
        else:
            B.append(0)
        print(str(B[-1]) + " seekTObin=" + str(w) + " and binTOseek="+ str(w))
        print(np.size(B)/reso)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(W, B, s=.5, color='black')
    plt.yticks([-1, 0, 1], ['FL Unengaged', 'Two-Phased', 'FL Always Engaged'])
    plt.xlabel("Seek <--> Binge Strength")
    plt.show()

linGraphs(150)

