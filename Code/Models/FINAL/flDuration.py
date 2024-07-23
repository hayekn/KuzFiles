from final import *
import ruptures as rpt
from scipy import integrate
grad = 1500
W = np.linspace(0, 5, 200)
t = np.linspace(0, 120, grad)
A = []
for w in W:
    y = xppaut_model(t, False, nacTOsetp=w)['Int'][6]
    d = np.gradient(y)
    algo = rpt.Dynp(model="l1", min_size=1, jump=.5).fit(d)
    result = algo.predict(n_bkps=1)[0]
    A.append(result)
A = [a*(120/grad) for a in A]
plt.rcParams.update({
        'font.size': 26
    })
fig, ax = plt.subplots(figsize=(10, 8))
n= 24
plt.xticks(fontsize=n)
plt.yticks(fontsize=n)
plt.ylim(0, 35)
ax.plot(W, A, color='black', linewidth=4)
ax.set_xlabel(r"$\bm{w_{N \to P}}$", labelpad=10)
ax.set_ylabel("Front-Loading Duration (min)", labelpad=10)
plt.savefig("Code/Models/FINAL/Pics/flDur.png", transparent=True, dpi=400)
