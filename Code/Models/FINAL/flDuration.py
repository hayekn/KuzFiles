from final import *
import ruptures as rpt
from scipy import integrate
grad = 1000
W = np.linspace(0, 5, 150)
t = np.linspace(0, 120, grad)
A = []
for w in W:
    y = xppaut_model(t, False, nacTOsetp=w)['Int'][6]
    d = np.gradient(y)
    algo = rpt.Dynp(model="l1", min_size=1, jump=.5).fit(d)
    result = algo.predict(n_bkps=1)[0]
    A.append(result)
A = [a*(120/grad) for a in A]

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(W, A, color='black')
ax.set_xlabel("Striatum -> Setp Weight")
ax.set_ylabel("Front-Loading Duration")
plt.show()