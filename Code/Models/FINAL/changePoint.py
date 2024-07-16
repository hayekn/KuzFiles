from final import *
import ruptures as rpt

grad = 1000

t = np.linspace(0, 120, grad)
y = xppaut_model(t, False, csTOvta=3)['Int'][6]
d = np.gradient(y)


algo = rpt.Dynp(model="l2", min_size=1, jump=1).fit(d)
result = algo.predict(n_bkps=1)
print(result)
rpt.display(d, result)


A = []
C = np.linspace(1, 3, grad)
for c in C:
    y = xppaut_model(t, False, csTOvta=c)['Int'][6]
    d = np.gradient(y)
    algo = rpt.Dynp(model="l2", min_size=1, jump=.5).fit(d)
    result = algo.predict(n_bkps=1)[0]
    A.append(result)
    print(result)
fig, ax = plt.subplots(figsize=(10, 8))
A = [a*(120/grad) for a in A]
ax.set_ylim(0, 45)
ax.set_ylabel("T (min)")
ax.set_xlabel("Sensitivity of DA to CS")
ax.plot(C, A, label='Change Point')



plt.show()

