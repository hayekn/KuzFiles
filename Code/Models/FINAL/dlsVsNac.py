from final import *
t = np.linspace(0, 120, 300)
nac = xppaut_model(t, False)['Int'][3]
dls = xppaut_model(t, False)['Int'][4]

result = dls/(dls+nac)
nacResult
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(t, result)
ax.set_ylabel("DLS Contribution To Consumption")
ax.set_xlabel("T (min)")
plt.show()


