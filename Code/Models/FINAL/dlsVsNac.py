from final import *
t = np.linspace(0, 120, 300)
nac = xppaut_model(t, False)['Int'][3]
dls = xppaut_model(t, False)['Int'][4]

result = dls/(dls+nac)

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(t, result, label="DLS")
ax.plot(t, 1-result, label="NAc")
plt.title("Striatal Contribution To Consumption")
ax.set_xlabel("T (min)")
plt.legend()
plt.show()


