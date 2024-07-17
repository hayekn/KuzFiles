from final import *
t = np.linspace(0, 120, 300)
nac = xppaut_model(t, False)['Int'][3]
dls = xppaut_model(t, False)['Int'][4]

result = dls/(dls+nac)
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(t, result, label="DLS")
ax.plot(t, 1-result, label="NAc")
ax.fill_between(t, 1, where=[(t >= 0) and (t <= 3) for t in t], color = 'grey', alpha = 0.15, linewidth = 0.05, label='CS') 
plt.title("Striatal Contribution To Consumption")
ax.set_xlabel("T (min)")
plt.legend()
plt.show()


