from final import *

plt.rcParams.update({
        "text.usetex": True,
        'font.size': 22
    })

t = np.linspace(0, 120, 300)
nac = xppaut_model(t, False)['Int'][3]
dls = xppaut_model(t, False)['Int'][4]

result = dls/(dls+nac)
fig, ax = plt.subplots(figsize=(8, 8))
ax.plot(t, result, label="DLS")

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

ax.plot(t, 1-result, label="NAc")
ax.fill_between(t, 1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label='CS') 
ax.set_xlabel("T (min)")
plt.legend()
plt.savefig("Code/Models/FINAL/Pics/dlsFig.png", transparent=False, dpi=400)


