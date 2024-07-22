from final import *

plt.rcParams.update({
        "text.usetex": True,
        'font.size': 28
    })

t = np.linspace(0, 120, 800)
nac = xppaut_model(t, False)['Int'][3]
dls = xppaut_model(t, False)['Int'][4]

result = 100*(dls/(dls+nac))

fig, ax = plt.subplots(1, 2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]}, figsize=(10, 9))
for i in range(2):
    ax[i].plot(t, result, label="DLS", linewidth=4, color='black')
    ax[i].plot(t, 100-result, marker='v', markersize=15, markevery=15, label="NAc", linewidth=4, color = 'black')


ax[0].set_ylim(-2, 100)
ax[0].set_xlim(-1, 40)
ax[1].set_xlim(115, 120)

n=26
ax[1].set_xticks([120], ['120'], fontsize=n)
ax[0].set_yticks(np.arange(0, 101, 20), list(map(str, np.arange(0, 101, 20))), fontsize=n)
ax[0].set_xticks([0, 10, 20, 30, 40], ['0', '10', '20', '30', '40'], fontsize=n)

ax[0].spines['right'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].tick_params(left = False) 
d = .015  
r = 0
s = 0.15
kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
ax[0].plot((1-d+r, 1+d+r), (-d, +d), **kwargs)
ax[0].plot((1-d+r, 1+d+r), (1-d, 1+d), **kwargs)
kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
ax[1].plot((-d-0.3+s, +d+s), (1-d, 1+d), **kwargs)
ax[1].plot((-d-0.3+s, +d+s), (-d, +d), **kwargs)
ax[0].set_xlabel('Time (min)')
ax[0].set_ylabel('Striatum Load (\%)')
ax[0].xaxis.set_label_coords(0.6, -0.09)

ax[0].fill_between(t, 100, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label='CS') 
ax[1].fill_between(t, 100, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label='CS') 
plt.legend()
plt.savefig("Code/Models/FINAL/Pics/dlsFig.png", transparent=True, dpi=400)


