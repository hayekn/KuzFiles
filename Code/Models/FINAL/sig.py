import numpy as np
import matplotlib.pyplot as plt

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))
noise = np.random.normal(0,.3, 150)

fig, ax = plt.subplots(figsize=(5, 5))
x = np.linspace(-10, 10, 150)
ax.plot(x, F(x+noise), color='black', label='Exc')
ax.plot(x, F(-x+noise), '--', color='black', label='Inh')
ax.set_xlabel("Incoming Projections")
ax.set_ylabel("Derivative of Firing Rate")
plt.legend()
plt.savefig("Code/Models/FINAL/Pics/sig.png", transparent=True, dpi=350)