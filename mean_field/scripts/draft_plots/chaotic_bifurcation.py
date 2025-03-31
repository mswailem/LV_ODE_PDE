# %%
import numpy as np
from matplotlib import pyplot as plt
# %%
# Importing data
data = np.loadtxt('../../output/bifurcation/draft_chaotic_bifurcation.dat')

alpha = data[:, 0]
a = data[:, 1]
# %%
plt.rcParams.update({'font.size': 20})  # You can change 14 to any size you prefer

plt.figure(figsize=(10, 8))
plt.scatter(alpha, a, c = 'black')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$a_p$', rotation=0, labelpad=20)
plt.savefig('../../plots/draft_chaos_bifurcation/k1_0.24_n_0.7.pdf', dpi=300)
plt.show()
