import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="Arial", size=16.0)

ns = 10
tf = 0.25

solution = np.loadtxt("solution.csv")

x = solution[0, :]
u = solution[1:, :]

fig2 = plt.figure("An example", figsize=(8, 6))
ax3 = fig2.add_subplot(1, 1, 1)

for i in range(ns-1):
    ax3.plot(x, u[i, :], lw=1, label=f"t = {tf*(i+1)/ns}")

ax3.set_xlabel("$x$")
ax3.set_ylabel("$u$")
ax3.set_xlim(0, 1)
ax3.legend(fontsize=14, loc=0, bbox_to_anchor=(0.3, 0.45, 0.5, 0.5))

fig2.tight_layout()
fig2.savefig("burgers_cds.pdf")


# plt.show()
