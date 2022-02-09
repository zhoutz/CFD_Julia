import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="Arial", size=16.0)

ns = 10
tf = 0.25

solution = np.loadtxt("solution_d.csv")
xd = solution[0, :]
ud = solution[1:, :]

solution = np.loadtxt("solution_p.csv")
xp = solution[0, :]
up = solution[1:, :]

fig = plt.figure("An example", figsize=(14, 6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for i in range(ns):
    ax1.plot(xd, ud[i, :], lw=1, label=f"t = {tf*(i+1)/ns}")

ax1.set_xlabel("$x$")
ax1.set_ylabel("$u$")
ax1.set_title("Dirichlet boundary")
ax1.set_xlim(0, 1)
ax1.legend(fontsize=14, loc=0, bbox_to_anchor=(0.55, 0.35, 0.5, 0.5))

for i in range(ns):
    ax2.plot(xp, up[i, :], lw=1, label=f"t = {tf*(i+1)/ns}")

ax2.set_xlabel("$x$")
ax2.set_ylabel("$u$")
ax2.set_title("Periodic boundary")
ax2.set_xlim(0, 1)
ax2.legend(fontsize=14, loc=0, bbox_to_anchor=(0.55, 0.35, 0.5, 0.5))


fig.tight_layout()
fig.savefig("crweno_2.pdf")


# plt.show()
