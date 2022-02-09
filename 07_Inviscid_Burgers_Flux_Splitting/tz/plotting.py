import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="Arial", size=16.0)

ns = 11
tf = 0.25

solution = np.loadtxt("solution_flux_split.csv")
xd = solution[0, :]
ud = solution[1:, :]

fig = plt.figure("An example", figsize=(8, 6))
ax1 = fig.add_subplot(1, 1, 1)

for i in range(ns):
    ax1.plot(xd, ud[i, :], lw=1, label=f"t = {tf*(i)/(ns-1)}")

ax1.set_xlabel("$x$")
ax1.set_ylabel("$u$")
ax1.set_title("Periodic boundary - flux splitting")
ax1.set_xlim(0, 1)
ax1.legend(fontsize=14, loc=0, bbox_to_anchor=(0.57, 0.3, 0.5, 0.5))

fig.tight_layout()
fig.savefig("burgers_flux_split.pdf")
