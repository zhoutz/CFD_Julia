import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="Arial", size=16.0)

final_field = np.loadtxt("field_final.csv", skiprows=1)

x = final_field[:, 0]
u_e = final_field[:, 1]
u_n = final_field[:, 2]
u_error = final_field[:, 3]

u_error = np.abs(u_error)

fig = plt.figure("FTCS", figsize=(14, 6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

ax1.plot(x, u_e, lw=4, ls="-", color="b", label="Exact solution")
ax1.plot(x, u_n, lw=4, ls="--", color="r", label="Crank-Nicolson solution")
ax1.set_xlabel("$x$")
ax1.set_ylabel("$u$")
ax1.set_title("Solution field")
ax1.set_xlim(-1, 1)
ax1.legend(fontsize=14, loc=0)

ax2.plot(x, u_error, marker="o", markeredgecolor="k",
         markersize=8, color="g", lw=4)
ax2.set_xlabel("$x$")
ax2.set_ylabel("$ϵ$")
ax2.set_title("Discretization error")
ax2.set_xlim(-1, 1)

fig.tight_layout()
fig.savefig("cn.pdf")

# plt.show()
