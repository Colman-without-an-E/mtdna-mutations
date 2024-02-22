import numpy as np
# from itertools import permutations
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
c = 2.5e-3
Nss = 10
kra = 2

def f(Y, t):
    w, m = Y
    rate = c * (Nss-w-m)
    return [c*(Nss-w-m)*w, c*(Nss-w-m+kra)*m]

w = np.arange(0, 16, 1)
m = np.arange(0, 16, 1)

W, M = np.meshgrid(w, m)

dW = c * (Nss-W-M) * W
dM = c * (Nss-W-M+kra) * M

fig, ax = plt.subplots(figsize = (4, 4))
ax.quiver(W, M, dW, dM)
ax.plot([0, Nss], [Nss+kra, 0], linestyle = "-.", color = "black", label = "$w + m = N_{ss} + k_{ra}$")

ts = np.linspace(0, 1000, 5001)
for w0, m0 in [[0,1], [1,0], [1,1], [1,2], [2,1],
               [0,15], [15,0], [15, 15], [7.5, 15], [15, 7.5]]:
        sols = odeint(f, [w0, m0], ts)
        ax.plot(sols[:,0], sols[:,1], color = "blue")
ax.arrow(0, 8, 0, 1, head_width = 0.25, color = "blue")
ax.arrow(8, 0, 1, 0, head_width = 0.25, color = "blue")
ax.arrow(3, 8.4, -1, 1.23, head_width = 0.25, color = "blue")
ax.arrow(5, 6, -1, 1.25, head_width = 0.25, color = "blue")
# ax.arrow(2, 4, 1, 2, head_width = 0.25, color = "blue")
# ax.arrow(4, 2, 2, 1, head_width = 0.25, color = "blue")
# ax.arrow(3, 3, 1, 1, head_width = 0.25, color = "blue")
ax.arrow(14, 0, -1, 0, head_width = 0.25, color = "blue")
ax.arrow(0, 14, 0, -1, head_width = 0.25, color = "blue")
# ax.arrow(5, 10, -1, -2, head_width = 0.25, color = "blue")
# ax.arrow(10, 5, -2, -1, head_width = 0.25, color = "blue")
# ax.arrow(8, 8, -1, -1, head_width = 0.25, color = "blue")
ax.legend()
ax.set_xlabel("w")
ax.set_ylabel("m")
plt.tight_layout()
# plt.savefig("../figures/ra_ode_quiver.png")
plt.show()

