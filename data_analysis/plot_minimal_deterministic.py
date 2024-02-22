import numpy as np
# from itertools import permutations
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
c = 2.5e-3
Nss = 10

def f(Y, t):
    w, m = Y
    rate = c * (Nss-w-m)
    return [rate*w, rate*m]

w = np.arange(0, 16, 1)
m = np.arange(0, 16, 1)

W, M = np.meshgrid(w, m)

rate = c * (Nss-W-M)
dW = rate * W
dM = rate * M

fig, ax = plt.subplots(figsize = (4, 4))
ax.quiver(W, M, dW, dM)
ax.plot([0, Nss], [Nss, 0], linestyle = "-.", color = "black", label = "$w + m = N_{ss}$")

ts = np.linspace(0, 300, 501)
for w0, m0 in [[0,1], [1,0], [1,1], [1,2], [2,1],
               [0,15], [15,0], [15, 15], [7.5, 15], [15, 7.5]]:
        sols = odeint(f, [w0, m0], ts)
        ax.plot(sols[:,0], sols[:,1], color = "blue")
ax.arrow(0, 5, 0, 1, head_width = 0.25, color = "blue")
ax.arrow(5, 0, 1, 0, head_width = 0.25, color = "blue")
ax.arrow(2, 4, 1, 2, head_width = 0.25, color = "blue")
ax.arrow(4, 2, 2, 1, head_width = 0.25, color = "blue")
ax.arrow(3, 3, 1, 1, head_width = 0.25, color = "blue")
ax.arrow(12, 0, -1, 0, head_width = 0.25, color = "blue")
ax.arrow(0, 12, 0, -1, head_width = 0.25, color = "blue")
ax.arrow(5, 10, -1, -2, head_width = 0.25, color = "blue")
ax.arrow(10, 5, -2, -1, head_width = 0.25, color = "blue")
ax.arrow(8, 8, -1, -1, head_width = 0.25, color = "blue")
ax.legend()
ax.set_xlabel("w")
ax.set_ylabel("m")
plt.tight_layout()
plt.savefig("../figures/minimal_ode_quiver.png")
plt.show()
