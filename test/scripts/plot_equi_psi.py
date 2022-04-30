import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(0.968251, 2.40975, 600)
y = np.linspace(-1.4425992, 1.4425992, 600)


ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

levels = np.linspace(0, 0.4, 8)

psi = np.loadtxt('Psi.dat')
ax1.contour(x, y, psi, levels=levels)
ax1.set_title('studentenherz/FOCUS')

psi = np.loadtxt('eqdsk_Psi.dat')
ax2.contour(x, y, psi, levels=levels)
ax2.set_title('cclauser/FOCUS')

ax1.set_aspect('equal')
ax2.set_aspect('equal')
plt.show()