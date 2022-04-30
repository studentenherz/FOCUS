import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('xys.dat')

x = np.linspace(0.968251, 2.40975, 600)
y = np.linspace(-1.4425992, 1.4425992, 600)


ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

u = np.loadtxt('Br.dat')
v = np.loadtxt('Bz.dat')
ax1.streamplot(x, y, u, v, density=1.4, color='g')
ax1.set_title('studentenherz/FOCUS')

u = np.loadtxt('eqdsk_Br.dat')
v = np.loadtxt('eqdsk_Bz.dat')
ax2.streamplot(x, y, u, v, density=1.4)
ax2.set_title('cclauser/FOCUS')

ax1.set_aspect('equal')
ax2.set_aspect('equal')
plt.show()