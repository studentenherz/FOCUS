import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('d3dpart_pol.dat')
t, r, theta, z, vr, vtheta, vz, Br, Btheta, Bz = np.transpose(data)

fig = plt.figure()
# ax = plt.axes(projection='3d', title='Trajectory')

# ax.plot3D(x, y, z)
# plt.show()

v0 = 0.6 * 1.84142e7

ax = plt.axes(title='Energy')
ax.plot(t, (vr**2 + vtheta**2 + vz**2)/ v0**2)
plt.show()

ax = plt.axes(title='$p_\\theta$')
ax.plot(t, vtheta / v0)
plt.show()

B = np.sqrt(Br**2 + Btheta**2 + Bz**2)

vtr = (vtheta * Bz - vz * Btheta) / B
vttheta = (vz * Br - vr * Bz) / B
vtz = (vr * Btheta - vtheta * Br) / B

vt2 = vtr**2 + vttheta**2 + vtz**2

ax = plt.axes(title='$\mu$')
ax.plot(t, (vt2 / B) / (vt2[0] / B[0]))
plt.show()

