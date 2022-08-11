import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('d3dpart_pol.dat')
t, r, theta, z, vr, vtheta, vz, Br, Btheta, Bz, Psi, Fpol = np.transpose(data)

x = r * np.cos(theta)
y = r * np.sin(theta)

fig = plt.figure()
ax = plt.axes(projection='3d', title='Trajectory')
ax.plot3D(x, y, z)
plt.show()

ax = plt.axes()
ax.plot(r, z)
ax.set_aspect('equal')
plt.show()

v0 = 0.6 * 1.84142e7
q_over_m = 9.58e7

ax = plt.axes(title='Energy')
ax.plot(t, (vr**2 + vtheta**2 + vz**2))
ax.plot(t, [v0**2 for _ in range(len(t))])
plt.show()

Ptheta = r * vtheta + q_over_m * Psi

ax = plt.axes(title=f'$p_\\theta = {np.mean(Ptheta):.3f} \pm {np.std(Ptheta):.3f}$')
ax.plot(t, Ptheta)
plt.show()

B = np.sqrt(Br**2 + Btheta**2 + Bz**2)

vtr = (vtheta * Bz - vz * Btheta) / B
vttheta = (vz * Br - vr * Bz) / B
vtz = (vr * Btheta - vtheta * Br) / B

vt2 = vtr**2 + vttheta**2 + vtz**2

ax = plt.axes(title='$\mu$')
ax.plot(t, (vt2 / B) / (vt2[0] / B[0]))
plt.show()