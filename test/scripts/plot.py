import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('d3dpart2_pol.dat')
t, r, theta, z, vr, vtheta, vz, Br, Btheta, Bz, Psi, Fpol = np.transpose(data)

x = r * np.cos(theta)
y = r * np.sin(theta)

fig = plt.figure()
ax = plt.axes(projection='3d', title='Trajectory')
ax.plot3D(x, y, z)
plt.show()

ax = plt.axes()

x = np.linspace(1.0116, 2.36566, 600)
y = np.linspace( -1.363, 1.348, 600)

u = np.loadtxt('Br.dat')
v = np.loadtxt('Bz.dat')
ax.streamplot(x, y, u, v, density=1.4, color='#ccc')

bdrydata = np.loadtxt('bdry.dat')
xlimit = bdrydata[:,0]
ylimit = bdrydata[:,1]
ax.plot(xlimit, ylimit, 'r-')

limdata = np.loadtxt('lim.dat')
xlimit = limdata[:,0]
ylimit = limdata[:,1]
ax.plot(xlimit, ylimit, 'k-')

ax.plot(r, z)

ax.set_aspect('equal')
plt.show()

v0 = (0.05**2 + 0.1**2)**0.5 * 1.84142e7
q_over_m = 9.58e7



E = (vr**2 + vtheta**2 + vz**2)
Em, dE = np.mean(E), np.std(E)
ax = plt.axes(title='Energy')
ax.plot(t, E)
ax.plot(t, [Em for _ in range(len(t))])
plt.show()

Ptheta = r * vtheta + q_over_m * Psi
Pmean, dP = np.mean(Ptheta), np.std(Ptheta)
ax = plt.axes(title=f'$p_\\theta = {Pmean:.3f} \pm {dP:.3f}$')
ax.plot(t, Ptheta)
ax.plot(t, [Pmean for _ in range(len(t))])
plt.show()

B = np.sqrt(Br**2 + Btheta**2 + Bz**2)

vtr = (vtheta * Bz - vz * Btheta) / B
vttheta = (vz * Br - vr * Bz) / B
vtz = (vr * Btheta - vtheta * Br) / B

vt2 = vtr**2 + vttheta**2 + vtz**2
mu = vt2 / B

mum, dmu = np.mean(mu), np.std(mu)
ax = plt.axes(title=f'$\mu = {mum:.3e} \pm {dmu:.3e}$')
ax.plot(t, mu)
ax.plot(t, [mum for _ in range(len(t))])
plt.show()