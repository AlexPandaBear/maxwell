print("----------------------------------------------------------------------")
print("---------------------------- MAXWELL TEST ----------------------------")
print("----------------------------------------------------------------------\n")




print("--- [ LOADING MODULES ] ----------------------------------------------")

print("Matplotlib", end=" ")

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

print("[ OK ]")

print("NumPy", end=" ")

import numpy as np
from itertools import combinations

print("[ OK ]")

print("Maxwell", end=" ")

mxll_path = "../build/install/lib/python"
import sys

if not (mxll_path in sys.path):
	sys.path.append(mxll_path)

import _maxwell as mxll

print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ LOADING PARAMETERS ] -------------------------------------------")

print("Physical constants", end=" ")

epsilon0 = 8.854187e-12
mu0 = np.pi*4.e-7

print("[ OK ]")

print("Time-stepping parameters", end=" ")

tMax = 1.e-9
nb_steps = 10
theta = 0
accuracy = 0.00001
max_nb_iterations = 100

print("[ OK ]")

print("Mesh parameters", end=" ")

Lx = 1.
Ly = 1.
Lz = 1.

nx = 4
ny = 4
nz = 4

print("[ OK ]")

print("Initial state", end=" ")

rho0_bg = 0.
j0_bg = mxll.Vec3D(0., 0., 0.)
E0_bg = mxll.Vec3D(0., 0., 0.)
B0_bg = mxll.Vec3D(0., 0., 0.)

wire_sklt = [mxll.Vec3D(0., 0.05, 0.025), mxll.Vec3D(1., 0.05, 0.025)]
wire_rad = 0.005
wire_current = 1.

print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ COMPUTING ] ----------------------------------------------------")

S = mxll.SM()
S.setConstants(epsilon0, mu0, 1.)
S.setSimulationParameters(tMax, nb_steps, theta, accuracy, max_nb_iterations, 0., Lx, 0., Ly, 0., Lz, nx, ny, nz)
S.defineInitialStateBackgroundValues(rho0_bg, j0_bg, E0_bg, B0_bg)
S.lockAllBoundaryNodes()

#S.addWire(wire_sklt, wire_rad, wire_current)
#S.addBoundaryCondition(...)

S.simulate()

print("--------------------------------------------------------- [ DONE ] ---\n")




M = np.asarray(S.getMesh())
T = np.asarray(S.getTime())

R = [np.asarray(S.getRho(i)) for i in range(nb_steps+1)]
J = [np.asarray(S.getJ(i)) for i in range(nb_steps+1)]
E = [np.asarray(S.getE(i)) for i in range(nb_steps+1)]
B = [np.asarray(S.getB(i)) for i in range(nb_steps+1)]

U = [np.asarray(S.getEnergyDensity(i)) for i in range(nb_steps+1)]
P = [np.asarray(S.getPoyntingVector(i)) for i in range(nb_steps+1)]



print("Plotting mesh...")

def update_mesh_plot(frame_number, mesh, ax):
	ids = S.getNodeIds(frame_number)

	for a,b in combinations([0, 1, 2, 3], 2):
		lines.append(ax.plot3D([M[ids[a],0], M[ids[b],0]], [M[ids[a],1], M[ids[b],1]], [M[ids[a],2], M[ids[b],2]], color='k'))

	plt.pause(0.5)

	for i in range(6):
		l = ax.lines.pop(-1)
		del l

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(S.getNbNodes()):
	ax.scatter(M[i,0], M[i,1], M[i,2], color='k')

lines = []
frames = range(S.getNbCells())
animate = anim.FuncAnimation(fig, update_mesh_plot, frames, fargs=(M, ax))

plt.show()



frames = range(nb_steps+1)



print("Plotting charge field...")

def update_charge_plot(frame, mesh, rho, j, ax):
	ax.cla()
	ax.set_title("Charge field plot at step {} (t = {})".format(frame, T[frame]))
	ax.scatter(mesh[:,0], mesh[:,1], mesh[:,2], c=rho[frame], cmap='jet')
	ax.quiver(mesh[:,0], mesh[:,1], mesh[:,2], j[frame][:,0], j[frame][:,1], j[frame][:,2], length=0.1, normalize=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
animate = anim.FuncAnimation(fig, update_charge_plot, frames, fargs=(M, R, J, ax))
#fig.colorbar(cm.ScalarMappable(norm=None, cmap='jet'), ax=ax)
plt.show()



print("Plotting electromagnetic field...")

def update_electromagnetic_plot(frame, mesh, E, B, ax):
	ax.cla()
	ax.set_title("Electromagnetic field plot at step {} (t = {})".format(frame, T[frame]))
	ax.quiver(mesh[:,0], mesh[:,1], mesh[:,2], E[frame][:,0], E[frame][:,1], E[frame][:,2], length=0.1, normalize=True)
	ax.quiver(mesh[:,0], mesh[:,1], mesh[:,2], B[frame][:,0], B[frame][:,1], B[frame][:,2], length=0.1, normalize=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
animate = anim.FuncAnimation(fig, update_electromagnetic_plot, frames, fargs=(M, E, B, ax))
#fig.colorbar(cm.ScalarMappable(norm=None, cmap='jet'), ax=ax)
plt.show()



print("Plotting energy density field...")

def update_energy_plot(frame, mesh, energy, poynting, ax):
	ax.cla()
	ax.set_title("Energy field plot at step {} (t = {})".format(frame, T[frame]))
	ax.scatter(mesh[:,0], mesh[:,1], mesh[:,2], c=energy[frame], cmap='jet')
	ax.quiver(mesh[:,0], mesh[:,1], mesh[:,2], poynting[frame][:,0], poynting[frame][:,1], poynting[frame][:,2], length=0.1, normalize=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
animate = anim.FuncAnimation(fig, update_energy_plot, frames, fargs=(M, U, P, ax))
#fig.colorbar(cm.ScalarMappable(norm=None, cmap='jet'), ax=ax)
plt.show()




print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")