print("----------------------------------------------------------------------")
print("---------------------------- MAXWELL TEST ----------------------------")
print("----------------------------------------------------------------------\n")



print("--- [ LOADING MODULES ] ----------------------------------------------")

print(" - Matplotlib", end=" ")

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

print("[ OK ]")

print(" - NumPy", end=" ")

import numpy as np
from itertools import combinations

print("[ OK ]")

print(" - Maxwell", end=" ")

mxll_path = "../build/install/lib/python"
import sys

if not (mxll_path in sys.path):
	sys.path.append(mxll_path)

import _maxwell as mxll

print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")



print("--- [ LOADING PARAMETERS ] -------------------------------------------")

print(" - Physical constants", end=" ")

epsilon0 = 8.854187e-12
mu0 = np.pi*4.e-7

print("[ OK ]")

print(" - Time-stepping parameters", end=" ")

tMax = 1.e-8
nb_steps = 10
theta = 0.5
accuracy = 0.01
max_nb_iterations = 2

print("[ OK ]")

print(" - Mesh parameters", end=" ")

Lx = 1.
Ly = 1.
Lz = 1.

nx = 3
ny = 3
nz = 3

print("[ OK ]")

print(" - Initial state", end=" ")

rho0_bg = 0.
j0_bg = mxll.Vec3D(0., 0., 0.)
E0_bg = mxll.Vec3D(0., 0., 0.)
B0_bg = mxll.Vec3D(0., 0., 0.)

wire_sklt = [mxll.Vec3D(0., 0.05, 0.025), mxll.Vec3D(1., 0.05, 0.025)]
wire_rad = 0.005
wire_current = 1.

print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")



print("--- [ COMUTATING ] ---------------------------------------------------")

print(" - Creation of simulation environment", end=" ")
S = mxll.SM()
print("[ OK ]")

print(" - Definition of the constants", end=" ")
S.setConstants(epsilon0, mu0)
print("[ OK ]")

print(" - Creation of the mesh and the time samples", end=" ")
S.setSimulationParameters(tMax, nb_steps, theta, accuracy, max_nb_iterations, 0., Lx, 0., Ly, 0., Lz, nx, ny, nz)
print("[ OK ]")

print(" - Definition of the initial state", end=" ")
S.defineInitialStateBackgroundValues(rho0_bg, j0_bg, E0_bg, B0_bg)
print("[ OK ]")

#S.addWire(wire_sklt, wire_rad, wire_current)
#S.addBoundaryCondition(...)

S.simulate()
print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")



M = np.asarray(S.getMesh())
E = [np.asarray(S.getEnergyDensity(i)) for i in range(nb_steps+1)]



#%% MESH PLOT

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



#%% ENERGY PLOT

print("Plotting energy density field...")

def update_energy_plot(frame_number, mesh, energy_field, ax):
	ax.cla()
	ax.scatter(mesh[:,0], mesh[:,1], mesh[:,2], c=energy_field[frame_number], cmap='jet')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

frames = range(nb_steps+1)
"""
if (nb_steps+1 > max_frames):
    frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]
"""
animate = anim.FuncAnimation(fig, update_energy_plot, frames, fargs=(M, E, ax))
fig.colorbar(cm.ScalarMappable(norm=None, cmap='jet'), ax=ax)
plt.show()
