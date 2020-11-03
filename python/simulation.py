import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as anim
import _maxwell as mxll



epsilon0 = 8.854187e-12
mu0 = np.pi*4.e-7



tMax = 1.e-8
nb_steps = 10
theta = 0.5
accuracy = 0.01
max_nb_iterations = 10



Lx = 1.
Ly = 1.
Lz = 1.

nx = 2
ny = 2
nz = 2



V0 = mxll.Vec3D(0., 0., 0.)

rho0_bg = 0.
j0_bg = V0
E0_bg = V0
B0_bg = V0



wire_sklt = [mxll.Vec3D(0., 0.05, 0.025), mxll.Vec3D(1., 0.05, 0.025)]
wire_rad = 0.005
wire_current = 1.


S = mxll.SM()
S.setConstants(epsilon0, mu0)
S.setSimulationParameters(tMax, nb_steps, theta, accuracy, max_nb_iterations, 0., Lx, 0., Ly, 0., Lz, nx, ny, nz)
print("yolo")
S.defineInitialStateBackgroundValues(rho0_bg, j0_bg, E0_bg, B0_bg)
#S.addWire(wire_sklt, wire_rad, wire_current)
#S.addBoundaryCondition(...)
S.simulate()

Mesh = S.getMesh()
Energy = S.getEnergyDensity(0)

print(type(S))
print(type(Mesh))
print(type(Energy))

M = np.asarray(Mesh)
E = np.asarray(Energy)

print(len(M))
print(len(E))


#%% MESH DISPLAY

print("Plotting mesh...")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(M[:,0])):
	ax.scatter(M[i,0], M[i,1], M[i,2], color='k')

plt.show()




#%% ENERGY PLOT

print("Plotting energy density field...")

def update_energy_plot(frame_number, mesh, energy_field, ax):
	ax.cla()
	for i in range(len(mesh)):
		c = energy_field[i]
		ax.scatter(mesh[i,0], mesh[i,1], mesh[i,2], color='k')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

frames = range(nb_steps+1)
"""
if (nb_steps+1 > max_frames):
    frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]
"""
animate = anim.FuncAnimation(fig, update_energy_plot, frames, fargs=(M, E, ax))
plt.show()
