import numpy as np
import matplotlib.pyplot as plt
import _maxwell as mxll



epsilon0 = 8.854187e-12
mu0 = np.pi*4.e-7



tMax = 1.e-8
nb_steps = 10
theta = 0.5
accuracy = 0.01
max_nb_iterations = 10



Lx = 1.
Ly = 0.1
Lz = 0.05

nx = 10
ny = 5
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
S.defineInitialStateBackgroundValues(rho0_bg, j0_bg, E0_bg, B0_bg)
S.addWire(wire_sklt, wire_rad, wire_current)
#S.addBoundaryCondition(...)
S.simulate()



M = S.getMesh()
E = S.getEnergyDensity()



#%% ENERGY PLOT

def update_energy_plot(frame_number, mesh, energy_field, ax):
    ax.cla()
    for node in mesh:
    	c = ...
	    ax.scatter(node.get_x(), node.get_y(), node.get_z(), color=c)

print("Plotting energy density field...")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

frames = range(nb_steps+1)
"""
if (nb_steps+1 > max_frames):
    frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]
"""
animate = animation.FuncAnimation(fig, update_energy_plot, frames, fargs=(M, E, ax))
plt.show()
