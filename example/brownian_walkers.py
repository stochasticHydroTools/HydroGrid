'''
Example to compute the 2D concentration of green and red particles
and to call HydroGrid.

The user should have a input file hydroGridOptions.nml
and call the function:

cc.calculate_concentration(outputname, lx, ly, ly_green_start, ly_green_end, mx, my, step, dt, num_particles, option, r_vectors)

with:
outputname = prefix of the output file *.hydroGridOptions.nml
lx = system's length along the x-axis
ly = system's length along the y-axis
ly_green_start = particles with y-coordinate  
                 ly_green_start < y < ly_green_end 
                 in the initial step are green other wise they are red.
ly_green_end = see above
mx = number of cells along the x-axis
my = number of cells along the y-axis
step = simulation step
dt = time step between calls to HydroGrid; in general it will be
     dt = dt_simulation * frequency_analyze_data
num_particles = number of particles
option = 0 analyze data, 1 print data, 2 print final data and free memory
r_vectors = array of dimensions (num_particles, 3) with the particle positions
'''


import numpy as np
import sys
import time

sys.path.append('./')

import calculateConcentration as cc


if __name__ == '__main__':

  # Set variables
  num_particles = 1024*16
  chi = 1 # Diffusion coefficient
  nsteps = 10000 # Number of steps
  sample = 0 # How often to output HydroGrid statistics, if zero only output stats at the end
  L = np.array([256.0, 256.0])
  cells = np.array([128, 128], dtype=int)  

  # Generate phase-separated particle configuration
  last_green = num_particles/2
  r_vectors = np.random.rand(num_particles,2)
  r_vectors[:,0] = r_vectors[:,0]*L[0] # x coordinates
  r_vectors[0:last_green,1] = r_vectors[0:last_green,1] * L[1]/2 # y coordinates in bottom half
  r_vectors[last_green+1:,1] = r_vectors[last_green+1:,1] * L[1]/2 + L[1]/2 # y coordinates in upper half   
  
  # Set Gaussian standard deviation along x, y and z
  dt = L[1]*L[1]/chi/1000.0 # Time step
  dx = np.array([1.0, 1.0])*np.sqrt(2*chi*dt)

  # Initialize HydroGrid library:
  cc.calculate_concentration("run", L[0], L[1], 0, last_green, int(cells[0]), int(cells[1]), 0, dt, num_particles, 0, r_vectors)

  print 'num_particles= ', num_particles, ' dt=', dt
  start = time.time()
  for step in range(nsteps+1):
    print 'step=', step

    # Update HydroGrid data
    cc.calculate_concentration("run", L[0], L[1], 0, last_green, int(cells[0]), int(cells[1]), step, dt, num_particles, 1, r_vectors)

    if sample > 0:
      if step % sample == 0:
        # Print HydroGrid data
        cc.calculate_concentration("run", L[0], L[1], 0, last_green, int(cells[0]), int(cells[1]), step, dt, num_particles, 2, r_vectors)

    # Random Gaussian displacement
    r_vectors += np.random.randn(r_vectors.shape[0], 2) * dx

  # Print final data and free memory
  cc.calculate_concentration("run", L[0], L[1], 0, last_green, int(cells[0]), int(cells[1]), step, dt, num_particles, 3, r_vectors)

  print 'total time = ', time.time() - start
