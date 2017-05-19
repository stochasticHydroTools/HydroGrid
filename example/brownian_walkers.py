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

sys.path.append('../src')

import calculateConcentration as cc



if __name__ == '__main__':

  # Set variables
  num_particles = 2086
  step = 1000
  L = np.array([256.0, 256.0, 0.0])
  r_vectors = np.random.rand(num_particles,3) * L
  cells = np.array([64, 64, 0], dtype=int)
  ly_green = np.array([-L[1] / 6.0, L[1] / 6.0])
  sample = 0
  # Set Gaussian standard deviation along x, y and z
  dx = np.array([1.0, 1.0, 0.0])

  print '#NUMBER PARTICLES ', num_particles
  start = time.time()
  for i in range(step):
    # Random Gaussian displacement
    r_vectors += np.random.randn(r_vectors.shape[0], 3) * dx

    # Update HydroGrid data
    cc.calculate_concentration("run", L[0], L[1], ly_green[0], ly_green[1], int(cells[0]), int(cells[1]), i, 1.0, num_particles, 0, r_vectors)

    if sample > 0:
      if i % sample == 0:
        # Print HydroGrid data
        cc.calculate_concentration("run", L[0], L[1], ly_green[0], ly_green[1], int(cells[0]), int(cells[1]), i, 1.0, num_particles, 1, r_vectors)

  # Print final data and free memory
  cc.calculate_concentration("run", L[0], L[1], ly_green[0], ly_green[1], int(cells[0]), int(cells[1]), i, 1.0, num_particles, 2, r_vectors)


  print 'total time = ', time.time() - start
