/*
  This code transforms the particles position to concentration field
  defined in a square lattice. This field is passed to the code
  HydroGrid.

  The concentration in a cell of volume (dx*dy)  is defined as 
  c = number_of_particles_inside_cell / (dx*dy);

  Therefore the average concentration in the system is
  c_avg = total_number_of_particles / (lx*ly)

  where (lx*ly) is the (2D-) volume of the system.
*/

#include <stdlib.h> 
#include <sstream>
#include <iostream>
#include <fstream>
#include "visit_writer.h"
#include "visit_writer.c"
#include "calculateConcentration.h"

#ifndef NOPYTHON
   #include <boost/python.hpp>
   #include <boost/python/numpy.hpp>
   namespace bp = boost::python;
   namespace np = boost::python::numpy;
#endif
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

void calculateConcentration(string outputname,
                            double lx,       // Domain x length
                            double ly,       // Domain y length
                            int green_start, // Start of "green" particles
                            int green_end,   // End of "green" particles
                            int mx,          // Grid size x
                            int my,          // Grid size y
                            int step,        // Step of simulation
                            double dt,       // Time interval between successive snapshots (calls to updateHydroGrid)
                            int np,          // Number of particles
                            int option,      // option = 0 (initialize), 1 (update), 2 (save), 3 (save+finalize), 4 (finalize only)
                            double *x_array, double *y_array){
  
  static int n_calls = 0;
  static double dx, dy, inverse_volume_cell;
  static double *c, *density, *velocity, *x_old, *y_old;

  dx = lx / mx;
  dy = ly / my;
  inverse_volume_cell = 1.0 / (dx * dy);
  
  c = new double [mx*my*2];
  density = new double [mx*my];
  velocity = new double [mx*my*2]; // vx and vy must be consecutive for each grid point
      
  if(option == 0) { // Initialize HydroGrid
    
    // In order to compute displacements we will need to 
    x_old = new double [np];
    y_old = new double [np];
    n_calls=0;
        
    callHydroGrid(0,
                  outputname,
                  c,
                  density,
                  velocity,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);  
  }
  else if(option == 1){ // Update hydrogrid data
    
    // Set concentration to zero
    for(int i=0; i < mx*my; i++){
      c[i] = 0;  
      c[mx*my+i] = 0;
      density[i] = 0;
      velocity[i] = 0;  
      velocity[mx*my+i] = 0;     
    }
  
    // Loop over particles and save as concentration
    for(int i=0;i<np;i++) {
      // Extract data
      
      double x = x_array[i];     
      double y = y_array[i];
      double v_x = 0;
      double v_y = 0;   

      if(n_calls > 0)
      {
         // We scale the Brownian "velocities" here assuming diffusive dynamics, not balistic:
         v_x = (x-x_old[i])/sqrt(dt);
         v_y = (y-y_old[i])/sqrt(dt);
      }
      else
      {
         v_x = 0;
         v_y = 0;
      }
      x_old[i]=x;
      y_old[i]=y;
      
      // Use PBC
      x = x - (int(x / lx + 0.5*((x>0)-(x<0)))) * lx;
      y = y - (int(y / ly + 0.5*((y>0)-(y<0)))) * ly;
	  
      // Find cell
      int jx   = int(x / dx + 0.5*mx) % mx;
      int jy   = int(y / dy + 0.5*my) % my;
      int icel = jx + jy * mx;

      // Is particle green or red
      if((i>=green_start)&&(i<=green_end)) { // Particle is green
        c[icel] += 1.0;
      }
      else{ // Particle is red
        c[mx*my+icel] += 1.0;
      }
      density[icel] += 1.0;
      velocity[icel] += v_x;
      velocity[mx*my+icel] += v_y;
      
    }

    // Scale concentration and density fields (but not velocities)
    for(int i=0; i < mx*my; i++){
      //if(density[i]>0.0) // Convert to "velocity"
      if(0) // Don't convert to "velocity", use "momentum"
      {
         velocity[i] = velocity[i]/density[i];
         velocity[mx*my+i] = velocity[mx*my+i]/density[i];
      }
      if(1) // Normalize by cell volume (not sure why we need to but it seems to be the right thing to do?)
      {
         velocity[i] = inverse_volume_cell*velocity[i];
         velocity[mx*my+i] = inverse_volume_cell*velocity[mx*my+i];
      }
      density[i] = inverse_volume_cell*density[i];
      c[i]       = inverse_volume_cell*c[i];
      c[mx*my+i] = inverse_volume_cell*c[mx*my+i];
    }

    if(n_calls == 1)
    {
      // Because we didn't have displacements the first time we called this code we reset now
      callHydroGrid(4, // Reset
                  outputname,
                  c,
                  density,
                  velocity,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step); 
    }

    // Call HydroGrid to update data but now add displacements
    n_calls++;
    callHydroGrid(1, // Update
               outputname,
               c,
               density,
               velocity,
               mx,
               my,
               lx,
               ly,
               dt,
               step); 
                  
                  
  }
  else if(option == 2){ // Call HydroGrid to print data
    callHydroGrid(3,
                  outputname,
                  c,
                  density,
                  velocity,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);  
  }  
  else if(option == 3){ // Call HydroGrid to print final data and free memory
    // Free HydroGrid
    callHydroGrid(2,
                  outputname,
                  c,
                  density,
                  velocity,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);
   delete[] x_old;
   delete[] y_old;                  
  }
  else if(option == 4){ // Call HydroGrid to print final data and free memory
    // Free HydroGrid
    callHydroGrid(5,
                  outputname,
                  c,
                  density,
                  velocity,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);
   delete[] x_old;
   delete[] y_old;                  
  }

  // Free memory
  delete[] c;
  delete[] density;
  delete[] velocity;

}

extern "C" { // Interoperable with C and Fortran
   void calculateConcentration_C(char *filename,
           double lx, // Domain x length
           double ly, // Domain y length
           int green_start, // Start of "green" particles
           int green_end, // End of "green" particles
           int mx, // Grid size x
           int my, // Grid size y
           int step, // Step of simulation
           double dt, // Time interval between successive snapshots (calls to updateHydroGrid)
           int np, // Number of particles
           int option, // option = 0 (initialize), 1 (update), 2 (save), 3 (save+finalize), 4 (finalize only)
           double *x_array, double *y_array){
    string outputname(filename); // Convert to C++ string     
    calculateConcentration(outputname,
           lx, ly, green_start, green_end, mx, my, step,
           dt, np, option, x_array, y_array);                            
                            
   }
}

// Python binding, if desired
//-----------------------------------------------------------
#ifndef NOPYTHON 

/*
  Wrapper to call calculateConcentration from python
 */
void calculateConcentrationPython(string outputname,
                                  double lx, 
                                  double ly, 
                                  int green_start,
                                  int green_end,
                                  int mx,
                                  int my,
                                  int step,
                                  double dt,
                                  int np,
                                  int option, // option = 0 (initialize), 1 (update), 2 (save), 3 (finalize)
                                  /*bp::numeric::array r_vectors*/
                                  np::ndarray r_vectors_np) // Should be an Nx3 numpy array
{

  int size;
  int dims = r_vectors_np.get_nd();
  if(dims == 1){
    size = r_vectors_np.shape(0);
  }
  else{
    size = r_vectors_np.shape(0) * r_vectors_np.shape(1);
  }
  int stride = size / np;
  double *r_vectors = reinterpret_cast<double *>(r_vectors_np.get_data());
  
  double* x = new double [np];
  double* y = new double [np];
  
  for(int i=0;i<np;i++){
    // Extract data
    x[i] = r_vectors[stride*i + 0];
    y[i] = r_vectors[stride*i + 1];
  }    
  
  calculateConcentration(outputname, lx, ly, 
                         green_start, green_end, mx, my, step,
                         dt, np, option, x, y);  
  delete[] x;
  delete[] y;                           
}

BOOST_PYTHON_MODULE(libCallHydroGrid)
{
  using namespace boost::python;

  // Initialize numpy
  Py_Initialize();
  np::initialize();
  def("calculate_concentration", calculateConcentrationPython);
}

#endif /* PYTHON */
//-----------------------------------------------------------
