/*
  This code transforms the particles position to concentration field
  defined in a square lattice. This field is passed to the code
  HydroGrid.

  The concentration in a cell of volume (dx*dy)  is defined as 
  c = number_of_particles_inside_cell / (dx*dy);

  Therefore the average concentration in the system is
  c_avg = total_number_of_particles / (lx*ly)

  where (lx*ly) is the (2D-) volume of the system.

  HOW TO COMPILE:
  1. Compile code HydroGrid.
  1. Use Makefile in fluam/tools/ to compile this code.

  HOW TO USE:
  1. Edit (if necessary) the input file hydroGridOptions.nml used by
     HydroGrid.
  2. Run the command
     spectra.exe run.particles outputname lx ly mx my nsteps sample [Npout] [Npin] [savefreq]

  with
  run.particles = file with the particles positions generated by fluam
  outputname = prefix for the output file. 
  lx = system's length along the x-axis
  ly = system's length along the y-axis
  mx = number of cells along the x-axis          
  my = number of cells along the y-axis
  nsteps = How many configurations to read from the file
  sample = analyze data only 1 of every "sample" steps.
  Npout = (default all particles) which particle to stop reading at
  Npin = (default 1) which particle to start reading from
  savefreq = (default -1) save data every savefreq. If savefreq=-1
             save data only at the end.
*/

#include <stdlib.h> 
#include <sstream>
#include <iostream>
#include <fstream>
#include "visit_writer.h"
#include "visit_writer.c"

#include <boost/python.hpp>
#include <math.h>
#include <stdio.h>
#include <vector>

namespace bp = boost::python;
using namespace std;

bool callHydroGrid(const int option,
                   const string outputname,
                   double *c,
                   double *density,
                   const int mx,
                   const int my,
                   const double lx,
                   const double ly,
                   const double dt,
                   const int step);

void calculateConcentration(string outputname,
                            double lx, 
                            double ly, 
                            double ly_green_start,
                            double ly_green_end,
                            int mx,
                            int my,
                            int step,
                            double dt,
                            int np,
                            int option,
                            bp::numeric::array r_vectors){
  
  // string outputname = "hola";
  static int init = 0;
  static double dx, dy, inverse_volume_cell;
  static double *y_init, *c, *density;

  if(init == 0){
    // Init function
    init = 1;
    dx = lx / mx;
    dy = ly / my;
    inverse_volume_cell = 1.0 / (dx * dy);
    y_init = new double [np];
    c = new double [mx*my*2];
    density = new double [mx*my];
  
    // Initialize HydroGrid
    callHydroGrid(0,
                  outputname,
                  c,
                  density,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);
    
    // Save initial y-coordinate
    for(int i=0;i<np;i++){
      bp::numeric::array r_vector_1 = bp::extract<bp::numeric::array>(r_vectors[i]);
      y_init[i] = bp::extract<double>(r_vector_1[1]);     
      y_init[i] = y_init[i] - (int(y_init[i] / ly + 0.5*((y_init[i]>0)-(y_init[i]<0)))) * ly;
    }
  }


  if(option == 0){ // Update hydrogrid data
    // Set concentration to zero
    for(int i=0; i < mx*my; i++){
      c[i] = 0;  
      c[mx*my+i] = 0;
      density[i] = 0;
    }
  
    // Loop over particles and save as concentration
    for(int i=0;i<np;i++){
      // Extract data
      bp::numeric::array r_vector_1 = bp::extract<bp::numeric::array>(r_vectors[i]);
      double x = bp::extract<double>(r_vector_1[0]);     
      double y = bp::extract<double>(r_vector_1[1]);     
      
      // Use PBC
      x = x - (int(x / lx + 0.5*((x>0)-(x<0)))) * lx;
      y = y - (int(y / ly + 0.5*((y>0)-(y<0)))) * ly;
	  
      // Find cell
      int jx   = int(x / dx + 0.5*mx) % mx;
      int jy   = int(y / dy + 0.5*my) % my;
      int icel = jx + jy * mx;

      // Is particle green or red
      if((y_init[i] < ly_green_start) or (y_init[i] > ly_green_end)){ // Particle is red
        c[mx*my+icel] += 1.0;
      }
      else{ // Particle is green
        c[icel] += 1.0;
      }
      density[icel] += 1.0;
    }

    // Scale concentration and density fields
    for(int i=0; i < mx*my; i++){
      density[i] = inverse_volume_cell*density[i];
      c[i]       = inverse_volume_cell*c[i];
      c[mx*my+i] = inverse_volume_cell*c[mx*my+i];
    }

    // Call HydroGrid to update data
    callHydroGrid(1,
                  outputname,
                  c,
                  density,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);  
  }
  else if(option == 1){ // Call HydroGrid to print data
    callHydroGrid(3,
                  outputname,
                  c,
                  density,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);  
  }  
  else if(option == 2){ // Call HydroGrid to print final data and free memory
    // Free HydroGrid
    callHydroGrid(2,
                  outputname,
                  c,
                  density,
                  mx,
                  my,
                  lx,
                  ly,
                  dt,
                  step);

    // Free memory
    delete[] c;
    delete[] density;
    delete[] y_init;
  }
}



BOOST_PYTHON_MODULE(calculateConcentration)
{
  using namespace boost::python;
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("calculate_concentration", calculateConcentration);
}
