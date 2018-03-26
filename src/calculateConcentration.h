/* Reduced C++ interface to HydroGrid for just working with concentrations/number densities of particles in Quasi2D */
bool callHydroGrid(const int option,
                   const string outputname,
                   double *c,
                   double *density,
                   double *velocity,
                   const int mx,
                   const int my,
                   const double lx,
                   const double ly,
                   const double dt,
                   const int step);

/* This function converts particle data to hydro data and calls callHydroGrid */
void calculateConcentration(std::string outputname,
   double lx, // Domain x length
   double ly, // Domain y length
   int green_start, // Start of "green" particles
   int green_end, // End of "green" particles
   int mx, // Grid size x
   int my, // Grid size y
   int step, // Step of simulation
   double dt, // Time interval between successive snapshots (calls to updateHydroGrid)
   int np, // Number of particles
   int option, // option = 0 (initialize), 1 (update), 2 (save), 3 (save+finalize), 4 (finalize only),
   double *x_array, double *y_array);

