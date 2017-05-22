# HydroGrid library for computing static and dynamic structure factors (spectra)


## Usage

First try to run the python example of diffusion in a binary mixture of uncorrelated Brownian walkers, which is in the folder/directory example.

* cd example
* cp ../src/MakefileHeader .
* Edit MakefileHeader to edit compilers/paths/etc
* make
* Edit hydroGridOptions.nml where all options for HydroGrid are set, and then run
* python brownian_walkers.py

In this example we have 
  num_particles = 1024*16
diffusing as independent Brownian walkers in the plane, with diffusion coefficient 1 in a domain of size 256^2 covered by 128^2 grid cells. The first half the particles is "green" (species 1) and the rest "red" (species 2). In the initial configuration all green particles are in the top half, y>Ly/2, and the red in the bottom half, y<=Ly/2, and they mix diffusively.

Every time step the example pyhon script calls the HydroGrid library to compute structure factors (spectra). Rigth now the code outputs all spectra at the end of the run, but if you change sample to a positive value, say 1000, it will output results every 1000 steps. Each time HydroGrid outputs statistics it resets all internal variables and starts anew. Since our configurations evolves dynamically, this is useful to output average statistics over a certain period of time where the statistics don't change much.

==================================
libHydroGrid
==================================

The main input file that controls libHydroGrid is the Fortran namelist file:
hydroGridOptions.nml
This is basically a text input file with ! for comments.
It contains two namelists, deliminated between

&hydroAnalysisOptions
   ...
/

The first namelist corresponds to doing analysis on the full 2D grid. When initialized it prints:

 Initializing HydroAnalysis for grid of size          128         128           1  nSpecies/isSingleFluid=           3 T
           6  Primitive variables indexes as: rho=           1  v=           2 -           3  T=           4  c=           5 -           6  s=           7 -           6

This code will also, however, project the grid along the y axes to get averages that are now just a function of x. When initialized this prints:

 Initializing HydroAnalysis for grid of size          128           1           1  nSpecies/isSingleFluid=           3 T
           6  Primitive variables indexes as: rho=           1  v=           2 -           3  T=           4  c=           5 -        
Spectra are then computed also for this 1D grid. Once can separately control the 2D and 1D analysis via the two separate namelists. I have selected the relevant options for you. The output files for the 2D grid have a prefix "run" set by:
   filePrefix = "run" ! Prefix for all file names
and the output files for the projected 1D grid have file names with "run_proj" since the second namelists sets:
   filePrefix = "run_proj" ! Prefix for all file names

HydroGrid numbers the different variables and can compute spectra of one variable or cross-correlations between variables. The numbering for the present 2D setup is:

! 1=rho (total density), 2-4 (unused = vx, vy, temperature), 5=rho_green, 6=rho_red

where rho here means number density of particles. So for us relevant variables are:
1=rho, 5=rho_green and 6=rho_red

==================================
Static Structure Factors
==================================

When computing structure factors, one can select as many *pairs* of variables as desired. For example, for the 2D grid we set:

   nStructureFactors = 1
   structureFactorPairs = "0 1 0 1" ! The zeros have to be here but mean nothing to us

which tells HydroGrid to compute the rho-rho correlation in Fourier space. The output file
run.S_k.pair=1.Re.dat
contains the static structure factor S(k) for the first pair of variables, i.e., here for rho-rho, real part. You can just use xmgrace to plot it. Everything is normalized including x and y axes. This spectrum should be flat even for true-2D hydro if the total density of particles is uniform.

For the projected 1D grid, we also compute separately spectra for green, red, and also correlations red-green:
   nStructureFactors = 4
   structureFactorPairs = "0 1 0 1, 0 5 0 5, 0 6 0 6, 0 5 0 6"
so the file 
run_proj.S_k.pair=2.Re.dat
contains the spectrum for the green particles (variable pair is (5,5)).
This is the one that is supposed to show giant fluctuations, especially in true 2D-hydro but also in quasi-2D hydro. Right now for uncorrelated walkers are spectra are flat.

==================================
Dynamic Structure Factors
==================================

HydroGrid can also compute dynamic structure factors, i.e., S(k,t) or S(k,w). It does that now but the documentation is not yet finished...

UNFINISHED here on


--------------------------

 nVariablesToFFT=           1  variablesToFFT:            0           1
 grid%structureFactorPairs=           1           1
 Min k =          -63         -63           0
 Max k =           64          64           0
 Recording dynamic structure factors with time interval Dt=   65.536000000000001    
 Keeping track of S(k,omega) and S(k,t) for           45  wavenumbers
...

 Tracking wavenumber index=           7  k_index=           2           0           0

====

    6  s=           7 -           6
 nVariablesToFFT=           3  variablesToFFT:            0           1           0           5           0           6
 grid%structureFactorPairs=           1           1           2           2           3           3           2           3
 Min k =          -63           0           0
 Max k =           64           0           0
 Recording dynamic structure factors with time interval Dt=   65.536000000000001    
 Keeping track of S(k,omega) and S(k,t) for           16  wavenumbers
...

 Tracking wavenumber index=           2  k_index=           2           0           0

=========

xmgrace run_proj.S_k_t.k={1,2,3,4,5}.Re.dat -log y -world 0 1e-3 2000 1

xmgrace run.S_k_t.k=7.Re.dat  run_proj.S_k_t.k=2.Re.dat # Same data

