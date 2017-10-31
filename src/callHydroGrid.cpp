
#include <stdlib.h> 
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;


extern "C" {
#include "HydroGrid.h"
}

// Meaning of option=0: Initialize, option=1: Add new data (update), 
// option=2: Finalize, option=3: Save data and finalize,
// option=4: Just reset without saving, option=5: Just finalize without saving data
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
                   const int step){

  // Define variables
  static int nCells[3];
  static double systemLength[3];

  if(option == 0){
    ifstream fileinput ("hydroGridOptions.nml");
    string word, wordfile;
    while(!fileinput.eof()){
      getline(fileinput,word);
      wordfile += word + "\n";
    }
    fileinput.close();
    string fileOutName = outputname + ".hydroGridOptions.nml";
    ofstream fileout(fileOutName.c_str());
    fileout << wordfile << endl;
    fileout.close();

    nCells[0] = mx;
    nCells[1] = my;
    nCells[2] = 1;
    systemLength[0] = lx;
    systemLength[1] = ly;
    systemLength[2] = 1;   // 0
    createHydroAnalysis_C(nCells,
                          3 /*nSpecies*/,
                          2 /*nVelocityDimensions*/,
                          1 /*isSingleFluid*/,
                          systemLength,
                          NULL /*heatCapacity*/,
                          dt /*time step*/,
                          0 /*nPassiveScalars*/,
                          1 /*structFactMultiplier*/,
                          1 /*project2D*/);
  }
  else if(option == 1){
    updateHydroAnalysisMixture_C(velocity /*velocities*/, density /*densities*/, c /*concentrations*/);
  }
  else if(option == 2){
    writeToFiles_C(-1); // Write to files
    destroyHydroAnalysis_C();
  }
  else if(option == 3){
    writeToFiles_C(step); // Write to files
  }
  else if(option == 4){
    resetHydroAnalysis_C(); // Write to files
  }
  else if(option == 5){
    destroyHydroAnalysis_C();
  }
  return 0;
}
