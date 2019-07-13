#include <iostream>
#include "PhaseSpace.h"
double map(double input, double inMin, double inMax, double outMin, double outMax);
void split_evolution(vector<PhaseSpace>& ev_pulses, double time);
void summing(vector<PhaseSpace> spaces, double graphingGrid[modelingXRange][modelingYRange]);
void summing(PhaseSpace space, double graphingGrid[modelingXRange][modelingYRange]);
void modeling(double modelMatrix[modelingXRange][modelingYRange]);
void pause();
void pixelSum(vector<double> &pixelArray, vector<vector<PhaseSpace>> spaces);
void specModeling(vector<double> &pixelArray);

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange], double grid3[modelingXRange][modelingYRange]);
double measureDeviation(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]);

void write_to_file(double modelMatrix[modelingXRange][modelingYRange]); //QUESTION: Why does the range HAVE to be defined for both the function in the h and cpp file?
void read_from_file(string file);