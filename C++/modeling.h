#include <iostream>
#include "phase_space.h"
double map(double input, double inMin, double inMax, double outMin, double outMax);
void split_evolution(vector<phase_space>& ev_pulses, double time);
void summing(vector<phase_space> spaces, double graphingGrid[modelingXRange][modelingYRange]);
void summing(phase_space space, double graphingGrid[modelingXRange][modelingYRange]);
void modeling(double modelMatrix[modelingXRange][modelingYRange]);
void pause();

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]);
double* t_g_s(double grid1[2][2], double grid2[2][2]);

void write_to_file(double modelMatrix[modelingXRange][modelingYRange]); //QUESTION: Why does the range HAVE to be defined for both the function in the h and cpp file?
void read_from_file(string file);