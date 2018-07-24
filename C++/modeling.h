#include <iostream>
#include "phase_space.h"
double map(double input, double inMin, double inMax, double outMin, double outMax);
void split_evolution(vector<phase_space>& ev_pulses, double time);
void summing(vector<phase_space> spaces);
void summing(phase_space space);
void modeling();
void pause_method();

void writing_to_file(double modelMatrix[modelingXRange][modelingYRange]); //QUESTION: Why does the range HAVE to be defined for both the function in the h and cpp file?
void reading_from_file();