#include <iostream>
#include "PhaseSpace.h"
void split_evolution(vector<PhaseSpace>& ev_pulses, double time);
void modeling(double modelMatrix[modelingXRange][modelingYRange]);
void pause();

//void analyzer(vector<vector<PhaseSpace>> spaces);
vector<vector<PhaseSpace>> analyzer(vector<vector<PhaseSpace>> spaces);

void specModeling(vector<double> &pixelArray);

void write_to_file(double modelMatrix[modelingXRange][modelingYRange]); //QUESTION: Why does the range HAVE to be defined for both the function in the h and cpp file?
void read_from_file(string file);