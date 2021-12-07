#include <iostream>
#include "PhaseSpace.h"
void split_evolution(vector<PhaseSpace>& ev_pulses, double time);
void modeling(double modelMatrix[modelingXRange][modelingYRange]);
void pause();

void deviationModeling(vector<pair<double, double>>& pixelArray, double xMin, double xMax, string title, string xAxis, string yAxis);
void specModeling(vector<double> &pixelArray, double xMin, double xMax);
void energyModeling(vector<double>& pixelArray);

void write_to_file(double modelMatrix[modelingXRange][modelingYRange]); //QUESTION: Why does the range HAVE to be defined for both the function in the h and cpp file?
void read_from_file(string file);