#pragma once
#include "PhaseSpace.h"
void print(PhaseSpace ps);
tuple<double, double, double, double, double, double> valid_variables_check(PhaseSpace ps);
bool longitudinal_area_conservation(PhaseSpace ps);

double check_energy_conservation(vector<vector<PhaseSpace>> shattered);

void summing(vector<PhaseSpace> spaces, double graphingGrid[modelingXRange][modelingYRange]);
void summing(PhaseSpace space, double graphingGrid[modelingXRange][modelingYRange]);

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange], double grid3[modelingXRange][modelingYRange]);
double measureDeviation(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]);
double measureDeviation(vector<double> base, vector<double> compare);//Compares two pixel sum vectors

vector<vector<PhaseSpace>> analyzer(vector<vector<PhaseSpace>> spaces);
vector<PhaseSpace> analyzer(vector<PhaseSpace> spaces);

pair<int, double> pixelSumHelper(vector<double>& pixelArray, vector<vector<PhaseSpace>> spaces, int p, double lowestXC, double xCDist);

void pixelSum(vector<double>& pixelArray, vector<vector<PhaseSpace>> spaces);
void pixelSumMulti(vector<double>& pixelArray, vector<vector<PhaseSpace>> spaces);
void pixelSum(vector<double>& pixelArray, vector<PhaseSpace> spaces);
void pixelSum(vector<double>& pixelArray, double highest, double lowest, vector<vector<double>> v);


double map(double input, double inMin, double inMax, double outMin, double outMax);

//bool transverse_area_conservation(PhaseSpace ps);