#pragma once
#include <iostream>
#include <vector>
#include "constants.h"
using namespace std;

class PhaseSpace{
private:
	double hWidth, hHeight, VzDist, zDist, chirp, b, pulseEnergy, intensityMultiplier,
	hDepth, hDepthVel, VxDist, xDist, chirpT, bT, VzC, zC, xC;
	//double spectroTable[][];
	//double modelingSpace[400][400];
	double totalEnergy, electronAmount;

public:
	PhaseSpace(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityMultiplierC,
		double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double VzCC, double zCC, double xCC);
	//PhaseSpace(PhaseSpace[splitNumber]);
	PhaseSpace(vector<PhaseSpace> spaces);
	vector<PhaseSpace> split();
	vector<PhaseSpace> PhaseSpace::shatter(vector<vector<double>> data);

	//Evolution methods
	PhaseSpace evolution(double time);
	PhaseSpace RFLens(double changeChirp);
	PhaseSpace mag_lens(double changeChirpT);
	PhaseSpace spectroscopy_function();

	//Debugging methods
	double get_split_intensity_multiplier(double numSections, double sectionNum, double hHeight, double hWidth);

	//Accessor methods
	double getHWidth();
	double getHHeight();
	double getVzDist();
	double getZDist();
	double getChirp();
	double getB();
	double getIntensityMultiplier();
	double getVzC();
	double getXC();
	double getZC();
};