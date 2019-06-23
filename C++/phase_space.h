#pragma once
#include <iostream>
#include <vector>
#include "constants.h"
using namespace std;


class phase_space{
private:
	double hWidth, hHeight, VzDist, zDist, chirp, b, pulseEnergy, intensityMultiplier,
	hDepth, hDepthVel, VxDist, xDist, chirpT, bT, VzC, zC, xC;
	//double spectroTable[][];
	//double modelingSpace[400][400];
	double totalEnergy, electronAmount;

public:
	phase_space();
	phase_space(double testConstructor);
	phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityMultiplierC,
		double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double VzCC, double zCC, double xCC);
	//phase_space(phase_space[splitNumber]);
	phase_space(vector<phase_space> spaces);
	vector<phase_space> split();
	vector<phase_space> phase_space::shatter();
	phase_space evolution(double time);
	phase_space RFLens(double changeChirp);
	phase_space mag_lens(double changeChirpT);
	phase_space spectroscopy_function();
	tuple<double, double, double, double, double, double> valid_variables_check();
	double get_split_intensity_multiplier(double numSections, double sectionNum, double hHeight, double hWidth);
	bool longitudinal_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
	bool transverse_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
	void print();

	//Accessor methods
	double hWidth_accessor();
	double hHeight_accessor();
	double VzDist_accessor();
	double zDist_accessor();
	double chirp_accessor();
	double b_accessor();
	double intensity_multiplier_accessor();
	double VzC_accessor();
	double xC_accessor();
	double zC_accessor();
};