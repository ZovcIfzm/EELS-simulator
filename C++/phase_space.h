#include <iostream>
#include <vector>
#include "external_variables.h"
#include "constants.h"
using namespace std;


class phase_space{
private:
	double hWidth, hHeight, VzDist, zDist, chirp, b, pulseEnergy, intensityMultiplier,
	hDepth, hDepthVel, VxDist, xDist, chirpT, bT, vZC, zC, xC;
	//double spectroTable[][];
	//double modelingSpace[400][400];
	double totalEnergy, electronAmount;

public:
	phase_space();
	phase_space(double testConstructor);
	phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityMultiplierC,
		double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC);
	//phase_space(phase_space[splitNumber]);
	phase_space(vector<phase_space> spaces);
	vector<phase_space> split();
	phase_space evolution(double time);
	phase_space RFLens(double changeChirp);
	phase_space mag_lens(double changeChirpT);
	phase_space spectroscopy_function();
	void valid_variables_check();
	double get_split_intensity_multiplier(double numSections, double sectionNum, double hHeight, double hWidth);
	bool longitudinal_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
	bool transverse_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
	void print();

	//DEBUGGING FUNCTions
	double hWidth_accessor();
};