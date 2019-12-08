#include "PhaseSpace.h"
#include "statistics.h"
#include <functional>


//Conservation Checking - Emittence based
bool longitudinal_area_conservation(PhaseSpace ps) {//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = ps.getHWidth() * ps.getVzDist();
	double consValue2 = ps.getHHeight() * ps.getZDist();
	if (consValue > 0.000499 && consValue < 0.000501 && consValue2 > 0.000499 && consValue2 < 0.000501) {//randomly decided range to account for data error
		return true;
	}
	else {
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;;
		return false;
	}
}

double check_energy_conservation(vector<vector<PhaseSpace>> shattered) {
	double totalEnergy = 0;
	for (int i = 0; i < shattered.size(); ++i) {
		for (int j = 0; j < shattered[i].size(); ++j) {
			auto grid = new double[modelingXRange][modelingYRange]{};
			summing(shattered[i][j], grid);
			for (int y = 0; y < modelingYRange; ++y) {
				for (int x = 0; x < modelingXRange; ++x) {
					totalEnergy += grid[x][y];
				}
			}
		}
	}
	return totalEnergy;
}

tuple<double, double, double, double, double, double> valid_variables_check(PhaseSpace ps) {
	tuple<double, double, double, double, double, double> response = make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	double hWidth_ = sqrt(1 / ((1 / (ps.getZDist() * ps.getZDist())) - (ps.getChirp() / ps.getVzDist()) * (ps.getChirp() / ps.getVzDist())));
	double hHeight_ = sqrt(1 / ((1 / pow(ps.getVzDist(), 2)) - pow(ps.getB() / ps.getZDist(), 2)));
	double zDist_ = sqrt(1 / ((1 / pow(ps.getHWidth(), 2)) + pow(ps.getChirp() / ps.getVzDist(), 2)));
	double VzDist_ = sqrt(1 / ((1 / pow(ps.getHHeight(), 2)) + pow((ps.getB() / ps.getZDist()), 2)));
	double chirp_ = ps.getB() * pow(ps.getVzDist() / ps.getZDist(), 2);
	double b_ = ps.getChirp() * pow(ps.getZDist() / ps.getVzDist(), 2);
	cout << "Checks" << endl;
	if (ps.getHWidth() / hWidth_ < 1.01 && ps.getHWidth() / hWidth_ > 0.99) get<0>(response) = 1;
	else	get<0>(response) = ps.getHWidth() / hWidth_;

	if (ps.getHHeight() / hHeight_ < 1.01 && ps.getHHeight() / hHeight_ > 0.99) get<1>(response) = 1;
	else	get<1>(response) = ps.getHHeight() / hHeight_;

	if (ps.getZDist() / zDist_ < 1.01 && ps.getZDist() / zDist_ > 0.99) get<2>(response) = 1;
	else	get<2>(response) = ps.getZDist() / zDist_;

	if (ps.getVzDist() / VzDist_ < 1.01 && ps.getVzDist() / VzDist_ > 0.99) get<3>(response) = 1;
	else	get<3>(response) = ps.getVzDist() / VzDist_;

	if (ps.getChirp() / chirp_ < 1.01 && ps.getChirp() / chirp_ > 0.99) get<4>(response) = 1;
	else	get<4>(response) = ps.getChirp() / chirp_;

	if (ps.getB() / b_ < 1.01 && ps.getB() / b_ > 0.99) get<5>(response) = 1;
	else	get<5>(response) = ps.getB() / b_;

	return response;
}

void print(PhaseSpace ps) {
	cout << "hWidth: " << ps.getHWidth() << endl;
	cout << "hHeight: " << ps.getHHeight() << endl;
	cout << "VzIntDist: " << ps.getVzDist() << endl;
	cout << "zIntDist: " << ps.getZDist() << endl;
	cout << "chirp: " << ps.getChirp() << endl;
	cout << "b: " << ps.getB() << endl;
	cout << "VzC: " << ps.getVzC() << "   " << "zC: " << ps.getZC() << endl;
	cout << "xC: " << ps.getXC() << endl;
	cout << "intensityMultiplier: " << ps.getIntensityMultiplier() << endl;
	//cout << "totalintensityMultiplier from PixelSum: " << sumUp());
	cout << endl;
	cout << "Longitudinal Emmittence Conserved: " << longitudinal_area_conservation(ps) << endl;
	cout << endl;
}

void summing(vector<PhaseSpace> spaces, double grid[modelingXRange][modelingYRange]) {
	double ySearchLB = -3.0 * splitNumber * spaces[0].getHHeight() - (-1 * spaces[spaces.size() - 1].getVzC());
	double ySearchUB = 3.0 * splitNumber * spaces[0].getHHeight() + (-1 * spaces[spaces.size() - 1].getVzC());
	double xSearchLB = -3.5 * spaces[0].getHWidth();
	double xSearchUB = 3.5 * spaces[0].getHWidth();

	double yPulseUB = 3.5 * spaces[0].getHHeight();
	double yPulseLB = -3.5 * spaces[0].getHHeight();
	double xPulseUB = 3.5 * spaces[0].getHWidth();
	double xPulseLB = -3.5 * spaces[0].getHWidth();
	
	for (PhaseSpace pulse : spaces) {
		//Convert to grid_integration parameters
		double xGridHalfRange = (xSearchUB - xSearchLB) / 2;
		double yGridHalfRange = (ySearchUB - ySearchLB) / 2;

		double xHalfRange = (xPulseUB - xPulseLB) / 2;
		double yHalfRange = (yPulseUB - yPulseLB) / 2;

		double xOffset = (xPulseUB + xPulseLB) / 2;
		double yOffset = (yPulseUB + yPulseLB) / 2;

		pulse.grid_integration(xHalfRange, yHalfRange, xOffset, yOffset, grid, xGridHalfRange, yGridHalfRange);
	}
}

void summing(PhaseSpace space, double grid[modelingXRange][modelingYRange]) {
	double ySearchLB = -3.5 * space.getHHeight();
	double ySearchUB = 3.5 * space.getHHeight();
	double xSearchLB = -3.5 * space.getHWidth();
	double xSearchUB = 3.5 * space.getHWidth();

	//Convert to grid_integration parameters
	double xHalfRange = (xSearchUB - xSearchLB) / 2;
	double yHalfRange = (ySearchUB - ySearchLB) / 2;

	//There are no offsets, phase space summing is centered around 0.
	double xOffset = 0;
	double yOffset = 0;

	space.grid_integration(xHalfRange, yHalfRange, xOffset, yOffset, grid, xHalfRange, yHalfRange);
}

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange], double grid3[modelingXRange][modelingYRange]) {
	for (int i = 0; i < modelingXRange; i++) {
		for (int j = 0; j < modelingYRange; j++) {
			grid3[i][j] = grid1[i][j] - grid2[i][j];
		}
	}
}

double measureDeviation(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]) {
	double deviation = 0.0;
	for (int i = 0; i < modelingXRange; i++) {
		for (int j = 0; j < modelingYRange; j++) {
			deviation += pow(grid1[i][j] - grid2[i][j], 2);
		}
	}
	return sqrt(deviation / (double(modelingXRange) * double(modelingYRange) - 1.0));
}

double measureDeviation(vector<double> base, vector<double> compare) {//Compares two pixel sum vectors
	if (base.size() != compare.size()) {
		return -1;//Error code for debugging
	}
	int size = base.size();
	double deviation = 0;
	for (int i = 0; i < size; ++i) {
		deviation += pow(base[i] - compare[i], 2);
	}
	return sqrt(deviation / size);
}

vector<vector<PhaseSpace>> analyzer(vector<vector<PhaseSpace>> spaces) {
	//for (vector<PhaseSpace> space : spaces) {
	//	for (PhaseSpace pulse : space) {
	//		pulse.spectroscopy_function();
	//	}
	//}

	/*for (vector<PhaseSpace> &spaceSet : spaces) {
		for (PhaseSpace &space: spaceSet){
			space.spectroscopy_function();
		}
	}*/
	vector<vector<PhaseSpace>> returnSpaces;
	for (int i = 0; i < spaces.size(); i++) {
		vector<PhaseSpace> spaceSet;
		for (int j = 0; j < spaces[0].size(); j++) {
			spaceSet.push_back(spaces[i][j].spectroscopy_function());
		}
		returnSpaces.push_back(spaceSet);
	}
	return returnSpaces;
}

vector<PhaseSpace> analyzer(vector<PhaseSpace> spaces) {
	vector<PhaseSpace> returnSpaces;
	for (int i = 0; i < spaces.size(); ++i) {
		returnSpaces.push_back(spaces[i].spectroscopy_function());
	}
	return returnSpaces;
}

void pixelSum(vector<double>& pixelArray, vector<vector<PhaseSpace>> spaces) {
	double lowestXC = spaces[spaces.size() - 1][spaces[0].size() - 1].getXC();
	double highestXC = spaces[0][0].getXC();
	int i = 0;
	for (int i = 0; i < spaces.size(); i++) {
		for (int j = 0; j < spaces[0].size(); j++) {
			pixelArray[int(map(spaces[i][j].getXC(), lowestXC, highestXC, 0, double(pixels) - 1) + 0.5)] += spaces[i][j].getXC() * spaces[i][j].getIntensityMultiplier();
		}
	}
}

void pixelSum(vector<double>& pixelArray, vector<PhaseSpace> spaces) {
	double lowestXC = spaces[spaces.size() - 1].getXC();
	double highestXC = spaces[0].getXC();
	double counter = 0;
	for (int i = 0; i < spaces.size(); i++) {//Iterate over pixels
		//for (int j = 0; j < spaces.size(); ++j) {//Iterate over spaces
			//if(-3*i*spaces[0].getVxDist() < )
			pixelArray[int(map(spaces[i].getXC(), lowestXC - 0.01, highestXC, 0, double(pixels) - 1) + 0.5)] += spaces[i].getIntensityMultiplier();
			//counter += spaces[i].getXC();// *spaces[i].getIntensityMultiplier();
		//}
	}
	//cout << counter << endl;
}

void pixelSum(vector<double>& pixelArray, double highest, double lowest, vector<vector<double>> v) {
	for (int i = 0; i < v.size(); i++) {
		pixelArray[int(map(v[i][0], lowest, highest, 0, double(pixels) - 1) + 0.5)] += v[i][1];
	}
}


/*
void print() {
	cout << "hWidth: " << HWidth << "  " << " hDepth: " << HDepth << endl;;
	cout << "hHeight: " << hHeight << " " << " hDepthVelocity: " << hDepthVel << endl;
	cout << "VzIntDist: " << VzDist << " VxIntDist: " << VxDist << endl;
	cout << "zIntDist: " << zDist << " " << " xIntDist: " << xDist << endl;
	cout << "chirp: " << chirp << " " << " chirpT: " << chirpT << endl;
	cout << "b: " << b << "     " << " bT: " << bT << endl;
	cout << "VzC: " << VzC << "   " << "zC: " << zC << endl;
	cout << "xC: " << xC << endl;
	cout << "pulseEnergy: " << pulseEnergy << endl;
	cout << "intensityMultiplier: " << intensityMultiplier << endl;
	//cout << "totalintensityMultiplier from PixelSum: " << sumUp());
	cout << endl;
	cout << "Longitudinal Emmittence Conserved: " << longitudinal_area_conservation(hWidth, VzDist, hHeight, zDist) << endl;
	cout << "Transverse Emmittence Conserved: " << transverse_area_conservation(hDepth, VxDist, hDepthVel, xDist) << endl;
	cout << endl;
}*/



double map(double input, double inMin, double inMax, double outMin, double outMax) {
	if (input < inMin) {
		input = inMin;
		cout << "under" << endl;
	}
	else if (input > inMax) {
		input = inMax;
		cout << "over" << endl;
	}
	return (input - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
}


/*
bool transverse_area_conservation(PhaseSpace ps) {//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = ps.getHDepth() * ps.getVxDist();
	double consValue2 = ps.getHDepthVel() * ps.getXDist();
	if (consValue > 0.00299 && consValue < 0.00301 && consValue2 > 0.00299 && consValue2 < 0.00301) {//randomly decided range to account for data error
		return true;
	}
	else {
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;
		return false;
	}
}*/