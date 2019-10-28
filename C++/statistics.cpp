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

void phase_space_integration(double grid[modelingXRange][modelingYRange], PhaseSpace ps, double xHalfRange, double xAccuracy, double xOffset, double yHalfRange, double yAccuracy, double yOffset) {
	double negTwohWidthsq = -2 * ps.getHWidth() * ps.getHWidth();
	double twoVzIntDistsq = 2 * ps.getVzDist() * ps.getVzDist();
	double twoPIhWidthVzIntDist = 2 * M_PI * (ps.getHWidth() * ps.getVzDist());
	double x = -xHalfRange + xAccuracy / 2;
	double y = -yHalfRange + yAccuracy / 2;
	while (y < yHalfRange) {
		while (x < xHalfRange) {
			if (x + ps.getZC() > -xHalfRange && x + ps.getZC() < xHalfRange && y + ps.getVzC() > yHalfRange && y + ps.getVzC() < yHalfRange) {
				grid[int(map(x, -xHalfRange, xHalfRange, 0, modelingXRange - 1) + 0.5)][int(map(y, -yHalfRange, yHalfRange, 0, modelingYRange - 1) + 0.5)] += (xAccuracy * yAccuracy * ps.intensity(x,y));
				valueHolder1 += (xAccuracy * yAccuracy * ps.intensity(x,y));
			}
			x += xAccuracy;
		}
		y += yAccuracy;
		x = -xHalfRange + xAccuracy / 2;
	}
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

	double accuracyY = (yPulseUB - yPulseLB) / 199;
	double accuracyX = (xPulseUB - xPulseLB) / 199;

	double negTwohWidthsq = -2 * spaces[0].getHWidth() * spaces[0].getHWidth();
	double twoVzIntDistsq = 2 * spaces[0].getVzDist() * spaces[0].getVzDist();
	double twoPIhWidthVzIntDist = 2 * M_PI * (spaces[0].getHWidth() * spaces[0].getVzDist());

	double x = xPulseLB + accuracyX / 2;
	double y = yPulseLB + accuracyY / 2;
	for (PhaseSpace pulse : spaces) {
		x = xPulseLB + accuracyX / 2;
		y = yPulseLB + accuracyY / 2;
		while (y < yPulseUB) {
			while (x < xPulseUB) {
				if (x + pulse.getZC() > xSearchLB && x + pulse.getZC() < xSearchUB && y + pulse.getVzC() > ySearchLB && y + pulse.getVzC() < ySearchUB) {
					
					valueHolder2 += pulse.getIntensityMultiplier() * (accuracyX * accuracyY * pulse.intensity(x, y));
				}
				x += accuracyX;
			}
			y += accuracyY;
			x = xPulseLB + accuracyX / 2;
		}
	}
}

void summing(PhaseSpace space, double grid[modelingXRange][modelingYRange]) {
	double ySearchLB = -3.5 * space.getHHeight();
	double ySearchUB = 3.5 * space.getHHeight();
	double xSearchLB = -3.5 * space.getHWidth();
	double xSearchUB = 3.5 * space.getHWidth();
	double accuracyY = (ySearchUB - ySearchLB) / 199;  // Cannot be 99 or else for some reason the middle (25th row out of 50) takes half that of all the other rows. Look into?
	double accuracyX = (xSearchUB - xSearchLB) / 199;
	double negTwohWidthsq = -2 * space.getHWidth() * space.getHWidth();
	double twoVzIntDistsq = 2 * space.getVzDist() * space.getVzDist();
	double twoPIhWidthVzIntDist = 2 * M_PI * (space.getHWidth() * space.getVzDist());
	double x = xSearchLB + accuracyX / 2;
	double y = ySearchLB + accuracyY / 2;
	while (y < ySearchUB) {
		while (x < xSearchUB) {
			if (x + space.getZC() > xSearchLB && x + space.getZC() < xSearchUB && y + space.getVzC() > ySearchLB && y + space.getVzC() < ySearchUB) {
				grid[int(map(x, xSearchLB, xSearchUB, 0, modelingXRange - 1) + 0.5)][int(map(y, ySearchLB, ySearchUB, 0, modelingYRange - 1) + 0.5)] += (accuracyX * accuracyY * space.intensity(x,y));
				//valueHolder += accuracyX*accuracyY*(1 / (2 * M_PI))*exp((x*x + y * y) / (-2));
				valueHolder1 += (accuracyX * accuracyY * space.intensity(x,y));
			}
			x += accuracyX;
		}
		y += accuracyY;
		x = xSearchLB + accuracyX / 2;
	}
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

/*void analyzer(vector<vector<PhaseSpace>> spaces) {
	//for (vector<PhaseSpace> space : spaces) {
	//	for (PhaseSpace pulse : space) {
	//		pulse.spectroscopy_function();
	//	}
	//}
	for (int i = 0; i < spaces.size(); i++) {
		for (int j = 0; j < spaces[0].size(); j++) {
			spaces[i][j].spectroscopy_function();
		}
	}
}*/

vector<vector<PhaseSpace>> analyzer(vector<vector<PhaseSpace>> spaces) {
	//for (vector<PhaseSpace> space : spaces) {
	//	for (PhaseSpace pulse : space) {
	//		pulse.spectroscopy_function();
	//	}
	//}
	vector<vector<PhaseSpace>> pulses;

	for (int i = 0; i < spaces.size(); i++) {
		vector<PhaseSpace> a;
		for (int j = 0; j < spaces[0].size(); j++) {
			a.push_back(spaces[i][j].spectroscopy_function());
		}
		pulses.push_back(a);
	}
	return pulses;
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
	for (int i = 0; i < spaces.size(); i++) {
		pixelArray[int(map(spaces[i].getXC(), lowestXC - 0.01, highestXC, 0, double(pixels) - 1) + 0.5)] += spaces[i].getXC() * spaces[i].getIntensityMultiplier();
		counter += spaces[i].getXC();// *spaces[i].getIntensityMultiplier();
	}
	cout << counter << endl;
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