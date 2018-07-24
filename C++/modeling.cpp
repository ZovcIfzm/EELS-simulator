#include "modeling.h"
#include "string"
#include <fstream>
#include <iomanip>
using namespace std;

double map(double input, double inMin, double inMax, double outMin, double outMax) {
	if (input < inMin) {
		input = inMin;
	}
	else if (input > inMax) {
		input = inMax;
	}
	return (input - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
}

void split_evolution(vector<phase_space>& ev_pulses, double time) {
	for (phase_space &pulse : ev_pulses) {
		pulse.evolution(10000);
	}
}

void summing(vector<phase_space> spaces) {
	if (printStarts){
		cout << "started summing multiple phase spaces" << endl;
	}
	double ySearchLB = -2 * spaces[0].hHeight_accessor();
	double ySearchUB = 2 * spaces[0].hHeight_accessor();
	double xSearchLB = -2 * spaces[0].hWidth_accessor();
	double xSearchUB = 2 * spaces[0].hWidth_accessor();
	double accuracyY = (ySearchUB - ySearchLB) / 99;
	double accuracyX = (xSearchUB - xSearchLB) / 99;
	double negTwohWidthsq = -2 * spaces[0].hWidth_accessor()*spaces[0].hWidth_accessor();
	double twoVzIntDistsq = 2 * spaces[0].VzDist_accessor()*spaces[0].VzDist_accessor();
	double twoPIhWidthsqVzIntDistsq = 2 * M_PI*(spaces[0].hWidth_accessor()*spaces[0].hWidth_accessor()*spaces[0].VzDist_accessor()*spaces[0].VzDist_accessor());
	for (phase_space pulse : spaces) {
		double x = xSearchLB + accuracyX / 2;
		double y = ySearchLB + accuracyY / 2;
		while (y < ySearchUB) {
			while (x < xSearchUB) {
				graphingMap[int(map(x+pulse.xC_accessor(), xSearchLB, xSearchUB, 0, modelingXRange - 1))][int(map(y+pulse.VzC_accessor(), ySearchLB, ySearchUB, 0, modelingYRange - 1))] += pulse.intensity_multiplier_accessor()*(accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthsqVzIntDistsq)) / 2000;
				//cout << map(x + pulse.xC_accessor(), xSearchLB, xSearchUB, 0, modelingXRange - 1) << " next: " << pulse.intensity_multiplier_accessor()*(accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthsqVzIntDistsq)) / 2000 << endl;
				x += accuracyX;
			}
			y += accuracyY;
			x = xSearchLB + accuracyX / 2;
		}
	}
	if (printEnds) {
		cout << "ended summing multiple phase spaces" << endl;
	}
}

void summing(phase_space space) {
	if (printStarts) {
		cout << "started summing single phase space" << endl;
	}
	double ySearchLB = -2 * space.hHeight_accessor();
	double ySearchUB = 2 * space.hHeight_accessor();
	double xSearchLB = -2 * space.hWidth_accessor();
	double xSearchUB = 2 * space.hWidth_accessor();
	double accuracyY = (ySearchUB - ySearchLB) / 99;
	double accuracyX = (xSearchUB - xSearchLB) / 99;
	double negTwohWidthsq = -2 * space.hWidth_accessor()*space.hWidth_accessor();
	double twoVzIntDistsq = 2 * space.VzDist_accessor()*space.VzDist_accessor();
	double twoPIhWidthsqVzIntDistsq = 2 * M_PI*(space.hWidth_accessor()*space.hWidth_accessor()*space.VzDist_accessor()*space.VzDist_accessor());
	double x = xSearchLB + accuracyX / 2;
	double y = ySearchLB + accuracyY / 2;
	while (y < ySearchUB) {
		while (x < xSearchUB) {
			graphingMap[int(map(x + space.xC_accessor(), xSearchLB, xSearchUB, 0, modelingXRange - 1))][int(map(y + space.VzC_accessor(), ySearchLB, ySearchUB, 0, modelingYRange - 1))] += (accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - space.chirp_accessor() * x)*(y - space.chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthsqVzIntDistsq));
			x += accuracyX;
		}
		y += accuracyY;
		x = xSearchLB + accuracyX / 2;
	}
	if (printEnds) {
		cout << "ended summing single phase space" << endl;
	}
}

void pause_method() {
	string random;
	cout << "Press enter to continue..." << endl;
	getline(cin, random);
}

void writing_to_file(double modelMatrix[modelingXRange][modelingYRange]) {//Untested
	ofstream modelFile;
	modelFile.open("simulation_model_matrix.txt");
	modelFile << "#X \t Y \t Z" << endl;
	for (int row = 0; row < modelingYRange; row++) {
		for (int col = 0; col < modelingXRange; col++) {
			modelFile << col << "\t" << row << "\t" << modelMatrix[col][row] << endl;
		}
	}
	modelFile.close();
}

void reading_from_file() {
}