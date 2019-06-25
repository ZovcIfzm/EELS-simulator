#include "modeling.h"
#include "string"
#include <fstream>
#include <iomanip>
#include "gnuplot_i.hpp"
using namespace std;

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

void split_evolution(vector<PhaseSpace>& ev_pulses, double time) {
	for (PhaseSpace &pulse : ev_pulses) {
		pulse.evolution(time);
	}
}

void modeling(double modelMatrix[modelingXRange][modelingYRange]) {
	if (printStarts)
		cout << "modeling started" << endl;
	vector <double> x, y, z;
	for (int row = 0; row < modelingYRange; row++) {
		for (int col = 0; col < modelingXRange; col++) {
			x.push_back(col);
			y.push_back(row);
			z.push_back(modelMatrix[col][row]);
		}
	}
	Gnuplot modeling("model");
	//modeling.set_grid();
	//modeling.set_samples(10);
	//modeling.set_contour("base");
	modeling.cmd("set dgrid3d 50,50");
	modeling.cmd("set style data lines");
	modeling.cmd("set pm3d");
	modeling.cmd("set palette rgb 23,28,3");
	write_to_file(modelMatrix);
	modeling.cmd("splot 'modeling_data.txt' pal");
	//modeling.set_hidden3d();
	//modeling.plot_xyz(x, y, z, "user-defined points 3d");
	pause();
	if(printEnds)
		cout << "modeling finished" << endl;
	modeling.cmd("exit");
}

void summing(vector<PhaseSpace> spaces, double grid[modelingXRange][modelingYRange]) {
	if (printStarts){
		cout << "summing multiple phase spaces started" << endl;
	}
	double ySearchLB = -3.0 * splitNumber * spaces[0].getHHeight() - (-1*spaces[spaces.size()-1].getVzC());
	double ySearchUB = 3.0 * splitNumber * spaces[0].getHHeight() + (-1*spaces[spaces.size()-1].getVzC());
	double xSearchLB = -3.5 * spaces[0].getHWidth();
	double xSearchUB = 3.5 * spaces[0].getHWidth();

	double yPulseUB = 3.5*spaces[0].getHHeight();
	double yPulseLB = -3.5*spaces[0].getHHeight();
	double xPulseUB = 3.5*spaces[0].getHWidth();
	double xPulseLB = -3.5*spaces[0].getHWidth();

	double accuracyY = (yPulseUB - yPulseLB) / 199;
	double accuracyX = (xPulseUB - xPulseLB) / 199;

	double negTwohWidthsq = -2 * spaces[0].getHWidth()*spaces[0].getHWidth();
	double twoVzIntDistsq = 2 * spaces[0].getVzDist()*spaces[0].getVzDist();
	double twoPIhWidthVzIntDist = 2 * M_PI*(spaces[0].getHWidth()*spaces[0].getVzDist());
	
	double x = xPulseLB + accuracyX / 2;
	double y = yPulseLB + accuracyY / 2; 
	for (PhaseSpace pulse : spaces) {
		x = xPulseLB + accuracyX / 2;
		y = yPulseLB + accuracyY / 2;
		while (y < yPulseUB) {
			while (x < xPulseUB) {
				if (x + pulse.getZC() > xSearchLB && x + pulse.getZC() < xSearchUB && y + pulse.getVzC() > ySearchLB && y + pulse.getVzC() < ySearchUB) {
					grid[int(map(x + pulse.getZC(), xSearchLB, xSearchUB, 0, modelingXRange - 1)+0.5)][int(map(y + pulse.getVzC(), ySearchLB, ySearchUB, 0, modelingYRange - 1)+0.5)] += pulse.getIntensityMultiplier()*(accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].getChirp() * x)*(y - spaces[0].getChirp() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
					valueHolder2 += pulse.getIntensityMultiplier() * (accuracyX * accuracyY * exp((x * x / (negTwohWidthsq)) - ((y - spaces[0].getChirp() * x) * (y - spaces[0].getChirp() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
				}
				x += accuracyX;
			}
			y += accuracyY;
			x = xPulseLB + accuracyX / 2;
		}
	}
	if (printEnds) {
		cout << "summing multiple phase spaces finished" << endl;
	}
}

void summing(PhaseSpace space, double grid[modelingXRange][modelingYRange]) {
	if (printStarts) {
		cout << "summing single phase space started" << endl;
	}
	double ySearchLB = -3.5 * space.getHHeight();
	double ySearchUB = 3.5 *space.getHHeight();
	double xSearchLB = -3.5 *space.getHWidth();
	double xSearchUB = 3.5 *space.getHWidth();
	double accuracyY = (ySearchUB - ySearchLB) / 199;  // Cannot be 99 or else for some reason the middle (25th row out of 50) takes half that of all the other rows. Look into?
	double accuracyX = (xSearchUB - xSearchLB) / 199;
	double negTwohWidthsq = -2 * space.getHWidth()*space.getHWidth();
	double twoVzIntDistsq = 2 * space.getVzDist()*space.getVzDist();
	double twoPIhWidthVzIntDist = 2 * M_PI*(space.getHWidth()*space.getVzDist());
	double x = xSearchLB + accuracyX / 2;
	double y = ySearchLB + accuracyY / 2;
	while (y < ySearchUB) {
		while (x < xSearchUB) {
			if (x + space.getZC() > xSearchLB && x + space.getZC() < xSearchUB && y + space.getVzC() > ySearchLB && y + space.getVzC() < ySearchUB) {
				grid[int(map(x, xSearchLB, xSearchUB, 0, modelingXRange - 1) + 0.5)][int(map(y, ySearchLB, ySearchUB, 0, modelingYRange - 1) + 0.5)] += (accuracyX*accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - space.getChirp() * x)*(y - space.getChirp() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
				//valueHolder += accuracyX*accuracyY*(1 / (2 * M_PI))*exp((x*x + y * y) / (-2));
				valueHolder1 += (accuracyX*accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - space.getChirp() * x)*(y - space.getChirp() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
			}
			x += accuracyX;
		}
		y += accuracyY;
		x = xSearchLB + accuracyX / 2;
	}
	if (printEnds) {
		cout << "summing single phase space finished" << endl;
	}
}

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange], double grid3[modelingXRange][modelingYRange]) {
	 for (int i = 0; i < modelingXRange; i++){
		 for (int j = 0; j < modelingYRange; j++) {
			grid3[i][j] = grid1[i][j] - grid2[i][j];
		 }
	 }
}

double measureDeviation(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]) {
	double deviation = 0.0;
	for (int i = 0; i < modelingXRange; i++) {
		for (int j = 0; j < modelingYRange; j++) {
			deviation += pow(grid1[i][j] - grid2[i][j],2);
		}
	}
	return sqrt(deviation/(double(modelingXRange)*double(modelingYRange) - 1.0));
}


void write_to_file(double grid[modelingXRange][modelingYRange]) {//Untested
	if (printStarts)
		cout << "writing to file started" << endl;
	ofstream modelFile;
	modelFile.open("modeling_data.txt");
	for (int row = 0; row < modelingYRange; row++) {
		for (int col = 0; col < modelingXRange; col++) {
			modelFile << col << "\t" << row << "\t" << grid[col][row] << endl;
		}
	}
	modelFile.close();
	if (printEnds)
		cout << "writing to file finished" << endl;
}

void read_from_file(string file) {
	if (printStarts)
		cout << "reading file started" << endl;
	ifstream modelFileRead(file);
	int x, y;
	double z;
	while (modelFileRead >> x >> y >> z) {
		graphingMap[x - 1][y - 1] = z;
	}
	if (printEnds)
		cout << "reading file finished" << endl;
}

void pause() {
	string random;
	cout << "Press enter to continue...";
	getline(cin, random);
}