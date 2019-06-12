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

void split_evolution(vector<phase_space>& ev_pulses, double time) {
	for (phase_space &pulse : ev_pulses) {
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
	modeling.cmd("splot 'modeling_data.txt' pal");
	//modeling.set_hidden3d();
	//modeling.plot_xyz(x, y, z, "user-defined points 3d");
	pause();
	if(printEnds)
		cout << "modeling finished" << endl;
	modeling.cmd("exit");
}

void summing(vector<phase_space> spaces, double graphingGrid[modelingXRange][modelingYRange]) {
	if (printStarts){
		cout << "summing multiple phase spaces started" << endl;
	}
	double ySearchLB = -1 * splitNumber * spaces[0].hHeight_accessor() - 2.5*spaces[0].hHeight_accessor();
	double ySearchUB = 1 * splitNumber * spaces[0].hHeight_accessor() + 2.5*spaces[0].hHeight_accessor();
	double xSearchLB = -3.5 * spaces[0].hWidth_accessor();
	double xSearchUB = 3.5 * spaces[0].hWidth_accessor();

	double yPulseUB = 3.5*spaces[0].hHeight_accessor();
	double yPulseLB = -3.5*spaces[0].hHeight_accessor();
	double xPulseUB = 3.5*spaces[0].hWidth_accessor();
	double xPulseLB = -3.5*spaces[0].hWidth_accessor();

	double accuracyY = (yPulseUB - yPulseLB) / 199;
	double accuracyX = (xPulseUB - xPulseLB) / 199;

	double negTwohWidthsq = -2 * spaces[0].hWidth_accessor()*spaces[0].hWidth_accessor();
	double twoVzIntDistsq = 2 * spaces[0].VzDist_accessor()*spaces[0].VzDist_accessor();
	double twoPIhWidthVzIntDist = 2 * M_PI*(spaces[0].hWidth_accessor()*spaces[0].VzDist_accessor());
	
	double x = xPulseLB + accuracyX / 2;
	double y = yPulseLB + accuracyY / 2; 
	for (phase_space pulse : spaces) {
		x = xPulseLB + accuracyX / 2;
		y = yPulseLB + accuracyY / 2;
		while (y < yPulseUB) {
			while (x < xPulseUB) {
				if (x + pulse.zC_accessor() > xSearchLB && x + pulse.zC_accessor() < xSearchUB && y + pulse.VzC_accessor() > ySearchLB && y + pulse.VzC_accessor() < ySearchUB) {
					graphingGrid[int(map(x + pulse.zC_accessor(), xSearchLB, xSearchUB, 0, modelingXRange - 1)+0.5)][int(map(y + pulse.VzC_accessor(), ySearchLB, ySearchUB, 0, modelingYRange - 1)+0.5)] += pulse.intensity_multiplier_accessor()*(accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
					valueHolder2 += pulse.intensity_multiplier_accessor()*(accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
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

void summing(phase_space space, double graphingGrid[modelingXRange][modelingYRange]) {
	if (printStarts) {
		cout << "summing single phase space started" << endl;
	}
	double ySearchLB = -3.5 * space.hHeight_accessor();
	double ySearchUB = 3.5 *space.hHeight_accessor();
	double xSearchLB = -3.5 *space.hWidth_accessor();
	double xSearchUB = 3.5 *space.hWidth_accessor();
	double accuracyY = (ySearchUB - ySearchLB) / 199;  // Cannot be 99 or else for some reason the middle (25th row out of 50) takes half that of all the other rows. Look into?
	double accuracyX = (xSearchUB - xSearchLB) / 199;
	double negTwohWidthsq = -2 * space.hWidth_accessor()*space.hWidth_accessor();
	double twoVzIntDistsq = 2 * space.VzDist_accessor()*space.VzDist_accessor();
	double twoPIhWidthVzIntDist = 2 * M_PI*(space.hWidth_accessor()*space.VzDist_accessor());
	double x = xSearchLB + accuracyX / 2;
	double y = ySearchLB + accuracyY / 2;
	while (y < ySearchUB) {
		while (x < xSearchUB) {
			if (x + space.zC_accessor() > xSearchLB && x + space.zC_accessor() < xSearchUB && y + space.VzC_accessor() > ySearchLB && y + space.VzC_accessor() < ySearchUB) {
				graphingGrid[int(map(x, xSearchLB, xSearchUB, 0, modelingXRange - 1) + 0.5)][int(map(y, ySearchLB, ySearchUB, 0, modelingYRange - 1) + 0.5)] += (accuracyX*accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - space.chirp_accessor() * x)*(y - space.chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
				//valueHolder += accuracyX*accuracyY*(1 / (2 * M_PI))*exp((x*x + y * y) / (-2));
				valueHolder += (accuracyX*accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - space.chirp_accessor() * x)*(y - space.chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist));
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

void grid_subtraction(double grid1[modelingXRange][modelingYRange], double grid2[modelingXRange][modelingYRange]) {
	 static double finalGrid[modelingXRange][modelingYRange];
	 for (int i = 0; i < modelingXRange; i++){
		 for (int j = 0; j < modelingYRange; j++) {
			graphingMap3[i][j] = grid1[i][j] - grid2[i][j];
		 }
	 }
}


 double* t_g_s(double grid1[2][2], double grid2[2][2]) {
	 static double finalGrid[2][2];
	 for (int i = 0; i < modelingXRange; i++) {
		 for (int j = 0; j < modelingYRange; j++) {
			 finalGrid[i][j] = grid1[i][j] - grid2[i][j];
		 }
	 }
	 return *finalGrid;
 }

void write_to_file(double modelMatrix[modelingXRange][modelingYRange]) {//Untested
	if (printStarts)
		cout << "writing to file started" << endl;
	ofstream modelFile;
	modelFile.open("modeling_data.txt");
	for (int row = 0; row < modelingYRange; row++) {
		for (int col = 0; col < modelingXRange; col++) {
			modelFile << col << "\t" << row << "\t" << modelMatrix[col][row] << endl;
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