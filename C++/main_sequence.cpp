#include <iostream>
#include <vector>
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
#include "modeling.h"
#include "data_processing.h"
#include "statistics.h"
#include "main_sequence.h"
#include <fstream>
//#include <iterator>
#include <map>
using namespace std;

double counter = 0.0;

PhaseSpace initialPulse(base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp, base_b, base_pulseEnergy, base_intensityMultiplier,
	base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0, 0.0);

PhaseSpace finalPulse = initialPulse;

void mainSequence() {
	if (loadData != true) {
		vector<pair<double,double>> deviations;
		std::map<double, double> chirps;
		//map<double, double> chirp;
		PhaseSpace modifiedPulse = initialPulse;

		vector<vector<double>> specimen;
		readSpec("Data files/hexogon BN-powder-eels.sl0", specimen);
		normalizeSpecimen(specimen);

		for (int q = 0; q < 50; ++q) {
			valueHolder5 = 0;
			initialPulse = modifiedPulse;
			cout << "beginning bT" << initialPulse.getBT() << endl;

			vector<PhaseSpace> splitPulses = initialPulse.split();
			vector<vector<PhaseSpace>> allPulses;

			for (PhaseSpace p : splitPulses) {
				allPulses.push_back(p.shatter(specimen));
			}

			double power = pow(2 * allPulses[0][0].getChirpT() / MAG_LENS_COEFFICIENT, 0.5);
			for (vector<PhaseSpace>& pulses : allPulses) {
				for (PhaseSpace& pulse : pulses) {
					//pulse.mag_lens(power);
			//		pulse.evolution(1.615E2);
				}
			}
			
			allPulses = analyzer(allPulses);
			double evolutionValue = 75*q;

			for (vector<PhaseSpace>& pulses : allPulses) {
				for (PhaseSpace& pulse : pulses) {
					pulse.mag_lens(power);
					pulse.evolution(evolutionValue);
				}
			}

			cout << "point reached: " << endl;
			vector<double> pixelArray(pixels);
			vector<double> base(pixels);
			pixelSumMulti(pixelArray, allPulses);
			pixelSum(base, -113.5, -625, specimen);
			//specModeling(pixelArray);
			//specModeling(base);
			cout << "deviation: " << measureDeviation(base, pixelArray) << endl;
			deviations.push_back({ evolutionValue, measureDeviation(base, pixelArray) });
			//cout << "chirps: " << allPulses[0][0].getChirpT() << endl;
			chirps.insert({ allPulses[0][0].getHDepthVel() / allPulses[0][0].getHDepth(), measureDeviation(base, pixelArray) });
			cout << "iteration: " << q << endl;
		}
				
		vector<pair<double, double>> chirpsInput;
		for (auto itr = chirps.begin(); itr != chirps.end(); ++itr) {
			chirpsInput.push_back({ itr->first, itr->second });
		}
		//deviationModeling(deviations, deviations.[0], deviations[deviations.size() - 1].first, "", "", "");
		deviationModeling(chirpsInput, chirpsInput[0].first, chirpsInput[chirpsInput.size()-1].first, "Slope vs. Deviation", "Slope", "Deviation");
	}
	else if (loadData) {//Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
		finalDataOutput();
	}

	cout << "Code ran to end" << endl;
	pause();
}

PhaseSpace returnInitialPS() {
	return initialPulse;
}

PhaseSpace returnFinalPS() {
	return finalPulse;
}

