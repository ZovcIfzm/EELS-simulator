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
using namespace std;

double counter = 0.0;

PhaseSpace initialPulse(base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp, base_b, base_pulseEnergy, base_intensityMultiplier,
	base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0, 0.0);

PhaseSpace finalPulse = initialPulse;

void mainSequence() {
	if (loadData != true) {
		PhaseSpace modifiedPulse = initialPulse;

		vector<vector<double>> specimen;
		readSpec("Data files/hexogon BN-powder-eels.sl0", specimen);
		normalizeSpecimen(specimen);

		
		vector<PhaseSpace> splitPulses = initialPulse.split();
		vector<vector<PhaseSpace>> allPulses;
		cout << "shatter count: ";
		int a = 0;
		for (PhaseSpace p : splitPulses) {
			cout << ++a << ", ";
			allPulses.push_back(p.shatter(specimen));
		}

		vector<PhaseSpace> singleShatteredPulse = initialPulse.shatter(specimen);

		//cout << "chirpT: " << allPulses[0][0].getChirpT() << endl;
		//cout << "HDepthVel: " << allPulses[0][0].getHDepthVel() << endl;
		//cout << "HDepth: " << allPulses[0][0].getHDepth() << endl;
		//cout << "bT: " << allPulses[0][0].getBT() << endl;

		double power = pow(2*allPulses[0][0].getChirpT() / MAG_LENS_COEFFICIENT, 0.5);
		for (vector<PhaseSpace> &pulses : allPulses) {
			for (PhaseSpace &pulse : pulses) {
				//pulse.mag_lens(power);
		//		pulse.evolution(1.615E2);
			}
		}
		double power2 = pow(2 * singleShatteredPulse[0].getChirpT() / MAG_LENS_COEFFICIENT, 0.5);
		for (PhaseSpace& pulse : singleShatteredPulse) {
			//pulse.mag_lens(power2);
		}

		double lowestXC;
		double highestXC;

		lowestXC = allPulses[allPulses.size() - 1][allPulses[0].size() - 1].getXC();
		highestXC = allPulses[0][0].getXC();

		cout << "chirpT: " << allPulses[0][0].getChirpT() << endl;
		cout << "true chirpT?: " << allPulses[0][0].getHDepthVel() / allPulses[0][0].getHDepth() << endl;
		cout << "VxDist: " << allPulses[0][0].getVxDist() << endl;
		cout << "HDepth: " << allPulses[0][0].getHDepth() << endl;
		cout << "bT: " << allPulses[0][0].getBT() << endl;
		cout << "lowest and highest xC: " << lowestXC << " | " << highestXC << endl;
	

		/*	TESTING EVOLUTION
		summing(initialPulse, graphingMap);
		double width = initialPulse.getHWidth();
		double height = initialPulse.getHHeight();
		double a = initialPulse.intensity_integration(width, height, 0, 0);

		initialPulse = initialPulse.evolution(1000000);
			
		double b = initialPulse.intensity_integration(width, height, 0, 0);

		cout << "lah: " << a - b << endl;
		*/

		//summing(initialPulse, graphingMap2);
		
		//modeling(graphingMap);
		//modeling(graphingMap2);

		//cout << "deviation: " << measureDeviation(graphingMap, graphingMap2) << endl;

		
		//cout << "intialPulse B: " << initialPulse.getB() << endl;
		//initialPulse = initialPulse.evolution(1E5);
		//cout << "intialPulse B: " << initialPulse.getB() << endl;
		
		singleShatteredPulse = analyzer(singleShatteredPulse);
		allPulses = analyzer(allPulses);

		lowestXC = allPulses[allPulses.size() - 1][allPulses[0].size() - 1].getXC();
		highestXC = allPulses[0][0].getXC();

		cout << "chirpT: " << allPulses[0][0].getChirpT() << endl;
		cout << "true chirpT?: " << allPulses[0][0].getHDepthVel() / allPulses[0][0].getHDepth() << endl;
		cout << "VxDist: " << allPulses[0][0].getVxDist() << endl;
		cout << "HDepth: " << allPulses[0][0].getHDepth() << endl;
		cout << "bT: " << allPulses[0][0].getBT() << endl;
		cout << "lowest and highest xC: " << lowestXC << " | " << highestXC << endl;

		for (vector<PhaseSpace>& pulses : allPulses) {
			for (PhaseSpace& pulse : pulses) {
				pulse.mag_lens(power);
				pulse.evolution(1.6164E3);
			}
		}

		for (PhaseSpace& pulse : singleShatteredPulse) {
			//pulse.evolution(1.6165E3);
		}

		lowestXC = allPulses[allPulses.size() - 1][allPulses[0].size() - 1].getXC();
		highestXC = allPulses[0][0].getXC();

		cout << "chirpT: " << allPulses[0][0].getChirpT() << endl;
		cout << "true chirpT?: " << allPulses[0][0].getHDepthVel() / allPulses[0][0].getHDepth() << endl;
		cout << "VxDist: " << allPulses[0][0].getVxDist() << endl;
		cout << "HDepth: " << allPulses[0][0].getHDepth() << endl;
		cout << "bT: " << allPulses[0][0].getBT() << endl;
		cout << "lowest and highest xC: " << lowestXC << " | " << highestXC << endl;

		cout << "point reached: " << endl;
		vector<double> pixelArray(pixels);
		vector<double> base(pixels);
		pixelSumMulti(pixelArray, allPulses);
		pixelSum(base, -113.5, -625, specimen);
		cout << pixelArray[0] << endl;
		specModeling(pixelArray);
		specModeling(base);
		cout <<"deviation: " << measureDeviation(base, pixelArray);
		
		/*
		for (vector<PhaseSpace> a : allPulses) {
			for (PhaseSpace p : a) {
				p.evolution(1000);
			}
		}*/
		/*print(allPulses[5][5]);
		for (int i = 0; i < allPulses.size(); i++) {
			for (int j = 0; j < allPulses[i].size(); j++) {
				allPulses[i][j] = allPulses[i][j].evolution(10000000000);
			}
		}
		print(allPulses[5][5]);*/

		//cout << "energy total: " << check_energy_conservation(allPulses);
		//vector<PhaseSpace> shatteredPulses = initialPulse.shatter(v);


		//summing(shatteredPulses, graphingMap3);
		//summing(initialPulse, graphingMap);
		//measureDeviation(graphingMap, graphingMap3);
		//psComparison(shatteredPulses[0], initialPulse);
		//grid_subtraction(graphingMap, graphingMap3, graphingMap2);
		//modeling(graphingMap);
		//modeling(graphingMap3);
		//for (int i = 0; i < modelingXRange; i++) {
		//	for (int j = 0; j < modelingYRange; j++) {
		//		valueHolder3 += graphingMap3[i][j];
		//	}
		//}

		//pixelSum(pixelArray, shatteredPulses);
		//specModeling(pixelArray);

	}
	else if (loadData) {//Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
		finalDataOutput();
	}
	// DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	//cout << testMax << endl;
	//cout << testMaxXCoordinate << endl;
	//cout << testMaxYCoordinate << endl;
	//cout << valueHolder << endl;

	cout << "Code ran to end" << endl;
	pause();
}

PhaseSpace returnInitialPS() {
	return initialPulse;
}

PhaseSpace returnFinalPS() {
	return finalPulse;
}

