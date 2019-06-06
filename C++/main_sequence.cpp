#include <iostream>
#include <vector>
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
#include "modeling.h"
#include "data_processing.h"
#include "main_sequence.h"
using namespace std;

double counter = 0.0;

phase_space initialPulse( base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp, base_b, base_pulseEnergy, base_intensityMultiplier,
	base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0);

phase_space finalPulse = initialPulse;

void mainSequence() {
	if (loadData != true) {
		if (printInitialPhaseSpace)
			initialPulse.print();
		phase_space modifiedPulse = initialPulse;
		//initialPulse.evolution(100);
		//initialPulse.print();
		modifiedPulse.evolution(10);
		vector<phase_space> pulses = initialPulse.split();
		//modifiedPulse.evolution(-10);
		//pulses[0].print();
		//split_evolution(pulses, 5000.0);
		//split_evolution(pulses, 0.0);
		//pulses[0].print();
		//phase_space recombinantSpace(pulses);
		finalPulse = modifiedPulse;
		psComparison(initialPulse, modifiedPulse);
		if (printFinalPhaseSpace){
			finalPulse.print();
		}
		if (modelFinalPhaseSpace) {
			modifiedPulse.print();
			summing(modifiedPulse);
			modeling(graphingMap);
		}
		if (modelSplitSpaces) {
			summing(pulses);
			modeling(graphingMap);
		}
		//if (modelRecombinantPhaseSpace)
	//	recombinantSpace.print();
		if (loadValueHolders) {
			cout << "ValueHolder1: " << valueHolder << endl;
			cout << "ValueHolder2:" << valueHolder2 << endl;
		}
		if (saveData)
			write_to_file(graphingMap);
	}
	else if (loadData)//Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
	if (openFinalData)
		finalDataOutput();
	// DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	//cout << testMax << endl;
	//cout << testMaxXCoordinate << endl;
	//cout << testMaxYCoordinate << endl;
	//cout << valueHolder << endl;

	cout << "Code ran to end" << endl;
	if(pauseClose)
		pause();
	//gnuplot6();
}

phase_space returnInitialPS() {
	return initialPulse;
}

phase_space returnFinalPS() {
	return finalPulse;
}

