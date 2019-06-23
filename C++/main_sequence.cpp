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
		phase_space modifiedPulse = initialPulse;
		//vector<phase_space> pulses = initialPulse.split();
		finalPulse = modifiedPulse;
		
		//summing(finalPulse, graphingMap);
		//summing(pulses, graphingMap2);
		
		//valueHolder4 = measureDeviation(graphingMap, graphingMap2);
		//vector<vector<double>> v;
		//readSpec("Data files/hexogon BN-powder-eels.sl0", v);
		summing(finalPulse, graphingMap);
		write_to_file(graphingMap);
		modeling(graphingMap);
	}
	else if (loadData)//Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
		//finalDataOutput();

	// DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	//cout << testMax << endl;
	//cout << testMaxXCoordinate << endl;
	//cout << testMaxYCoordinate << endl;
	//cout << valueHolder << endl;

	cout << "Code ran to end" << endl;
	if(pauseClose)
		pause();
}

phase_space returnInitialPS() {
	return initialPulse;
}

phase_space returnFinalPS() {
	return finalPulse;
}

