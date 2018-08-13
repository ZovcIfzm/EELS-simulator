#include <iostream>
#include <vector>
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
#include "modeling.h"
using namespace std;

//TODO
//Test limits of roundoff error for variable checking
//
//

int main(){
	if (loadData != true) {
		phase_space initialPulse(base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp, base_b, base_pulseEnergy, base_intensityMultiplier,
			base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0);

		if (printInitialPhaseSpace)
			initialPulse.print();

		initialPulse.print();
		//initialPulse.evolution(100);
		//initialPulse.print();
		initialPulse.RFLens(1);
		vector<phase_space> pulses = initialPulse.split();
		//pulses[0].print();
		//split_evolution(pulses, 5000.0);
		//split_evolution(pulses, 0.0);
		//pulses[0].print();
		//phase_space recombinantSpace(pulses);

		if (modelFinalInitialPhaseSpace){
			initialPulse.print();
			summing(initialPulse);
			modeling(graphingMap);
		}
		if (modelSplitSpaces) {
			summing(pulses);
			modeling(graphingMap);
		}
			//if (modelRecombinantPhaseSpace)
		//	recombinantSpace.print();
		if (saveData)
			write_to_file(graphingMap);
	}
	else if (loadData)//Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
	// DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	//cout << testMax << endl;
	//cout << testMaxXCoordinate << endl;
	//cout << testMaxYCoordinate << endl;
	//cout << valueHolder << endl;

	cout << "Code ran to end" << endl;
	pause();
	//gnuplot6();
}