#include <iostream>
#include <vector>
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
#include "modeling.h"
using namespace std;

/*TODO
Test limits of roundoff error for variable checking

*/

int main(){
	phase_space initialPulse (base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp,  base_b, base_pulseEnergy, base_intensityMultiplier, 
		base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0);
	
	if (printInitialPhaseSpace)
		initialPulse.print();

	//initialPulse.print();
	//initialPulse.evolution(1);
	//initialPulse.print();
	//initialPulse.RFLens(100);
	vector<phase_space> pulses = initialPulse.split();
	summing(initialPulse);
	phase_space recombinantSpace(pulses);
	
	if (printFinalInitialPhaseSpace)
		initialPulse.print();
	if (printRecombinantPhaseSpace)
		recombinantSpace.print();
	if (saveData)
		writing_to_file(graphingMap);
	
	/* DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	cout << testMax << endl;
	cout << testMaxXCoordinate << endl;
	cout << testMaxYCoordinate << endl;
	cout << valueHolder << endl;
	*/
	cout << "Code ran to end" << endl;
	pause_method();
	//gnuplot6();
}