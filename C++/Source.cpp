#include <iostream>
#include <vector>
#include "phase_space.h"
#include "test_class.h"
#include "using_gnuplot.h"
using namespace std;

/*TODO
Test limits of roundoff error for variable checking

*/

int main(){
	phase_space initialPulse (base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp,  base_b, base_pulseEnergy, base_intensityMultiplier, 
		base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0);
	
	//initialPulse.print();
	//initialPulse.evolution(1);
	//initialPulse.print();
	//initialPulse.RFLens(1);

	initialPulse.print();
	//initialPulse.evolution(10000);
	vector<phase_space> pulses = initialPulse.split();
	/*for(phase_space pulse : pulses)
	{
		pulse.evolution(10000);
		cout <<  pulse.hWidth_accessor() << endl;
		//pulse.print();
	}*/
	for(int i = 0; i < splitNumber; i++)
	{
		pulses[i].evolution(10000);
	}
	/*pulses[0].evolution(10000);
	pulses[1].evolution(10000);
	pulses[2].evolution(10000);
	pulses[0].print();
	pulses[1].print();
	pulses[2].print();*/
	//initialPulse.RFLens(100);
	phase_space recombinantSpace(pulses);
	recombinantSpace.print();
	/* DEBUGGING in conjuction with code inside get_split_intensity_multiplier that assigns values to the printed variables
	cout << testMax << endl;
	cout << testMaxXCoordinate << endl;
	cout << testMaxYCoordinate << endl;
	cout << valueHolder << endl;
	*/
	cout << "Code ran to end" << endl;
	gnuplot6();
}