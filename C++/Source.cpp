#include <iostream>
#include <vector>
#include "phase_space.h"
#include "test_class.h"
using namespace std;

int main(){
	phase_space initialPulse (1.82534513746, 0.00069511761, 0.00027392079, 0.71930270898, 0.00035, 2413.46744543, 100, 1.0, 0.60844837915, 0.00494641054, 0.00493057439, 0.60650040415, 0.00065, 9.83513928219, 0.0, 0.0, 0.0);
	
	//initialPulse.print();
	//initialPulse.evolution(1);
	//initialPulse.print();
	//initialPulse.RFLens(1);

	initialPulse.print();
	vector<phase_space> pulses = initialPulse.split();
	for(int i = 0; i < pulses.size(); i++){
		pulses[i].print();
		//cout << pulses[i].intensityMultiplierAccessor();
	}
	phase_space recombinantSpace(pulses);
	recombinantSpace.print();
	//pulses[0].print();
	//pulses[1].print();
	//pulses[2].print();
	//pulses[3].print();
	//pulses[4].print();

	cout << "Code ran to end" << endl;
}