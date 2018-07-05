#include <iostream>
#include <vector>
#include "phase_space.h"
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
using namespace std;


/*TODO
Test limits of roundoff error for variable checking

*/
double mapping(double input, double lowerBound, double upperBound, double outputLowerBound, double outputUpperBound) {
	return (   (input - lowerBound)/(upperBound - lowerBound)*(outputUpperBound - outputLowerBound)   );
}

void split_evolution(vector<phase_space>& ev_pulses, double time) {
	for (int i = 0; i < splitNumber; i++)
	{
		ev_pulses[i].evolution(time);
	}
}

void summing(vector<phase_space> spaces) {
	double ySearchLB = -5.803*spaces[0].VzC_accessor();
	double ySearchUB = 5.803*spaces[splitNumber-1].VzC_accessor();
	double xSearchLB = -5.803*spaces[0].hWidth_accessor();
	double xSearchUB = 5.803*spaces[0].hWidth_accessor();
	double accuracyY = (ySearchUB - ySearchLB) / 99;
	double accuracyX = (xSearchUB - xSearchLB) / 99;
	//double x = xSearchLB;
	//double y = ySearchLB;
	double x = xSearchLB;
	double y = ySearchLB;
	double intensityMultiplier = 0;
	double negTwohWidthsq = -2 * spaces[0].hWidth_accessor()*spaces[0].hWidth_accessor();
	double twoVzIntDistsq = 2 * spaces[0].VzDist_accessor()*spaces[0].VzDist_accessor();
	double twoPIhWidthsqVzIntDistsq = 2 * M_PI*(spaces[0].hWidth_accessor()*spaces[0].hWidth_accessor()*spaces[0].VzDist_accessor()*spaces[0].VzDist_accessor());
	while (y < ySearchUB) {
		while (x < xSearchUB) {
			intensityMultiplier += accuracyX * accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthsqVzIntDistsq);
			if (accuracyX*accuracyY*exp((x*x / (negTwohWidthsq)) - ((y - spaces[0].chirp_accessor() * x)*(y - spaces[0].chirp_accessor() * x) / (twoVzIntDistsq))) / (twoPIhWidthsqVzIntDistsq) > testMax) {
				/*DEBUGGING to test how testMax (which is the maximum value for each of the coordinates that intensityMultiplier checks for)
				compares to valueHolder (the value at [0,0]- the absolute maximum of the gaussian function)
				testMax = accuracyX*accuracyY*exp((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
				testMaxXCoordinate = x;
				testMaxYCoordinate = y;
				valueHolder = accuracyX*accuracyY*exp((0/(negTwohWidthsq))-((0-chirp*0)*(0-chirp*0)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
				*/
			}
			x += accuracyX;
		}
		y += accuracyY;
		x = xSearchLB + accuracyX / 2;
	}
}

int main(){
	phase_space initialPulse (base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp,  base_b, base_pulseEnergy, base_intensityMultiplier, 
		base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0);
	
	//initialPulse.print();
	//initialPulse.evolution(1);
	//initialPulse.print();
	//initialPulse.RFLens(100);

	initialPulse.print();
	vector<phase_space> pulses = initialPulse.split();
	split_evolution(pulses, 10000);
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