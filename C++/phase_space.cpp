#include <iostream>
#include "phase_space.h"
#include <math.h>
#include <functional>
#include <vector>
#include "external_variables.h"
#include "constants.h"
#include "using_gnuplot.h"
//#include "cpp_dec_float.hpp"
//#include "gnuplot_i.hpp"
//#include "gnuplot-iostream.h"
//#include <boost/tuple/tuple.hpp>
using namespace std;

phase_space::phase_space(){//For debugging only
	cout << "empty phase space initialized" << endl;
}

phase_space::phase_space(double test){//For debugging only
	b = test;
	cout << "test phase space initialized, b: " << b << endl;
}

phase_space::phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityMultiplierC,
						 double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double VzCC, double zCC, double xCC)
: hWidth(hWidthC), hHeight(hHeightC), VzDist(VzDistC), zDist(zDistC), chirp(chirpC), b(bC), pulseEnergy(pulseEnergyC), intensityMultiplier(intensityMultiplierC),
hDepth(hDepthC), hDepthVel(hDepthVelC), VxDist(VxDistC), xDist(xDistC), chirpT(chirpTC), bT(bTC), VzC(VzCC), zC(zCC), xC(zCC) {}

phase_space::phase_space(vector<phase_space> spaces):hWidth(0),hHeight(0),VzDist(0),zDist(0),chirp(0),b(0),pulseEnergy(0),intensityMultiplier(0),hDepth(0),hDepthVel(0),VxDist(0),xDist(0),chirpT(0),bT(0),VzC(0),zC(0),xC(0){
	for (int i = 0; i < splitNumber; i++) {
		VzC += spaces[i].VzC*spaces[i].intensityMultiplier;
		//DEBUGGING cout << spaces[i].intensityMultiplier << endl;
		zC += spaces[i].zC*spaces[i].intensityMultiplier;
		xC += spaces[i].xC*spaces[i].intensityMultiplier;
		intensityMultiplier += spaces[i].intensityMultiplier;
		pulseEnergy += spaces[i].pulseEnergy;
		//spaces[i].print();
	}
	//this.VzC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.VzC"(spaces))*this.intensityRatio;
	//int maxSize = sizeof(spaces[]);
	
	//for (phase_space space : spaces[splitNumber]){
	//	VzC += space.VzC*space.intensityMultiplier;
	//	zC += space.zC*space.intensityMultiplier;
	//	xC += space.xC*space.intensityMultiplier;
	//}
	//this.zC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.zC"(spaces));



		//Recalculation from single variable change		
	//HWIDTH BASED RECALCULATION
	//hWidth = spaces[0].hWidth; (Shown not to be accurate/gives wrong results)
	hWidth = originalHWidth + (spaces[0].hWidth-originalHWidth)*splitNumber;
	hHeight = spaces[0].hHeight * splitNumber; //If you add a * 1/chirp here sometimes it doesn't process it for some reason, a 1/chirp isn't needed here anyway but its an odd mystery why its only sometimes processed
	zDist = spaces[0].zDist;
	VzDist = zDist*hHeight/hWidth;
	chirp = VzDist*sqrt((1/pow(zDist,2))-(1/pow(hWidth,2)));
	b = chirp*pow(zDist/VzDist,2);
	hDepth = spaces[0].hDepth;
	hDepthVel = spaces[0].hDepthVel;
	VxDist = spaces[0].VxDist;
	xDist = spaces[0].xDist;
	chirpT = spaces[0].chirpT;
	bT = spaces[0].bT;

	//pulseEnergy = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.totalPulseEnergy"(spaces));
	//intensityMultiplier = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.intensityRatio"(spaces));
}
vector<phase_space> phase_space::split(){
	vector<phase_space> splitSpaces;
	double splitVzDist = 0;
	double splitChirp = 0;
	double splitB = 0;
	originalHWidth = hWidth;
	//phaseSpaces.length = to!int(spaces);
	double intensityMultipliers [splitNumber];
	for (int i = 0; i < splitNumber; i++){
		intensityMultipliers[i] = get_split_intensity_multiplier(splitNumber, i+1, hHeight, hWidth);
		//intensityMultipliers[i] = 0.01;
	}
	splitVzDist = VzDist/splitNumber;
	splitChirp = VzDist*sqrt((1/pow(zDist,2))-(1/pow(hWidth,2)));
	splitB = chirp*pow(zDist/VzDist,2);
	//i, ref elem; phaseSpaces
	for(int j = 0; j < splitNumber; j++){
		splitSpaces.push_back(phase_space(hWidth, hHeight/splitNumber, splitVzDist, zDist, splitChirp, splitB, pulseEnergy*intensityMultipliers[j], intensityMultipliers[j],
																		   hDepth, hDepthVel, VxDist, xDist, chirpT, bT, hHeight-(hHeight*2/splitNumber)*(double(j)+0.5), zC, xC));
	}
	phaseSpaces += splitNumber;
	return splitSpaces;
}

double phase_space::get_split_intensity_multiplier(double numSections, double sectionNum, double hHeight, double hWidth){
	//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
	//search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space (equal to 1)	
	double ySearchLB = -5.803*hHeight + ((5.803*hHeight*2.0/numSections)*(sectionNum-1));
	double ySearchUB = 5.803*hHeight - (5.803*hHeight*2.0/numSections)*(numSections - sectionNum);
	double xSearchLB = -5.803*hWidth;
	double xSearchUB = 5.803*hWidth;
	double accuracyY = (ySearchUB-ySearchLB)/99;
	double accuracyX = (xSearchUB-xSearchLB)/99;
	//double x = xSearchLB;
	//double y = ySearchLB;
	double x = xSearchLB+accuracyX/2;
	double y = ySearchLB+accuracyY/2;
	double intensityMultiplier = 0;
	double negTwohWidthsq = -2*hWidth*hWidth;
	double twoVzIntDistsq = 2*VzDist*VzDist;
	double twoPIhWidthsqVzIntDistsq = 2*M_PI*(hWidth*hWidth*VzDist*VzDist);
	if(numSections==sectionNum){
		ySearchUB += 0.000000000001;
	}
	while(y < ySearchUB-0.000000000001){
		while(x < xSearchUB){
			intensityMultiplier += accuracyX*accuracyY*exp((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
			/*DEBUGGING to test how testMax (which is the maximum value for each of the coordinates that intensityMultiplier checks for) 
				compares to valueHolder (the value at [0,0]- the absolute maximum of the gaussian function)
			if(accuracyX*accuracyY*exp((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq) > testMax){
				testMax = accuracyX*accuracyY*exp((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
				testMaxXCoordinate = x;
				testMaxYCoordinate = y;
				valueHolder = accuracyX*accuracyY*exp((0/(negTwohWidthsq))-((0-chirp*0)*(0-chirp*0)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
			}*/
			x += accuracyX;
		}
		y += accuracyY;
		x = xSearchLB+accuracyX/2;
	}
	//DEBUGGING	cout << "x: " << x + accuracyX*49 << endl; //Not 50 because x already starts on the first one, so plus 50 would make it the 51st slot
	return intensityMultiplier/2000;		
}

phase_space phase_space::evolution(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
	b += time;
	bT += time;
	if(chirp>0){
		VzDist = sqrt(1/((1/pow(hHeight,2))+pow((b/zDist),2)));
	}
	if(chirpT>0){
		VxDist = sqrt(1/((1/pow(hDepthVel,2))+pow((bT/xDist),2)));
	}
	if(chirp<0){
		zDist = b/(sqrt((1/pow(VzDist,2))-(1/pow(hHeight,2))));
	}
	if(chirpT<0){
		xDist = bT/(sqrt((1/pow(VxDist,2))-(1/pow(hDepthVel,2))));
	}
	chirp = b*pow(VzDist/zDist,2);
	chirpT = bT*pow(VxDist/xDist,2);
	hWidth = sqrt(1/((1/pow(zDist,2))-pow(chirp/VzDist,2)));
	hDepth = sqrt(1/((1/pow(xDist,2))-pow(chirpT/VxDist,2)));
	return *this;
}

phase_space phase_space::RFLens(double changeChirp){
	chirp += changeChirp;
	zDist = sqrt(1/((1/pow(hWidth,2))+pow(chirp/VzDist,2)));
	b = chirp*pow(zDist/VzDist,2);
	hHeight = sqrt(1/((1/pow(VzDist,2))-pow(b/zDist,2)));
	return *this;
}

phase_space phase_space::mag_lens(double changechirpT){
	chirpT += changechirpT;
	xDist = sqrt(1/((1/pow(hDepth,2))+pow(chirpT/VxDist,2)));
	bT = chirpT*pow(xDist/VxDist,2);
	hDepthVel = sqrt(1/((1/pow(VxDist,2))-pow(bT/xDist,2)));
	return *this;
}

phase_space phase_space::spectroscopy_function(){
	xC = xC + 7172.99042634*VzC;
	return *this;
}
//Conservation Checking - Emittence based
bool phase_space::longitudinal_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW){//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = hDimensionW*distH;
	double consValue2 = hDimensionH*distW;
	if(consValue > 0.000499 && consValue < 0.000501  && consValue2 > 0.000499 && consValue2 < 0.000501){//randomly decided range to account for data error
		return true;
	}
	else{
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;;
		return false;
	}
}

bool phase_space::transverse_area_conservation(double hDimensionW, double intDistH, double hDimensionH, double intDistW){//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = hDimensionW*intDistH;
	double consValue2 = hDimensionH*intDistW;
	if(consValue > 0.00299 && consValue < 0.00301  && consValue2 > 0.00299 && consValue2 < 0.00301){//randomly decided range to account for data error
		return true;
	}
	else{
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;
		return false;
	}
}

void phase_space::valid_variables_check(){
	cout << "Checks" << endl;
	cout << "hWidth: " << hWidth << endl;
	cout << "hWidth 1: " << sqrt(1/((1/(zDist*zDist)) - (chirp/VzDist)*(chirp/VzDist))) << endl;
}

void phase_space::print(){
	cout << "hWidth: " << hWidth << "  " <<    " hDepth: " << hDepth << endl;;
	cout << "hHeight: " << hHeight << " " <<   " hDepthVelocity: " << hDepthVel << endl;
	cout << "VzIntDist: " << VzDist <<    " VxIntDist: " << VxDist << endl;
	cout << "zIntDist: " << zDist << " " << " xIntDist: " << xDist << endl;
	cout << "chirp: " << chirp << " " <<       " chirpT: " << chirpT << endl;
	cout << "b: " << b << "     " <<           " bT: " << bT << endl;
	cout << "VzC: " << VzC << "   " <<         "zC: " << zC << endl;
	cout << "xC: " << xC << endl;
	cout << "pulseEnergy: " << pulseEnergy << endl;
	cout << "intensityMultiplier: " << intensityMultiplier << endl;
	//cout << "totalintensityMultiplier from PixelSum: " << sumUp());
	cout << endl;
	cout << "Longitudinal Emmittence Conserved: " << longitudinal_area_conservation(hWidth, VzDist, hHeight, zDist) << endl;
	cout << "Transverse Emmittence Conserved: " << transverse_area_conservation(hDepth, VxDist, hDepthVel, xDist) << endl;
	cout << endl;
	//valid_variables_check();
}

//DEBUGGING FUNCTIONS

//Accessor methods
double phase_space::hWidth_accessor(){return hWidth;}
double phase_space::hHeight_accessor(){return hHeight;}
double phase_space::VzDist_accessor(){return VzDist;}
double phase_space::zDist_accessor(){return zDist;}
double phase_space::chirp_accessor(){return chirp;}
double phase_space::b_accessor(){return b;}
double phase_space::intensity_multiplier_accessor(){return intensityMultiplier;}
double phase_space::VzC_accessor(){return VzC;}
double phase_space::zC_accessor(){return zC;}
double phase_space::xC_accessor() { return xC; }