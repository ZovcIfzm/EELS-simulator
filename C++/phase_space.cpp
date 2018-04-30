#include <iostream>
#include "phase_space.h"
#include <math.h>
using namespace std;

phase_space::phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityRatioC,
double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC) 
: hWidth(hWidthC), hHeight(hHeightC), VzDist(VzDistC), zDist(zDistC), chirp(chirpC), b(bC), pulseEnergy(pulseEnergyC), intensityRatio(intensityRatioC),
hDepth(hDepthC), hDepthVel(hDepthVelC), VxDist(VxDistC), xDist(xDistC), chirpT(chirpTC), bT(bTC), vZC(vZCC), zC(zCC), xC(zCC) {}

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
		xC = xC + 7172.99042634*vZC;
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
void phase_space::print(){
		cout << "hWidth: " << hWidth << "  " <<    " hDepth: " << hDepth << endl;;
		cout << "hHeight: " << hHeight << " " <<   " hDepthVelocity: " << hDepthVel << endl;
		cout << "VzIntDist: " << VzDist <<    " VxIntDist: " << VxDist << endl;
		cout << "zIntDist: " << zDist << " " << " xIntDist: " << xDist << endl;
		cout << "chirp: " << chirp << " " <<       " chirpT: " << chirpT << endl;
		cout << "b: " << b << "     " <<           " bT: " << bT << endl;
		cout << "VzC: " << vZC << "   " <<         "zC: " << zC << endl;
		cout << "xC: " << xC << endl;
		cout << "pulseEnergy: " << pulseEnergy << endl;
		cout << "intensityRatio: " << intensityRatio << endl;
		//cout << "totalIntensityRatio from PixelSum: " << sumUp());
		cout << endl;
		cout << "Longitudinal Emmittence Conserved: " << longitudinal_area_conservation(hWidth, VzDist, hHeight, zDist) << endl;
		cout << "Transverse Emmittence Conserved: " << transverse_area_conservation(hDepth, VxDist, hDepthVel, xDist) << endl;
		cout << endl;
	}
};