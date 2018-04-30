#include "phase_space.h"
#include <math.h>
phase_space::phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityRatioC,
double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC) 
: hWidth(hWidthC), hHeight(hHeightC), VzDist(VzDistC), zDist(zDistC), chirp(chirpC), b(bC), pulseEnergy(pulseEnergyC), intensityRatio(intensityRatioC),
hDepth(hDepthC), hDepthVel(hDepthVelC), VxDist(VxDistC), xDist(xDistC), chirpT(chirpTC), bT(bTC), vZC(vZCC), zC(zCC), xC(zCC) {}
//	phase_space::hWidth = hWidthC;
//	hHeight, VzIntDist, zIntDist, chirp, b, pulseEnergy, intensityRatio,
//hDepth, hDepthVel, VxDist, xDist, chirpT, bT, vZC, zC, xC
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
};