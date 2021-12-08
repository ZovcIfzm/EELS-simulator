#include <iostream>
#include "PhaseSpace.h"
#include <math.h>
#include <vector>
#include "statistics.h"
//#include "data_processing.h"
//#include "using_gnuplot.h"
//#include "cpp_dec_float.hpp"
//#include "gnuplot_i.hpp"
//#include "gnuplot-iostream.h"
//#include <boost/tuple/tuple.hpp>
using namespace std;

PhaseSpace::PhaseSpace(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityMultiplierC,
					   double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double VzCC, double zCC, double VxCC, double xCC)
	: hWidth(hWidthC), hHeight(hHeightC), VzDist(VzDistC), zDist(zDistC), chirp(chirpC), b(bC), pulseEnergy(pulseEnergyC), intensityMultiplier(intensityMultiplierC),
	  hDepth(hDepthC), hDepthVel(hDepthVelC), VxDist(VxDistC), xDist(xDistC), chirpT(chirpTC), bT(bTC), VzC(VzCC), zC(zCC), VxC(VxCC), xC(zCC) {}

PhaseSpace::PhaseSpace(vector<PhaseSpace> spaces) : hWidth(0), hHeight(0), VzDist(0), zDist(0), chirp(0), b(0), pulseEnergy(0), intensityMultiplier(0), hDepth(0), hDepthVel(0), VxDist(0), xDist(0), chirpT(0), bT(0), VzC(0), zC(0), VxC(0), xC(0)
{
	for (int i = 0; i < splitNumber; i++)
	{
		VzC += spaces[i].VzC * spaces[i].intensityMultiplier;
		//DEBUGGING cout << spaces[i].intensityMultiplier << endl;
		zC += spaces[i].zC * spaces[i].intensityMultiplier;
		xC += spaces[i].xC * spaces[i].intensityMultiplier;
		intensityMultiplier += spaces[i].intensityMultiplier;
		pulseEnergy += spaces[i].pulseEnergy;
		//spaces[i].print();
	}
	//this.VzC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.VzC"(spaces))*this.intensityRatio;
	//int maxSize = sizeof(spaces[]);

	//for (PhaseSpace space : spaces[splitNumber]){
	//	VzC += space.VzC*space.intensityMultiplier;
	//	zC += space.zC*space.intensityMultiplier;
	//	xC += space.xC*space.intensityMultiplier;
	//}
	//this.zC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.zC"(spaces));

	//Recalculation from single variable change
	//HWIDTH BASED RECALCULATION
	//hWidth = spaces[0].hWidth; (Shown not to be accurate/gives wrong results)
	hWidth = originalHWidth + (spaces[0].hWidth - originalHWidth) * splitNumber;
	hHeight = spaces[0].hHeight * splitNumber; //If you add a * 1/chirp here sometimes it doesn't process it for some reason, a 1/chirp isn't needed here anyway but its an odd mystery why its only sometimes processed
	zDist = spaces[0].zDist;
	VzDist = zDist * hHeight / hWidth;
	chirp = VzDist * sqrt((1 / pow(zDist, 2)) - (1 / pow(hWidth, 2)));
	b = chirp * pow(zDist / VzDist, 2);
	hDepth = spaces[0].hDepth;
	hDepthVel = spaces[0].hDepthVel;
	VxDist = spaces[0].VxDist;
	xDist = spaces[0].xDist;
	chirpT = spaces[0].chirpT;
	bT = spaces[0].bT;

	//pulseEnergy = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.totalPulseEnergy"(spaces));
	//intensityMultiplier = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.intensityRatio"(spaces));
}

vector<PhaseSpace> PhaseSpace::split()
{
	vector<PhaseSpace> splitSpaces;
	originalHWidth = hWidth;
	//phaseSpaces.length = to!int(spaces);
	double intensityMultipliers[splitNumber];
	for (int i = 0; i < splitNumber; i++)
	{
		intensityMultipliers[i] = get_intensity(splitNumber, i + 1);
	}
	double splitHHeight = 0;
	double splitVzDist = 0;
	double splitB = 0;
	double splitZDist = 0;
	double splitChirp = 0;
	splitHHeight = hHeight / splitNumber;
	splitVzDist = VzDist / splitNumber;
	//b = chirp * pow(zDist / VzDist, 2);
	splitB = zDist * sqrt((1 / pow(splitVzDist, 2)) - (1 / pow(splitHHeight, 2)));
	splitZDist = splitB / (sqrt((1 / pow(splitVzDist, 2)) - (1 / pow(splitHHeight, 2))));
	splitChirp = splitVzDist * sqrt((1 / pow(zDist, 2)) - (1 / pow(hWidth, 2)));
	//i, ref elem; phaseSpaces

	for (int j = 0; j < splitNumber; j++)
	{
		splitSpaces.push_back(PhaseSpace(hWidth, splitHHeight, splitVzDist, splitZDist, splitChirp, splitB, pulseEnergy * intensityMultipliers[j], intensityMultipliers[j],
										 hDepth, hDepthVel, VxDist, xDist, chirpT, bT, hHeight - (hHeight * 2 / splitNumber) * (double(j) + 0.5), (hHeight - (hHeight * 2 / splitNumber) * (double(j) + 0.5)) / chirp, VxC, xC));
	}
	phaseSpaces += splitNumber;
	return splitSpaces;
}

vector<PhaseSpace> PhaseSpace::shatter(vector<vector<double>> spectroTable)
{
	int spaces = spectroTable.size();
	vector<PhaseSpace> shatteredPulses;
	double newVzDist = VzDist / spaces;
	double newChirp = newVzDist * sqrt((1 / pow(zDist, 2)) - (1 / pow(hWidth, 2)));
	double newB = newChirp * pow(zDist / newVzDist, 2);
	for (int i = 0; i < spaces; i++)
	{
		//cout << spectroTable[i][1] << endl;
		shatteredPulses.push_back(PhaseSpace(hWidth, hHeight / spaces, newVzDist, zDist, newChirp, newB, pulseEnergy * spectroTable[i][1] / baseTotal, intensityMultiplier * spectroTable[i][1],
											 hDepth, hDepthVel, VxDist, xDist, chirpT, bT, hHeight + (spectroTable[i][0] / 1117), zC, VxC, xC));
	}
	return shatteredPulses;
}

double PhaseSpace::intensity(double x, double y)
{
	double negTwohWidthsq = -2 * hWidth * hWidth;
	double twoVzIntDistsq = 2 * VzDist * VzDist;
	double twoPIhWidthVzIntDist = 2 * M_PI * (hWidth * VzDist);
	return exp((x * x / (negTwohWidthsq)) - ((y - chirp * x) * (y - chirp * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist);
}

double PhaseSpace::get_intensity(double numSections, double sectionNum)
{
	//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
	//search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space (equal to 1)
	double ySearchLB = -catchFactor * hHeight + ((catchFactor * hHeight * 2.0 / numSections) * (sectionNum - 1));
	double ySearchUB = catchFactor * hHeight - (catchFactor * hHeight * 2.0 / numSections) * (numSections - sectionNum);
	double xSearchLB = -catchFactor * hWidth;
	double xSearchUB = catchFactor * hWidth;

	//Convert to intensity_integration parameters
	double xHalfRange = (xSearchUB - xSearchLB) / 2;
	double yHalfRange = (ySearchUB - ySearchLB) / 2;

	double xOffset = (xSearchUB + xSearchLB) / 2;
	double yOffset = (ySearchUB + ySearchLB) / 2;

	return intensity_integration(xHalfRange, yHalfRange, xOffset, yOffset);
}

double PhaseSpace::x_integration(double xLeftLim, double xRightLim)
{
	double accuracyY = 2 * hHeight / 199;
	double accuracyX = (xRightLim - xLeftLim) / 199;
	double x = xLeftLim;
	double y = -hHeight + accuracyY / 2;

	double intensityValue = 0;

	while (y < hHeight)
	{
		while (x < xRightLim)
		{
			intensityValue += accuracyX * accuracyY * intensity(x - xC, y);
			x += accuracyX;
		}
		y += accuracyY;
		double x = xLeftLim;
	}
	return intensityValue;
}

double PhaseSpace::intensity_integration(double xHalfRange, double yHalfRange, double xOffset, double yOffset)
{
	double accuracyY = 2 * yHalfRange / 199;
	double accuracyX = 2 * xHalfRange / 199;
	double x = -xHalfRange + xOffset + accuracyX / 2;
	double y = -yHalfRange + yOffset + accuracyY / 2;

	double intensityValue = 0;

	while (y < yHalfRange + yOffset)
	{
		while (x < xHalfRange + xOffset)
		{
			intensityValue += accuracyX * accuracyY * intensity(x, y);
			x += accuracyX;
		}
		y += accuracyY;
		x = -xHalfRange + xOffset + accuracyX / 2;
	}
	return intensityValue;
}

void PhaseSpace::grid_integration(double xHalfRange, double yHalfRange, double xOffset, double yOffset, double grid[modelingXRange][modelingYRange], double xGridHalfRange, double yGridHalfRange)
{
	double accuracyY = 2 * yHalfRange / 199;
	double accuracyX = 2 * xHalfRange / 199;
	double x = -xHalfRange + xOffset + accuracyX / 2;
	double y = -yHalfRange + yOffset + accuracyY / 2;

	double intensityValue = 0;

	while (y < yHalfRange + yOffset)
	{
		while (x < xHalfRange + xOffset)
		{
			grid[int(map(x, -xGridHalfRange, xGridHalfRange, 0, double(modelingXRange) - 1) + 0.5)][int(map(y, -yGridHalfRange, yGridHalfRange, 0, double(modelingYRange) - 1) + 0.5)] += accuracyX * accuracyY * intensity(x, y);
			x += accuracyX;
		}
		y += accuracyY;
		x = -xHalfRange + xOffset + accuracyX / 2;
	}
}

PhaseSpace PhaseSpace::evolution(double dist)
{ //--To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
	//A divided by 1E6 was found in D code. Reason is unknown.
	double postTime = dist / 164.35;
	b += postTime;
	bT += postTime;

	if (chirp > 0)
	{
		VzDist = sqrt(1 / ((1 / pow(hHeight, 2)) + pow((b / zDist), 2)));
	}
	if (chirpT > 0)
	{
		VxDist = sqrt(1 / ((1 / pow(hDepthVel, 2)) + pow((bT / xDist), 2)));
	}
	if (chirp < 0)
	{
		zDist = b / (sqrt((1 / pow(VzDist, 2)) - (1 / pow(hHeight, 2))));
	}
	if (chirpT < 0)
	{
		xDist = bT / (sqrt((1 / pow(VxDist, 2)) - (1 / pow(hDepthVel, 2))));
	}

	//VzDist = sqrt(1/((1/pow(hHeight,2))+pow((b/zDist),2)));
	//VxDist = sqrt(1/((1/pow(hDepthVel,2))+pow((bT/xDist),2)));

	chirp = b * pow(VzDist / zDist, 2);
	chirpT = bT * pow(VxDist / xDist, 2);
	hWidth = sqrt(1 / ((1 / pow(zDist, 2)) - pow(chirp / VzDist, 2)));
	hDepth = sqrt(1 / ((1 / pow(xDist, 2)) - pow(chirpT / VxDist, 2)));

	zC += VzC * postTime;
	xC += hDepthVel * postTime;
	return *this;
}

PhaseSpace PhaseSpace::RFLens(double power)
{
	double tempSlope = VzC / zC;
	tempSlope -= sqrt(power) * RF_LENS_COEFFICIENT;
	VzC = tempSlope * zC;

	chirp -= sqrt(power) * RF_LENS_COEFFICIENT;
	zDist = sqrt(1 / ((1 / pow(hWidth, 2)) + pow(chirp / VzDist, 2)));
	b = chirp * pow(zDist / VzDist, 2);
	hHeight = sqrt(1 / ((1 / pow(VzDist, 2)) - pow(b / zDist, 2)));
	return *this;
}

PhaseSpace PhaseSpace::mag_lens(double power)
{
	//A divided by 1E12 for power was found in D code. Reason is unknown.
	//double tempSlope = VxC / xC;
	//tempSlope -= pow(power, 2) * MAG_LENS_COEFFICIENT;
	//VxC = tempSlope * xC;
	//xC = VxC / tempSlope;

	chirpT -= pow(power, 2) * MAG_LENS_COEFFICIENT;
	xDist = sqrt(1 / ((1 / pow(hDepth, 2)) + pow(chirpT / VxDist, 2)));
	bT = chirpT * pow(xDist / VxDist, 2);
	hDepthVel = sqrt(1 / ((1 / pow(VxDist, 2)) - pow(bT / xDist, 2)));
	return *this;
}

PhaseSpace PhaseSpace::spectroscopy_function()
{
	xC = xC + 7172.99042634 * VzC;
	return *this;
}

//Accessor methods
double PhaseSpace::getHWidth() { return hWidth; }
double PhaseSpace::getHHeight() { return hHeight; }
double PhaseSpace::getVzDist() { return VzDist; }
double PhaseSpace::getZDist() { return zDist; }
double PhaseSpace::getChirp() { return chirp; }
double PhaseSpace::getB() { return b; }
double PhaseSpace::getIntensityMultiplier() { return intensityMultiplier; }
double PhaseSpace::getVzC() { return VzC; }
double PhaseSpace::getZC() { return zC; }

//Longitudinal accessor methods
double PhaseSpace::getHDepth() { return hDepth; }
double PhaseSpace::getHDepthVel() { return hDepthVel; }
double PhaseSpace::getVxDist() { return VxDist; }
double PhaseSpace::getXDist() { return xDist; }
double PhaseSpace::getChirpT() { return chirpT; }
double PhaseSpace::getBT() { return bT; }
double PhaseSpace::getVxC() { return VxC; }
double PhaseSpace::getXC() { return xC; }