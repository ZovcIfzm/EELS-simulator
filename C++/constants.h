#pragma once

//Defined variables (no need to be included in cpp file, they only need to be written once)
//Debugging
	//Printing
#define printStarts false
#define printEnds false
#define printInitialPhaseSpace true
#define printFinalPhaseSpace true
#define modelFinalPhaseSpace false
#define modelRecombinantPhaseSpace false
#define modelSplitSpaces false
#define saveData true
#define loadData false
#define loadValueHolders true
#define openFinalData false
#define pauseClose true

//Miscellaneous defined variables
#define splitNumber 1001 //Needs to be int, to be used as array size & needs to be odd in order for doubleintensityMultiplier towork (to sum up the center correctly, with even it's unsymmetric for some reason)
#define M_PI 3.14159265358979323846264338328
#define RF_LENS_COEFFICIENT 0.0002232305
#define MAG_LENS_COEFFICIENT 2800000000
#define modelingXRange 51 //Is 51 so that there is a middle number i.e. 0-24, 25, 26-50
#define modelingYRange 51 //Is 51 so that there is a middle number i.e. 0-24, 25, 26-50
#define pixels 500
#define magFactor 1.1

//Phase space parameters
#define  base_hWidth 1.82534513746
#define  base_hHeight 0.00069511761
#define  base_VzDist 0.00027392079
#define  base_zDist 0.71930270898
#define  base_chirp 0.00035
#define  base_b 2413.46744543
#define  base_pulseEnergy 100.0
#define  base_intensityMultiplier 1.0
#define  base_hDepth 0.60844837915
#define  base_hDepthVelocity 0.00494641054
#define  base_VxDist 0.00493057439
#define  base_xDist 0.60650040415
#define  base_chirpT 0.00065
#define  base_bT 9.83513928219
#define  base_vZC 0.0
#define  base_zC 0.0
#define  base_xC 0.0

//1/(2*3.14159*1.82534513746^2*0.00027392079^2)*e^(	(-(x^2)/(2*1.82534513746))-(y-0.00035*x)^2/(2*0.71930270898))

//Splitting parameters
#define catchFactor 4 //factor needed to be multipled to hWidth/hHeight during summing to get almost 100% of the population of the 2d gaussian

extern double originalHWidth;
extern double originalHWidth2;
extern int phaseSpaces;
extern double testMax;
extern double testMaxXCoordinate;
extern double testMaxYCoordinate;
extern double valueHolder1;
extern double valueHolder2;
extern double valueHolder3;
extern double valueHolder4;
extern double valueHolder5;
extern double valueHolder6;
extern double centralVzReference;

extern double graphingMap[modelingXRange][modelingYRange];
extern double graphingMap2[modelingXRange][modelingYRange];
extern double graphingMap3[modelingXRange][modelingYRange];

//shattered pulse total: 2890661135.000000

/*
PhaseSpace secondPulse(double a, double b, double c, double d, double e,
double f, double g, double h, double i, double j,
double k, double l, double m, double n, double o,
double p, double q);*/