#include <iostream>
#include <vector>
#include "test_class.h"
#include "using_gnuplot.h"
#include "constants.h"
#include "modeling.h"
#include "data_processing.h"
#include "statistics.h"
#include "main_sequence.h"
#include <fstream>
//#include <iterator>
#include <map>

#include <mpi.h>

using namespace std;

double counter = 0.0;

PhaseSpace initialPulse(base_hWidth, base_hHeight, base_VzDist, base_zDist, base_chirp, base_b, base_pulseEnergy, base_intensityMultiplier,
						base_hDepth, base_hDepthVelocity, base_VxDist, base_xDist, base_chirpT, base_bT, 0.0, 0.0, 0.0, 0.0);

PhaseSpace finalPulse = initialPulse;

void mainSequence()
{
	//Initialize MPI environment
	MPI_Init(&argc, &argv);

	//Obtain ID and world size
	int P;
	int ID;
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);

	MPI_Win win;
	int *shared;
	if (ID == 0)
	{
		MPI_Win_allocate_shared(sizeof(double) * splitNumber, sizeof(double), MPI_INFO_NULL,
								MPI_COMM_WORLD, &shared, &win);
	}
	else
	{
		int disp_unit;
		MPI_Aint ssize;
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL,
								MPI_COMM_WORLD, &shared, &win);
		MPI_Win_shared_query(win, 0, &ssize, &disp_unit, &shared);
	}

	if (loadData != true)
	{
		vector<pair<double, double>> deviations;
		std::map<double, double> chirps;
		//map<double, double> chirp;
		PhaseSpace modifiedPulse = initialPulse;

		vector<vector<double>> specimen;
		readSpec("Data files/hexogon BN-powder-eels.sl0", specimen);
		normalizeSpecimen(specimen);

		for (int q = 0; q < 50; ++q)
		{
			valueHolder5 = 0;
			initialPulse = modifiedPulse;

			vector<PhaseSpace> splitPulses = initialPulse.split();
			vector<vector<PhaseSpace>> allPulses;

			for (PhaseSpace p : splitPulses)
			{
				allPulses.push_back(p.shatter(specimen));
			}

			double power = pow(2 * allPulses[0][0].getChirpT() / MAG_LENS_COEFFICIENT, 0.5);

			allPulses = analyzer(allPulses);

			double evolutionValue = 75 * q;
			for (vector<PhaseSpace> &pulses : allPulses)
			{
				for (PhaseSpace &pulse : pulses)
				{
					pulse.mag_lens(power);
					pulse.evolution(evolutionValue);
				}
			}

			vector<double> pixelArray(pixels);
			vector<double> base(pixels);
			pixelSum(pixelArray, allPulses);
			pixelSum(base, -113.5, -625, specimen);
			energyModeling(pixelArray);
			//specModeling(base, -625, -113.5);

			deviations.push_back({allPulses[0][0].getChirpT(), measureDeviation(base, pixelArray)});
			chirps.insert({allPulses[0][0].getHDepth(), measureDeviation(base, pixelArray)});
		}

		vector<pair<double, double>> chirpsInput;
		for (auto itr = chirps.begin(); itr != chirps.end(); ++itr)
		{
			chirpsInput.push_back({itr->first, itr->second});
		}
		deviationModeling(deviations, deviations[0].first, deviations[deviations.size() - 1].first, "Slope vs Deviation", "Slope", "Deviation (eV)");
		deviationModeling(chirpsInput, chirpsInput[0].first, chirpsInput[chirpsInput.size() - 1].first, "Width vs. Deviation", "Width (micrometer)", "Deviation (eV)");
	}
	else if (loadData)
	{ //Redundant - else would do just as fine as else if, but else if makes the logic easier to understand
		read_from_file("modeling_data.txt");
		finalDataOutput();
	}

	cout << "Code ran to end" << endl;
	pause();
}

PhaseSpace returnInitialPS()
{
	return initialPulse;
}

PhaseSpace returnFinalPS()
{
	return finalPulse;
}
