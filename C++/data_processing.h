#include <iostream>
#include "PhaseSpace.h"

void finalDataOutput();
void psComparison(PhaseSpace space, PhaseSpace space2);
void readSpec(string filename, vector<vector<double>>& v);
void outputPhaseSpace(ofstream& file, PhaseSpace pulse, string name);
void normalizeSpecimen(vector<vector<double>>& specimen);