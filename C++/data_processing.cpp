#include <iostream>
#include "data_processing.h"
#include <fstream>
#include <sstream>
#include <Windows.h>
#include "main_sequence.h"
#include <tuple>
#include "constants.h"
#include <string>

void readSpec(string filename, vector<vector<double>> &v) {
	ifstream dataFile(filename);

	ofstream dataOutput;
	dataOutput.open("testing.txt");

	string w1, w2;
	double v1, v2, v3;
	dataFile >> w1 >> w2;//removing the first two lines (first line is one character, second has no spaces)
	cout << "words: " << w1 << " " << w2 << endl;
	while (dataFile >> w1 >> w2 >> v3) {
		w1.pop_back();
		w2.pop_back();
		v1 = atof(w1.c_str());
		v2 = atof(w2.c_str());

		vector <double> values;
		values.push_back(v1);
		values.push_back(v2);
		values.push_back(v3);
		v.push_back(values);
		dataOutput << values[0] << ", " << values[1] << endl;
	}
	dataOutput.close();
}

void psComparison(phase_space space, phase_space space2) {
	int counter = 0;
	if (space.hWidth_accessor() / space2.hWidth_accessor() > 1.01 || space.hWidth_accessor() / space2.hWidth_accessor() < 0.99) {
		cout << space.hWidth_accessor() / space2.hWidth_accessor() << "% hWidth divergence" << endl;
		counter++;
	}
	if (space.hHeight_accessor() / space2.hHeight_accessor() > 1.01 || space.hHeight_accessor() / space2.hHeight_accessor() < 0.99) {
		cout << space.hHeight_accessor() / space2.hHeight_accessor() << "% hHeight divergence" << endl;
		counter++;
	}
	if (space.VzDist_accessor() / space2.VzDist_accessor() > 1.01 || space.VzDist_accessor()/space2.VzDist_accessor() < 0.99) {
		cout << space.VzDist_accessor() / space2.VzDist_accessor() << "% VzDist divergence" << endl;
		counter++;
	}
	if (space.zDist_accessor() / space2.zDist_accessor() > 1.01 || space.zDist_accessor() / space2.zDist_accessor() < 0.99) {
		cout << space.zDist_accessor() / space2.zDist_accessor() << "% zDist divergence" << endl;
		counter++;
	}
	if (space.chirp_accessor() / space2.chirp_accessor() > 1.01 || space.chirp_accessor() / space2.chirp_accessor() < 0.99) {
		cout << space.chirp_accessor() / space2.chirp_accessor() << "% chirp divergence" << endl;
		counter++;
	}
	if (space.b_accessor() / space2.b_accessor() > 1.01 || space.b_accessor() / space2.b_accessor() < 0.99) {
		cout << space.b_accessor() / space2.b_accessor() << "% b divergence" << endl;
		counter++;
	}
	if (counter == 0)
		cout << "phase space dimensions identical within 1%" << endl;
	else
		counter = 0;
}

//Conservation Checking - Emittence based
bool phase_space::longitudinal_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW) {//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = hDimensionW * distH;
	double consValue2 = hDimensionH * distW;
	if (consValue > 0.000499 && consValue < 0.000501  && consValue2 > 0.000499 && consValue2 < 0.000501) {//randomly decided range to account for data error
		return true;
	}
	else {
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;;
		return false;
	}
}

bool phase_space::transverse_area_conservation(double hDimensionW, double intDistH, double hDimensionH, double intDistW) {//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
	double consValue = hDimensionW * intDistH;
	double consValue2 = hDimensionH * intDistW;
	if (consValue > 0.00299 && consValue < 0.00301  && consValue2 > 0.00299 && consValue2 < 0.00301) {//randomly decided range to account for data error
		return true;
	}
	else {
		cout << "Area1: " << consValue << "  Area2: " << consValue2 << endl;
		return false;
	}
}

void outputPhaseSpace(ofstream& file, phase_space pulse, string name) {
	tuple<double, double, double, double, double, double> response = pulse.valid_variables_check();
	file <<
		"	<tr>" << endl <<
		"		<th> " << name.c_str() << " </th>" << endl <<
		"		<td> " << pulse.hWidth_accessor() << " </td>" << endl <<
		"		<td> " << pulse.hHeight_accessor() << " </td>" << endl <<
		"		<td> " << pulse.VzDist_accessor() << " </td>" << endl <<
		"		<td> " << pulse.zDist_accessor() << " </td>" << endl <<
		"		<td> " << pulse.chirp_accessor() << " </td>" << endl <<
		"		<td> " << pulse.b_accessor() << " </td>" << endl <<
		"		<td> " << pulse.longitudinal_area_conservation(pulse.hWidth_accessor(), pulse.VzDist_accessor(), pulse.hHeight_accessor(), pulse.zDist_accessor()) << " </th>" << endl <<
		"	</tr>" << endl <<
		"	<tr> " << endl <<
		"		<th> ValidityRatio </th>" << endl; 
	if (get<0>(response) = get<1>(response) = get<2>(response) = get<3>(response) = get<4>(response) = get<5>(response) = 1) {
		file << "		<td colspan=\"7\"> all variables valid </td>" << endl <<
						"	</tr>" << endl;
	}
	else {
		file << "		<td> " << get<0>(response) << " </td>" << endl <<
			"		<td> " << get<1>(response) << " </td>" << endl <<
			"		<td> " << get<2>(response) << " </td>" << endl <<
			"		<td> " << get<3>(response) << " </td>" << endl <<
			"		<td> " << get<4>(response) << " </td>" << endl <<
			"		<td> " << get<5>(response) << " </td>" << endl <<
			"	</tr>" << endl;
	}
}

void finalDataOutput() {
	if (printStarts)
		cout << "writing to final data file started" << endl;
	ofstream finalOutputFile;
	finalOutputFile.open("final_output.html");
	finalOutputFile <<
		"<!DOCTYPE html>" << endl <<
		"<html>" << endl <<
		"<head>" << endl <<
		"<title> Simulator data </title>" << endl <<
		"<link rel=\"stylesheet\" href=\"final_output_file.css\">" << endl <<
		"</head>" << endl <<
		"<body>" << endl <<
		"<h1> Simulator data </h1>" << endl <<
		"<table style=\"width=100%\"> " << endl <<
		"	<tr>" << endl <<
		"		<th> Name </th>" << endl <<
		"		<th> hWidth </th>" << endl <<
		"		<th> hHeight </th>" << endl <<
		"		<th> VzDist </th>" << endl <<
		"		<th> zDist </th>" << endl <<
		"		<th> chirp </th>" << endl <<
		"		<th> b </th> " << endl <<
		"		<th> longitudinal conservation </th> " << endl <<
		"	</tr>" << endl;
	if (printInitialPhaseSpace) {
		outputPhaseSpace(finalOutputFile, returnInitialPS(), "initialPulse");
	}
	if (printFinalPhaseSpace) {
		outputPhaseSpace(finalOutputFile, returnFinalPS(), "finalPulse");
	}
	finalOutputFile <<
		"</table>" << endl <<
		"<br> <br>" << endl <<
		"<table style = \"width=100%\">" << endl <<
		"	<tr>" << endl <<
		"		<td> " << "label1" << " </td>" << endl <<
		"		<td> " << "label2" << " </td>" << endl <<
		"		<td> " << "label3" << " </td>" << endl <<
		"		<td> " << "deviation" << " </td>" << endl <<
		"		<td> " << "label5" << " </td>" << endl <<
		"		<td> " << "label6" << " </td>" << endl <<
		"	</tr>" << endl <<
		"	<tr>" << endl <<
		"		<td> "<< valueHolder1 << " </td>" << endl <<
		"		<td> "<< valueHolder2 << " </td>" << endl <<
		"		<td> " << valueHolder3 << " </td>" << endl <<
		"		<td> " << valueHolder4 << " </td>" << endl <<
		"		<td> " << valueHolder5 << " </td>" << endl <<
		"		<td> " << valueHolder6 << " </td>" << endl <<
		"	</tr>" << endl <<
		"</table>" << endl <<
		"</body>" << endl <<
		"</html>" << endl;
	finalOutputFile.close();
	if (printEnds)
		cout << "writing to final data file finished" << endl;

	ShellExecute(NULL, "open", "final_output.html",
		NULL, NULL, SW_SHOWNORMAL);
}