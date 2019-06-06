#include <iostream>
#include "data_processing.h"
#include <fstream>
#include <Windows.h>
#include "main_sequence.h";
#include <tuple>;

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
	if (get<0>(response) == get<1>(response) == get<2>(response) == get<3>(response) == get<4>(response) == get<5>(response) == 1) {
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
	finalOutputFile.open("finalOutput.html");
	finalOutputFile <<
		"<!DOCTYPE html>" << endl <<
		"<html>" << endl <<
		"<head>" << endl <<
		"<title> Simulator data </title>" << endl <<
		"<link rel=\"stylesheet\" href=\"finalOutputFile.css\">" << endl <<
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
		"</body>" << endl <<
		"</html>" << endl;
	finalOutputFile.close();
	if (printEnds)
		cout << "writing to final data file finished" << endl;

	ShellExecute(NULL, "open", "finalOutput.html",
		NULL, NULL, SW_SHOWNORMAL);
}