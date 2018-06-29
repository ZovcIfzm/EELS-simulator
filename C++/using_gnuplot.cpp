
#include <iostream>
#include <string>
#include "using_gnuplot.h"
#include "gnuplot_i.hpp"
using namespace std;

void gnuplot6(){
	Gnuplot g6("lines");
	cout << "window 6: splot with contour" << endl;
	g6.set_isosamples(60).set_contour();
	g6.unset_surface().plot_equation3d("sin(x)*sin(y)+4");
	g6.set_surface().replot();
	pause_method();
}

void pause_method(){
	string rando;
	cout << "Press enter to continue..." << endl;
	getline(cin, rando);
}