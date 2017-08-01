import std.stdio, std.array, std.algorithm, std.conv;
class PhaseSpace
{
	double data[][];
    this(double v, double pos){
    	data ~= [v];
    	data ~= [pos];
	}
	void currentPhaseSpace(){
		writeln(text(data));
	}
	void goForward(double dt){ //just here temporarily
		foreach (ref items; data) { // overwrite value
    		foreach (ref item; items)  item *= dt;
		}
	}
}
void main()
{
	auto space = new PhaseSpace(2.001, 5.009); //Random numbers
	space.currentPhaseSpace();
    space.goForward(4);
    space.currentPhaseSpace();
}