import std.stdio, std.array, std.algorithm, std.conv, std.math;
class PhaseSpace{
	double sigma, zeta, pos, v, a, b; //Can sigma and zeta be streamlined?
	double[][] data;
    this(double v, double pos){
    	data ~= [v];
    	data ~= [pos];
	}
	void currentPhaseSpace(){
		writeln(text(data));
	}
		//Phase Space Evolution//
//1
	void goForward(double dt){ //just here temporarily
		foreach (ref items; data) { // overwrite value
    		foreach (ref item; items)  item *= dt;
		}
	}
//2
	void reverseChirpA(double change){//a is a value, not 1/-1 shown in equation 1 on the sheet
		this.a = a + change;
	}
	void reverseChirpB(double change){//Is there a way to make one method but able to change both a and b for the specific object?
		this.b = b + change;
	}
//Conservation Checking
	bool checkAreaConservation(double sigma, double zeta, double pos, double v, double a){	//Form A and B use different variables for sigma, zeta, and a.
		//double consAValue = exp(-pow(-a,2)/(2*(pow(zeta,2))))/(2*PI*pow(sigma*zeta,2));		//Need to figure out how to use pos and v (Z and Vz) this equation sets them to 0)
		double consAValue = exp((pow(zeta,2)/(-2*pow(sigma,2)))-pow(v-a,2)/(2*(pow(zeta,2))))/(2*PI*pow(sigma*zeta,2));	//Form A and B use different variables for sigma, zeta, and a. 
		if(consAValue > 0.9 && consAValue < 1.1){//How much range will need to be decided, in reality the value actually goes above 1!
			writeln("Area is conserved");
			return true;
		}
		else{
			writeln("Area is not conserved");
			writeln("Area:", consAValue);
			return false;
		}
	bool checkVolumeConservation(){//Two checkAreaConservations but with different variables for the 1st and 2nd axis.
		/*if(checkAreaConservation() &&checkAreaConservation()){
			writeln("Volume is conserved");*/
			return true;
		}
	}
}

void main(){
	auto space = new PhaseSpace(2.001, 5.009); //Random numbers
	space.currentPhaseSpace();
    space.goForward(4);
    space.currentPhaseSpace();
	space.checkAreaConservation(1,1,1,1,1);
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
}