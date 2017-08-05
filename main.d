import std.stdio, std.array, std.algorithm, std.conv, std.math;
class PhaseSpace{
	double sigmaZ, zetaZ, pos, a;
	double sigmaVz, zetaVz, Vz, b, tau;
	double[][] data;
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
	void freeExpansion(){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.zetaZ = sqrt(1/((1/pow(sigmaVz,2))+pow((b/tau),2)));
		this.a = b*pow(zetaZ,2)/pow(tau,2);
		this.sigmaZ = sqrt(1/((1/pow(tau,2))-pow(a/zetaZ,2)));
	}
	void opticalManipulationA(double change){
		this.a = a + change;
	}
	void opticalMantipulationB(double change){//Is there a way to make one method but able to change both a and b for the specific object?
		this.b = b + change;
	}
/*Conservation Checking
	Form A and B use different variables for sigma, zeta, and a.
	consAValue needs to be integrated! 
	This process only needs to be run as a sort of debug tool for the programmer- the program might not necessarily need it
	-Should we have auto checkAreaConservation? 
	-(More processing would be needed but as different variables are inputed the user can be warned if the program doesn't work and then it can be fixed)
*/
	bool checkAreaConservation(double sigma, double zeta, double pos, double v, double a){
		double consAValue = exp((pow(zeta,2)/(-2*pow(sigma,2)))-pow(v-a,2)/(2*(pow(zeta,2))))/(2*PI*pow(sigma*zeta,2));	//needs to be integrated
		if(consAValue > 0.9 && consAValue < 1.1){//randomly decided range to account for integration error 
			writeln("Area is conserved");
			return true;
		}
		else{
			writeln("Area is not conserved");
			writeln("Area:", consAValue);
			return false;
		}
	bool checkVolumeConservation(){//Two checkAreaConservations but with different variables for form A and form B
		/*if(checkAreaConservation() &&checkAreaConservation()){
			writeln("Volume is conserved");*/
			return true;
		}
	}
}

void main(){
	auto space = new PhaseSpace(2.001, 5.009); //Random numbers
	space.currentPhaseSpace();
    //space.goForward(4);
    space.currentPhaseSpace();
	space.checkAreaConservation(1,1,1,1,1);
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
}