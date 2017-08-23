import std.stdio, std.array, std.algorithm, std.conv, std.math, arsd.simpledisplay, script; //std.datetime;
/*		Math Conversion
sigmaZ = width;		the width along the z axis
sigmaVz = height;	the height along the Vz axis
zetaZ = VzIntDist;	the distance between the two Vz intercepts				//I think since it is between z (distance) and Vz(relative difference in velocity)
psi = zIntDist;		the distance between the two z intercepts				//Its important to write z and Vz, instead of x and y
a = chirp;			the chirp. chirp = slope
z = dist;			the distance from the origin-from the center of mass
Vz = vel;			the difference in velocity from center of mass
b = b;				the relationship between the chirp, zIntDist, and VzIntDist. (b=a(tau/zeta)^2)

Objectives
Finish opticalManipulation method
Fragment Phase Space
*/
double totalPulseEnergy, electronAmount, 
	distCT1, distT12, distT2L1, distL1A, distAT3, distT3S, distST4, distT4L2, distL2A, distAT5, distT5C;

class PhaseSpace{
	double width, height, VzIntDist, zIntDist, chirp, b;
	this(double widthC, double heightC, double VzIntDistC, double zIntDistC, double chirpC, double bC){
		width=widthC, height=heightC, VzIntDist=VzIntDistC, zIntDist=zIntDistC, chirp=chirpC, b=bC;
	}

	PhaseSpace freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		this.VzIntDist = sqrt(1/((1/pow(height,2))+pow((b/zIntDist),2)));
		writeln(VzIntDist);
		this.chirp = b*pow(VzIntDist,2)/pow(zIntDist,2);
		writeln(chirp);
		writeln(b);
		this.width = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
		writeln(width);
		return this;
	}
	PhaseSpace[] split(int spaces){
		PhaseSpace[] phaseSpaces;
		for(int i = 0; i<spaces; i++){
			phaseSpaces ~= new PhaseSpace(this.height/this.chirp, this.height/spaces*(i+1), this.VzIntDist, this.zIntDist, this.chirp, this.b);
		}
		return phaseSpaces;
	}
	PhaseSpace opticalManipulation(double changeChirp){
		this.chirp += changeChirp;
		return this;
	}

	PhaseSpace modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(6*width), to!int(6*height)); 
		{// introduce sub-scope
			auto painter = window.draw(); // begin drawing
			double x, y;
			x = -3*width, y = -3*height;
			while(y < 3*height){
				while(x < 3*width){
					double h = 1000000000000*exp((-1*pow(x,2)/(2*pow(width,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(width*VzIntDist,2));
					if(h>71){ //Value at F(300,260)
						painter.outlineColor = Color.white;
						painter.fillColor = Color.white;					
						if(h>3861){ //Value at F(100,87)
							painter.outlineColor = Color.blue;
							painter.fillColor = Color.blue;
						}
						if(h>5618){//Value at F(50,43)
							painter.outlineColor = Color.green;
							painter.fillColor = Color.green;
						}
						if(h>6170){//Value at F(25,21)
							painter.outlineColor = Color.yellow;
							painter.fillColor = Color.yellow;
						}
						if(h>6334){//Value at F(10,9)
							painter.outlineColor = Color.red;
							painter.fillColor = Color.red;
						}
						painter.drawLine(Point(to!int(x+(width*3)), to!int(-y+(height*3))), Point(to!int(x+(width*3)+1), to!int(-y+(height*3)+1)));
					}
					x += accuracy; //accuracy
				}
				x = -3*width;
				y += accuracy;
			}
		} // end scope, calling `painter`'s, drawing to the screen.
		window.eventLoop(0);// handle events
		return this;
	}

	//Conservation Checking
	bool checkAreaConservation(double widthHeight, double intDist){
		double consValue = widthHeight*intDist;
		if(consValue > 0.9999 && consValue < 1.0001){//randomly decided range to account for data error 
			writeln("Area is conserved");
			return true;
		}
		else{
			writeln("Area is not conserved");
			writeln("Area:", consValue);
			return false;
		}
	}
}


void main(){
	auto test = new Script("test.xml");
	test.run();
	/*auto initialPulse = new PhaseSpace(100,86.6,50,50,0.866,0.866);
	//auto initialPulse = new PhaseSpace(100, 87, 50, 50, 0.866, 0.866);
	writeln("Modeling Initial Pulse...");
	initialPulse.modelPhaseSpace(1);
	writeln("How long is free expansion?");
	string input = stdin.readln();
	initialPulse.freeExpansion(parse!double(input)).modelPhaseSpace(1);
	auto splitPhases = initialPulse.split(3);
	foreach (PhaseSpace space; splitPhases) {
    	space.modelPhaseSpace(1);
	}*/
	/*writeln("How much did the optical lens alter the chirp?");
	input = stdin.readln();
	initialPulse.opticalManipulation(parse!double(input)); Need to finish
	initialPulse.modelPhaseSpace(1);*/
	writeln("End of Program, enter anything to continue");
	input = stdin.readln();
}


/*Important things to know

We will need to create a phase space simulation that not only has the shape of the phase space
but and VERY IMPORTANTLY
has the data gathered from hitting the specimen in each point.

Objectives
Create a phase space simulator that simulates the phase space and collects data from hitting the specimen
in order to collect the data it needs to shatter
recombine the phase spaces at the end and collect the data

Time Table
Due September 19
Things to do
Finish free expansion
Figure out shattering
Recenter center of mass individually for each phase space
Figure out data containment
Figure out how to recombine

LARGE PART NOT YET FIGURED OUT Figure out how to optimize the microscope through combinations of optics
Write an at most 18 pages (not including references) paper about our research
Total 1600 projects
Score top 300 for semi finalist (top 18.75%)
Score top 60 for regional finalist (top 3.75%)
Score top 6 (in team category) for national finalist (top 0.375%) top (0.65% of team projects)
Score first (in team category) (top 0.0625%) top (0.125% of team projects)
*/