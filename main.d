import std.stdio, std.array, std.algorithm, std.conv, std.math, arsd.simpledisplay, script; //std.datetime;
/*		Math Conversion
sigmaZ = width;		the width along the z axis
sigmaVz = height;	the height along the Vz axis
zetaZ = VzIntDist;	the distance between the two Vz intercepts				//I think since it is between z (distance) and Vz(relative difference in velocity)
tau = zIntDist;		the distance between the two z intercepts				//Its important to write z and Vz, instead of x and y
a = chirp;			the chirp. chirp = slope
b = b;				the relationship between the chirp, zIntDist, and VzIntDist. (b=a(tau/zeta)^2)
z = dist;			the distance from the origin-from the center of mass
Vz = vel;			the difference in velocity from center of mass

Objectives
Finish opticalManipulation method
Fragment Phase Space
*/
double totalEnergy = 100, electronAmount, 
	distCT1, distT12, distT2L1, distL1A, distAT3, distT3S, distST4, distT4L2, distL2A, distAT5, distT5C;

class PhaseSpace{
	double width, height, VzIntDist, zIntDist, chirp, b, totalPulseEnergy = 0, intensityRatio = 0;
	//double emmittence - ability to find other varaibles in terms of these three
	//double emmittenceD - same but for depth
	double depth, depthVelocity, VxIntDist, xIntDist, chirpD, bD;
	this(double widthC, double heightC, double VzIntDistC, double zIntDistC, double chirpC, double bC, double totalPulseEnergyC, double intensityRatioC){
		// double depthC, double depthVelocityC, double VxIntDistC, double xIntDistC, double chirpDC, double bDC){
		width=widthC, height=heightC, VzIntDist=VzIntDistC, zIntDist=zIntDistC, chirp=chirpC, b=bC, totalPulseEnergy=totalPulseEnergyC, intensityRatio=intensityRatioC;
		//depth=depthC, depthVelocity=depthVelocityC, VxIntDist=VxIntDistC, xIntDist=xIntDistC, chirpD=chirpDC, bD=bDC;
	}
	this(PhaseSpace[] spaces){
		this.VzIntDist = spaces[0].VzIntDist;
		this.zIntDist = spaces[0].zIntDist;
		this.chirp = spaces[0].chirp;
		this.b = spaces[0].b;
		foreach(PhaseSpace space; spaces){
			this.totalPulseEnergy += space.totalPulseEnergy;
			this.intensityRatio += space.intensityRatio;
		}
		this.width = spaces[0].width * spaces.length;
		this.height = spaces[0].height * spaces.length * (1/chirp);

	}
	void printPhaseSpace(){
		writeln("width: ", width);
		writeln("height: ", height);
		writeln("VzIntDist: ", VzIntDist);
		writeln("zIntDist: ", zIntDist);
		writeln("chirp: ", chirp);
		writeln("b: ", b);
		writeln("totalPulseEnergy: ", totalPulseEnergy);
		writeln("intensityRatio: ", width);
	}
	PhaseSpace freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		if(chirp>0){
			this.VzIntDist = sqrt(1/((1/pow(height,2))+pow((b/zIntDist),2)));
			writeln(VzIntDist);
		}
		if(chirp<0){
			this.zIntDist = b/(sqrt((1/pow(VzIntDist,2))-(1/pow(height,2))));
			writeln(zIntDist);
		}
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
			double intensityRatio = getSplitIntensityRatio(1,spaces,i+1, this.height, this.width);
			phaseSpaces ~= new PhaseSpace((this.height/this.chirp)/spaces, this.height/spaces, this.VzIntDist, this.zIntDist, this.chirp, this.b, this.totalPulseEnergy*intensityRatio, intensityRatio);
		}
		return phaseSpaces;
	}
	PhaseSpace opticalManipulation(double changeChirp){
		this.chirp += changeChirp;
		return this;
	}

	PhaseSpace modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(6*width), to!int(6*height)); 
		{// introduce sub-scope;
			auto painter = window.draw(); // begin drawing
			double x, y;
			x = -3*width, y = -3*height;
			while(y < 3*height){
				while(x < 3*width){
					double h = this.intensityRatio*1000000000000*exp((-1*pow(x,2)/(2*pow(width,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(width*VzIntDist,2));
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

	double getSplitIntensityRatio(double accuracy, int numSections, double sectionNum, double height, double width){
		//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
		//search with xSearch & ySearch = +- 5.803*width or height to get the total intensity of the phase space (equal to 1)
		double xSearchLB; double xSearchUB; 
		double ySearchLB; double ySearchUB;	
		//xSearchLB = -width + (width*2/numSections)*(sectionNum-1);
		//xSearchUB = width - (width*2/numSections)*(numSections - sectionNum);
		ySearchLB = -5.803*height + ((5.803*height*2/numSections)*(sectionNum-1));
		ySearchUB = 5.803*height - (5.803*height*2/numSections)*(numSections - sectionNum);
		xSearchLB = -5.803*width;
		xSearchUB = 5.803*width;
		//ySearchLB = -5.803*height;
		//ySearchUB = 5.803*height;

		auto x = xSearchLB;
		auto y = ySearchLB;
		//double x = -width;
		double intensityRatio = 0;
		if(numSections==sectionNum){
			ySearchUB += 1;
		}
		while(y < ySearchUB-1){
			while(x < xSearchUB){
				intensityRatio += pow(accuracy,2)*exp((-1*pow(x,2)/(2*pow(width,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(width*VzIntDist,2));
				x += accuracy;
			}
			y += accuracy;
			x = xSearchLB;
		}
		intensityRatio = intensityRatio*5000;
		return intensityRatio;		
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
	/*auto initialPulse = new PhaseSpace(100,86.6,50,50,0.866,0.866, 100, 1);
	initialPulse.printPhaseSpace();
	double finInfo = 0;
	int numSectionsC = 1005;//Max amount before it stops working (cant do above 1005)
	for(int j=1;j<numSectionsC+1;j++){
		writeln((initialPulse.getSplitIntensityRatio(1,numSectionsC,j,initialPulse.height, initialPulse.width)));
		finInfo += (initialPulse.getSplitIntensityRatio(1,numSectionsC,j,initialPulse.height, initialPulse.width));
	}
	writeln((initialPulse.getSplitIntensityRatio(1,numSectionsC,1,initialPulse.height, initialPulse.width)));
	writeln(finInfo);*/
	/*writeln("Modeling Initial Pulse...");
	initialPulse.modelPhaseSpace(1);
	writeln("How long is free expansion?");
	string input = stdin.readln();
	initialPulse.freeExpansion(parse!double(input)).modelPhaseSpace(1);
	auto splitPhases = initialPulse.split(3);
	foreach (PhaseSpace space; splitPhases) {
    	space.modelPhaseSpace(1);
		writeln(space.intensityRatio);
	}
	new PhaseSpace(splitPhases).printPhaseSpace();
	writeln("How much did the optical lens alter the chirp?");
	input = stdin.readln();
	initialPulse.opticalManipulation(parse!double(input)); Need to finish
	initialPulse.modelPhaseSpace(1);*/
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
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