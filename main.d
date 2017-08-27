import std.stdio, std.array, std.algorithm, std.conv, std.math, std.parallelism, arsd.simpledisplay, script; //std.datetime;
/*		Math Conversion
sigmaZ = hWidth;		the hWidth along the z axis
sigmaVz = hHeight;	the hHeight along the Vz axis
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
int count = 0;
class PhaseSpace{
	double hWidth, hHeight, VzIntDist, zIntDist, chirp, b, totalPulseEnergy = 0, intensityRatio = 0;
	//double emmittence - ability to find other varaibles in terms of these three
	//double emmittenceD - same but for hDepth
	double hDepth, hDepthVelocity, VxIntDist, xIntDist, chirpD, bD;
	this(double hWidthC, double hHeightC, double VzIntDistC, double zIntDistC, double chirpC, double bC, double totalPulseEnergyC, double intensityRatioC){
		// double hDepthC, double hDepthVelocityC, double VxIntDistC, double xIntDistC, double chirpDC, double bDC){
		hWidth=hWidthC, hHeight=hHeightC, VzIntDist=VzIntDistC, zIntDist=zIntDistC, chirp=chirpC, b=bC, totalPulseEnergy=totalPulseEnergyC, intensityRatio=intensityRatioC;
		//hDepth=hDepthC, hDepthVelocity=hDepthVelocityC, VxIntDist=VxIntDistC, xIntDist=xIntDistC, chirpD=chirpDC, bD=bDC;
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
		this.hWidth = spaces[0].hWidth * spaces.length;
		this.hHeight = spaces[0].hHeight * spaces.length * (1/chirp);

	}
	void printPhaseSpace(){
		writeln("hWidth: ", hWidth);
		writeln("hHeight: ", hHeight);
		writeln("VzIntDist: ", VzIntDist);
		writeln("zIntDist: ", zIntDist);
		writeln("chirp: ", chirp);
		writeln("b: ", b);
		writeln("totalPulseEnergy: ", totalPulseEnergy);
		writeln("intensityRatio: ", intensityRatio);
		writeln("");
	}
	PhaseSpace freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		if(chirp>0){
			this.VzIntDist = sqrt(1/((1/pow(hHeight,2))+pow((b/zIntDist),2)));
			writeln(VzIntDist);
		}
		if(chirp<0){
			this.zIntDist = b/(sqrt((1/pow(VzIntDist,2))-(1/pow(hHeight,2))));
			writeln(zIntDist);
		}
		this.chirp = b*pow(VzIntDist,2)/pow(zIntDist,2);
		writeln(chirp);
		writeln(b);
		this.hWidth = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
		writeln(hWidth);
		return this;
	}
	PhaseSpace[] split(int spaces){
		PhaseSpace[] phaseSpaces;
		double spacesD = (to!double(spaces));
		foreach (i; taskPool.parallel(new int[spaces])) {
			double intensityRatio = getSplitIntensityRatio(1005/spacesD, spaces, i, this.hHeight, this.hWidth);
			phaseSpaces ~= new PhaseSpace((this.hHeight/this.chirp)/spacesD, this.hHeight/spacesD, this.VzIntDist, this.zIntDist, this.chirp, this.b, this.totalPulseEnergy*intensityRatio, intensityRatio);
		}
		count = spaces;
		return phaseSpaces;
	}
	PhaseSpace opticalManipulation(double changeChirp){
		this.chirp += changeChirp;
		return this;
	}
	PhaseSpace modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(6*hWidth), to!int(6*hHeight)); 
		{// introduce sub-scope;
			auto painter = window.draw(); // begin drawing
			double x, y;
			x = -3*hWidth, y = -3*hHeight;
			while(y < 3*hHeight){
				while(x < 3*hWidth){
					double h = this.intensityRatio*1000000000000*exp((-1*pow(x,2)/(2*pow(hWidth,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(hWidth*VzIntDist,2));
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
						painter.drawLine(Point(to!int(x+(hWidth*3)), to!int(-y+(hHeight*3))), Point(to!int(x+(hWidth*3)+1), to!int(-y+(hHeight*3)+1)));
					}
					x += accuracy; //accuracy
				}
				x = -3*hWidth;
				y += accuracy;
			}
		} // end scope, calling `painter`'s, drawing to the screen.
		window.eventLoop(0);// handle events
		return this;
	}
	double getSplitIntensityRatio(double accuracy, double numSections, double sectionNum, double hHeight, double hWidth){
		//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
		//search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space (equal to 1)	
		double ySearchLB = -5.803*hHeight + ((5.803*hHeight*2/numSections)*(sectionNum-1));
		double ySearchUB = 5.803*hHeight - (5.803*hHeight*2/numSections)*(numSections - sectionNum);
		double xSearchLB = -5.803*hWidth;
		double xSearchUB = 5.803*hWidth;
		double x = xSearchLB;
		double y = ySearchLB;
		double intensityRatio = 0;
		if(numSections==sectionNum){
			ySearchUB += 1;
		}
		while(y < ySearchUB-0.0001){
			while(x < xSearchUB){
				intensityRatio += pow(accuracy,2)*exp((-1*pow(x,2)/(2*pow(hWidth,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(hWidth*VzIntDist,2));
				x += accuracy;
			}
			y += accuracy;
			x = xSearchLB;
		}
		return intensityRatio*5000;		
	}
	//Conservation Checking
	bool checkAreaConservation(double hWidthHeight, double intDist){
		double consValue = hWidthHeight*intDist;
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
	writeln("Total fragmentated phase spaces",count);//For testing purposes
	/*auto initialPulse = new PhaseSpace(100,86.6,50,50,0.866,0.866, 100, 1);
	initialPulse.printPhaseSpace();
	double finInfo = 0;
	int numSectionsC = 1005;//Max amount before it stops working (cant do above 1005)
	for(int j = 1;j<numSectionsC+1; j++){
	//	writeln((initialPulse.getSplitIntensityRatio(1,numSectionsC,j,initialPulse.hHeight, initialPulse.hWidth)));
	//	writeln(j);
		finInfo += (initialPulse.getSplitIntensityRatio(1,numSectionsC,j,initialPulse.hHeight, initialPulse.hWidth));
	}
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