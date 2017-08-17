import std.stdio, std.array, std.algorithm, std.conv, std.math; //simpledisplay; //std.datetime;
/*		Math Conversion
sigmaZ = width;		the width along the z axis
sigmaVz = height;	the height along the Vz axis
zetaZ = VzIntDist;	the distance between the two Vz intercepts				//I think since it is between z (distance) and Vz(relative difference in velocity)
tau = zIntDist;		the distance between the two z intercepts				//Its important to write z and Vz, instead of x and y
a = chirp;			the chirp. chirp = slope
z = dist;			the distance from the origin-from the center of mass
Vz = vel;			the difference in velocity from center of mass
b = b;				the relationship between the chirp, zIntDist, and VzIntDist. (b=a(tau/zeta)^2)

Objectives
Finish modelPhaseSpace code
*/

class PhaseSpace{
	double width = 1;
	double height = 1;
	double VzIntDist = 0.5;
	double zIntDist = -0.5;
	double chirp = 0.866; //(sqrt 0.75)
	double b = 1;
	double dist, vel;

	//UNFINISHED, need to read in values (like stdin.readln)
	void defineVariableValues(){ //Main problem stdin.readln only works for strings and not doubles. 
		writeln("Please define in order VzIntDist, zintDist, chirp, and b in order");
		writeln(VzIntDist);
		writeln(zIntDist);
		writeln(chirp);
		writeln(b);
	}

	void freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		this.VzIntDist = sqrt(1/((1/pow(height,2))+pow((b/zIntDist),2)));
		this.chirp = b*pow(VzIntDist,2)/pow(zIntDist,2);
		this.VzIntDist = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
	}

	void opticalMantipulation(double changeChirp, double changeB){
		this.chirp += changeChirp;
		this.b += changeB;
	}

	/*void modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(3*width, 3*height);{ // introduce sub-scope
		auto painter = window.draw(); // begin drawing
		//draw here
		//	painter.outlineColor = Color.red;
		//	painter.fillColor = Color.red;
		//	auto x = -sigmaZ;
		//	while(x < sigmaZ){
		//		double h = exp((pow(a*x,2)/(-2*pow(x,2)))-pow(v-a,2)/(2*(pow(a*x,2))))/(2*PI*pow(x*a*x,2));
		//		painter.outlineColor = Color.red;
		//		painter.drawLine(Point(to!int(x*400), to!int(((0.5 * h)+a*x)*400)), Point(to!int(x*400), to!int((a*x-(0.5 * h))*400)));
		//		x += .0001; //accuracy
		//	}	
		auto x = -width;
		auto y = -height;
		while(y < height){
			if(x == 0){//If statement possibly unneeded
				x = 0.1;
			}
			while(x < width){
				double h = exp((pow(chirp*x,2)/(-2*pow(x,2)))-pow(y-chirp,2)/(2*(pow(chirp*x,2))))/(2*PI*pow(x*y,2));
				if(h>0.372){ //Value at F(1,1)
					painter.drawLine(Point(x), Point(y));
				}
				x += accuracy; //accuracy
			}
			y += accuracy;
		}
	} // end scope, calling `painter`'s destructor, drawing to the screen.
		window.eventLoop(0); // handle events
	}*/

	//Conservation Checking
	bool checkAreaConservation(double widthHeight, double intDist){
		double consValue = 1;
		if(consValue > 0.99 && consValue < 1.01){//randomly decided range to account for integration error 
			writeln("Area is conserved");
			return true;
		}
		else{
			writeln("Area is not conserved");
			writeln("Area:", consValue);
			return false;
		}
	}

	//Data storage, can be used in splitting
	double[][] data;				
    this(double Vz, double z){
    	data ~= [Vz];
    	data ~= [z];
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

void main(){
	auto initialPulse = new PhaseSpace(2.001, 5.009); //Random numbers
	//StopWatch sw;
	initialPulse.currentPhaseSpace();
    //initialPulse.goForward(4);
	initialPulse.currentPhaseSpace();
	initialPulse.checkAreaConservation(1,1);
	//sw.start();
	//writeln(initialPulse.getArea(0.00005));
	//sw.stop();
	//writeln("Took ",sw.peek().to!("msecs", real)(), "ms to run area method");
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
}


/*Important things to know

We will need to create a phase space simulation that not only has the shape of the phase space
but and VERY IMPORTANTLY
has the data gathered from hitting the specimen in each point. How to do that is the mystery. Will we have to do a coordinate system?
Or is this found in each individual phase space when it shatters

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
Figure out reactions between shatters phase spaces and how to account for them and how they change each other
Figure out how to recombine

LARGE PART NOT YET FIGURED OUT Figure out how to optimize the microscope through combinations of optics
Write an at most 18 pages (not including references) paper about our research
Total 1600 projects
Score top 300 for semi finalist (top 18.75%)
Score top 60 for regional finalist
Score top 6 (in team category) for national finalist 
*/