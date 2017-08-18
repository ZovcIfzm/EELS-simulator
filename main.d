import std.stdio, std.array, std.algorithm, std.conv, std.math, arsd.simpledisplay;//simpledisplay; //std.datetime;
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
	//Initial Conditions of initial pulse
	double width = 100;
	double height = 87;
	double VzIntDist = 50;
	double zIntDist = -50;
	double chirp = 0.866; //(sqrt 0.75)
	double b = 0.866;
	double dist, vel;

	void defineVariableValues(){
		writeln("Please define in order VzIntDist, zintDist, chirp, and b in order");
		VzIntDist = to!double(stdin.readln());
		writeln(VzIntDist);
		zIntDist = to!double(stdin.readln());
		writeln(zIntDist);
		chirp = to!double(stdin.readln());
		writeln(chirp);
		b = to!double(stdin.readln());
		writeln(b);
	}

	void freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		this.VzIntDist = sqrt(1/((1/pow(height,2))+pow((b/zIntDist),2)));
		writeln(VzIntDist);
		this.chirp = b*pow(VzIntDist,2)/pow(zIntDist,2);
		writeln(chirp);
		writeln(b);
		this.width = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
		writeln(width);
	}

	void opticalMantipulation(double changeChirp, double changeB){
		this.chirp += changeChirp;
		this.b += changeB;
	}

	void modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(6*width), to!int(6*height)); 
		{// introduce sub-scope
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
			auto x = -3*width;
			auto y = -3*height;
			painter.outlineColor = Color.yellow;
			painter.fillColor = Color.yellow;

			painter.drawLine(Point(to!int(60), to!int(60)), Point(to!int(50), to!int(50)));
			while(y < 3*height){
				while(x < 3*width){//Issues with scaling for some reason ( should have to be 7772 at 10^12, but is only 7772 at 10^14, further it doesn't reach values it should closer to center
					double h = 10000000000000*exp((-1*pow(x,2)/(2*pow(width,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(width*VzIntDist,2));
					if(h>7722){ //Value at F(100,87)
						painter.outlineColor = Color.blue;
						painter.fillColor = Color.blue;
						if(h>11236){
							painter.outlineColor = Color.green;
							painter.fillColor = Color.green;
						}
						if(h>12339){
							painter.outlineColor = Color.yellow;
							painter.fillColor = Color.yellow;
						}
						if(h>12668)
						{
							painter.outlineColor = Color.red;
							painter.fillColor = Color.red;
						}
						painter.drawLine(Point(to!int(x+300), to!int(-y+300)), Point(to!int(x+301), to!int(-y+301)));
					}
					x += accuracy; //accuracy
				}
				x = -3*width;
				y += accuracy;
			}
		} // end scope, calling `painter`'s destructor, drawing to the screen.
		window.eventLoop(0);// handle events
	}
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
	initialPulse.freeExpansion(1);
	initialPulse.modelPhaseSpace(1);
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