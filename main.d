import std.stdio, std.array, std.algorithm, std.conv, std.math, std.parallelism, std.range, std.string, std.file, arsd.simpledisplay, script;
double[][] readSpec(string filename){
	double[][] toReturn;
	toReturn.length = 2;
	auto lines = readText(filename).splitLines();
	foreach(string line; lines){
			string[] nums = line.split(",");
			if(isNumeric(nums[0].strip())){
				toReturn[0] ~= to!double(nums[0].strip());
				toReturn[1] ~= to!double(nums[1].strip());
			}
	}
	double sum = 0;
	foreach(double dat; toReturn[1])
		sum += dat;
	foreach (i, ref elem; toReturn[1]) {
		toReturn[1][i] = elem/sum;
	}
	return toReturn;
}
double[][] spectroTable;
double totalEnergy = 100E3, electronAmount = 100E3;
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}
int count = 0;
class PhaseSpace{
	double hWidth, hHeight, VzIntDist, zIntDist, chirp, b, totalPulseEnergy = 0, intensityRatio = 0;
	double hDepth, hDepthVelocity, VxIntDist, xIntDist, chirpT, bT;
	double vZC = 0, zC = 0;
	this(double hWidthC, double hHeightC, double VzIntDistC, double zIntDistC, double chirpC, double bC, double totalPulseEnergyC, double intensityRatioC,
		double hDepthC, double hDepthVelocityC, double VxIntDistC, double xIntDistC, double chirpTC, double bTC, double vZCC, double zCC){
		hWidth=hWidthC, hHeight=hHeightC, VzIntDist=VzIntDistC, zIntDist=zIntDistC, chirp=chirpC, b=bC, totalPulseEnergy=totalPulseEnergyC, intensityRatio=intensityRatioC,
		hDepth=hDepthC, hDepthVelocity=hDepthVelocityC, VxIntDist=VxIntDistC, xIntDist=xIntDistC, chirpT=chirpTC, bT=bTC, vZC=vZCC, zC=zCC;
	}
	this(PhaseSpace[] spaces){
		this.VzIntDist = spaces[0].VzIntDist;
		this.zIntDist = spaces[0].zIntDist;
		this.chirp = spaces[0].chirp;
		this.b = spaces[0].b;
		this.hDepth = spaces[0].hDepth;
		this.hDepthVelocity = spaces[0].hDepthVelocity;
		this.VxIntDist = spaces[0].VxIntDist;
		this.xIntDist = spaces[0].xIntDist;
		this.chirpT = spaces[0].chirpT;
		this.bT = spaces[0].bT;
		this.totalPulseEnergy = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.totalPulseEnergy"(spaces));
		this.intensityRatio = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.intensityRatio"(spaces));
		this.hWidth = spaces[0].hWidth * spaces.length;
		this.hHeight = spaces[0].hHeight * spaces.length; //If you add a * 1/chirp here sometimes it doesn't process it for some reason, a 1/chirp isn't needed here anyway but its an odd mystery why its only sometimes processed
		this.vZC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.vZC"(spaces));
		this.zC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.zC"(spaces));
	}
	PhaseSpace[] split(long spaces){
		PhaseSpace[] phaseSpaces;
		phaseSpaces.length = to!int(spaces);
		double spacesD = 1/(to!double(spaces));
		double[] intensityRatios;
		intensityRatios.length = to!int(spaces);
		foreach (i, ref elem; parallel(intensityRatios)) {
			elem = getSplitintensityRatio(1005*spacesD, spaces, i+1, this.hHeight, this.hWidth);
		}
		foreach (i, ref elem; phaseSpaces) {
			elem = new PhaseSpace((this.hHeight/this.chirp)*spacesD, this.hHeight*spacesD, this.VzIntDist, this.zIntDist, this.chirp, this.b, this.totalPulseEnergy*intensityRatios[i], intensityRatios[i],
										  this.hDepth, this.hDepthVelocity, this.VxIntDist, this.xIntDist, this.chirpT, this.bT, this.hHeight-(hHeight*2/spaces)*(i+0.5), this.hWidth-(hWidth*2/spaces)*(i+0.5));
		}
		count += spaces;
		return phaseSpaces;
	}
	PhaseSpace modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(6*hWidth), to!int(6*hHeight)); 
		{// introduce sub-scope;
			auto painter = window.draw(); // begin drawing
			double x, y;
			x = -3*hWidth, y = -3*hHeight;
			while(y < 3*hHeight){
				while(x < 3*hWidth){
					double map( double x, double in_min,double in_max, double out_min, double out_max)
					{
						return (100*(x - in_min) / (in_max - in_min) * (in_max - in_min) / (out_max - out_min) + out_min); 
					}
					double h = this.intensityRatio*1E12*exp((-1*pow(x,2)/(2*pow(hWidth,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(hWidth*VzIntDist,2));
					double a, b, c, d;
					a = 3861.0;
					b = 6366.2;
					c = 0.0;
					d = 100.0;
					painter.outlineColor(Color( to!int(map(h, a, b, c, d)), 0, 0, 255));
					painter.drawLine(Point(to!int(x+(hWidth*3)), to!int(-y+(hHeight*3))), Point(to!int(x+(hWidth*3)+1), to!int(-y+(hHeight*3)+1)));
				
					x += accuracy;
				}
				x = -3*hWidth;
				y += accuracy;
			}
		} // end scope, calling `painter`'s, drawing to the screen.
		window.eventLoop(0);// handle events
		return this;
	}
	double getSplitintensityRatio(double accuracy, double numSections, double sectionNum, double hHeight, double hWidth){
		//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
		//search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space (equal to 1)	
		double ySearchLB = -5.803*hHeight + ((5.803*hHeight*2/numSections)*(sectionNum-1));
		double ySearchUB = 5.803*hHeight - (5.803*hHeight*2/numSections)*(numSections - sectionNum);
		double xSearchLB = -5.803*hWidth;
		double xSearchUB = 5.803*hWidth;
		double x = xSearchLB;
		double y = ySearchLB;
		double intensityRatio = 0;
		double negTwohWidthsq = -2*hWidth*hWidth;
		double twoVzIntDistsq = 2*VzIntDist*VzIntDist;
		double twoPIhWidthsqVzIntDistsq = 2*PI*(hWidth*hWidth*VzIntDist*VzIntDist);
		if(numSections==sectionNum){
			ySearchUB += 0.0001;
		}
		while(y < ySearchUB-0.0001){
			while(x < xSearchUB){
				intensityRatio += 112*accuracy*exp1((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
				x += 112;
			}
			y += accuracy;
			x = xSearchLB;
		}
		return intensityRatio*5000;		
	}
	PhaseSpace freeExpansion(double time){//To deal with processing we might need to make our own math functions. (less/more digits of accuracy)
		this.b += time;
		this.bT += time;
		if(chirp>0){
			this.VzIntDist = sqrt(1/((1/pow(hHeight,2))+pow((b/zIntDist),2)));
		}
		if(chirpT>0){
			this.VxIntDist = sqrt(1/((1/pow(hDepthVelocity,2))+pow((bT/xIntDist),2)));
		}
		if(chirp<0){
			this.zIntDist = b/(sqrt((1/pow(VzIntDist,2))-(1/pow(hHeight,2))));
		}
		if(chirpT<0){
			this.xIntDist = bT/(sqrt((1/pow(VxIntDist,2))-(1/pow(hDepthVelocity,2))));
		}
		this.chirp = b*pow(VzIntDist,2)/pow(zIntDist,2);
		this.chirpT = bT*pow(VxIntDist,2)/pow(xIntDist,2);
		this.hWidth = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
		this.hDepth = sqrt(1/((1/pow(xIntDist,2))-pow(chirpT/VxIntDist,2)));
		return this;
	}
	PhaseSpace RFLens(double changeChirp){
		this.chirp += changeChirp;
		return this;
	}
	PhaseSpace magLens(double changechirpT){
		this.chirpT += changechirpT;
		return this;
	}
	//Conservation Checking - Emittence based
	bool checkAreaConservation(double hDimensionW, double intDistH, double hDimensionH, double intDistW){//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
		double consValue = hDimensionW*intDistH;
		double consValue2 = hDimensionH*intDistW;
		if(consValue > 4329 && consValue < 4331  && consValue2 > 4329 && consValue2 < 4331){//randomly decided range to account for data error
			return true;
		}
		else{
			writeln("Area1: ", consValue, "  Area2: ", consValue2);
			return false;
		}
	}
	void printPhaseSpace(){
		writeln("hWidth: ", hWidth, "  ",    " hDepth: ", hDepth);
		writeln("hHeight: ", hHeight, " ",   " hDepthVelocity: ", hDepthVelocity);
		writeln("VzIntDist: ", VzIntDist,    " VxIntDist: ", VxIntDist);
		writeln("zIntDist: ", zIntDist, " ", " xIntDist: ", xIntDist);
		writeln("chirp: ", chirp, " ",       " chirpT: ", chirpT);
		writeln("b: ", b, "     ",           " bT: ", bT);
		writeln("VzC: ", vZC, "   ",         "zC: ", zC);
		writeln("totalPulseEnergy: ", totalPulseEnergy);
		writeln("intensityRatio: ", intensityRatio);
		writeln("");
		writeln("Longitudinal Emmittence Conserved: ", checkAreaConservation(hWidth, VzIntDist, hHeight, zIntDist));
		writeln("Transverse Emmittence Conserved: ", checkAreaConservation(hDepth, VxIntDist, hDepthVelocity, xIntDist));
		writeln("");
	}
}
void electronLossCalculations(PhaseSpace[][] spaces){
	for(int i = 0;i<spaces.length;i++){
		for(int j = 0;spaces[i].length;j++){
			spaces[i][j].totalPulseEnergy -= spectroTable[1][j];
		}
	}
}
void main(){
	spectroTable = readSpec("hexogon BN-powder-eels.sl0");
	auto test = new Script("test.xml");
	test.run();
	writeln("Total Fragmentated Phase Spaces: ", count);
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
}