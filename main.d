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
double originalHWidth, originalHWidth2;
double[][] spectroTable;
double totalEnergy = 100E3, electronAmount = 100E3;
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}
double map( double x, double in_min, double in_max, double out_min, double out_max){
	if(x < in_min){
		x = in_min;
	}
	return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min; 
}
int count = 0;
double[100] pixel = 0;//!!Temporary!!
class PhaseSpace{
	double hWidth, hHeight, VzIntDist, zIntDist, chirp, b, totalPulseEnergy = 0, intensityRatio = 0;
	double hDepth, hDepthVelocity, VxIntDist, xIntDist, chirpT, bT;
	double vZC = 0, zC = 0, xC = 0;
	this(double hWidthC, double hHeightC, double VzIntDistC, double zIntDistC, double chirpC, double bC, double totalPulseEnergyC, double intensityRatioC,
		double hDepthC, double hDepthVelocityC, double VxIntDistC, double xIntDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC){
		hWidth=hWidthC, hHeight=hHeightC, VzIntDist=VzIntDistC, zIntDist=zIntDistC, chirp=chirpC, b=bC, totalPulseEnergy=totalPulseEnergyC, intensityRatio=intensityRatioC,
		hDepth=hDepthC, hDepthVelocity=hDepthVelocityC, VxIntDist=VxIntDistC, xIntDist=xIntDistC, chirpT=chirpTC, bT=bTC, vZC=vZCC, zC=zCC, xC=xCC;
	}
	this(PhaseSpace[] spaces){
		this.hWidth = originalHWidth + (spaces[0].hWidth-originalHWidth)*spaces.length;
		this.hHeight = spaces[0].hHeight * spaces.length; //If you add a * 1/chirp here sometimes it doesn't process it for some reason, a 1/chirp isn't needed here anyway but its an odd mystery why its only sometimes processed
		//this.vZC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.vZC"(spaces))*this.intensityRatio;
		foreach(PhaseSpace space; spaces){
			this.vZC += space.vZC*space.intensityRatio;
			this.zC += space.zC*space.intensityRatio;
			this.xC += space.xC*space.intensityRatio;
		}
		//this.zC = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.zC"(spaces));
		this.zIntDist = spaces[0].zIntDist;
		this.VzIntDist = zIntDist*hHeight/hWidth;
		this.chirp = VzIntDist*sqrt((1/pow(zIntDist,2))-(1/pow(hWidth,2)));
		this.b = this.chirp*pow(this.zIntDist/this.VzIntDist,2);
		this.hDepth = spaces[0].hDepth;
		this.hDepthVelocity = spaces[0].hDepthVelocity;
		this.VxIntDist = spaces[0].VxIntDist;
		this.xIntDist = spaces[0].xIntDist;
		this.chirpT = spaces[0].chirpT;
		this.bT = spaces[0].bT;
		this.totalPulseEnergy = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.totalPulseEnergy"(spaces));
		this.intensityRatio = taskPool.reduce!"a + b"(0.0, std.algorithm.map!"a.intensityRatio"(spaces));
	}

	PhaseSpace[] split(long spaces){
		originalHWidth = hWidth;
		PhaseSpace[] phaseSpaces;
		phaseSpaces.length = to!int(spaces);
		double spacesD = (to!double(spaces));
		double[] intensityRatios;
		intensityRatios.length = to!int(spaces);
		foreach (i, ref elem; parallel(intensityRatios)) {
			elem = getSplitIntensityRatio(spacesD, i+1, hHeight, hWidth);
		}
		VzIntDist = VzIntDist/spacesD;
		chirp = VzIntDist*sqrt((1/pow(zIntDist,2))-(1/pow(hWidth,2)));
		b = chirp*pow(zIntDist/VzIntDist,2);
		foreach (i, ref elem; phaseSpaces) {
			elem = new PhaseSpace(this.hWidth, this.hHeight/spacesD, this.VzIntDist, this.zIntDist, this.chirp, this.b, this.totalPulseEnergy*intensityRatios[i], intensityRatios[i],
										  this.hDepth, this.hDepthVelocity, this.VxIntDist, this.xIntDist, this.chirpT, this.bT, this.hHeight-(hHeight*2/spacesD)*(i+0.5), this.zC, this.xC);
		}
		count += spaces;
		return phaseSpaces;
	}
	PhaseSpace[] shatter(){
		originalHWidth = hWidth;
		double spaces = spectroTable[0].length;
		PhaseSpace[] phaseSpaces2;
		phaseSpaces2.length = to!int(spaces);
		VzIntDist = VzIntDist/spaces;
		chirp = VzIntDist*sqrt((1/pow(zIntDist,2))-(1/pow(hWidth,2)));
		b = chirp*pow(zIntDist/VzIntDist,2);
		foreach (i, ref elem; phaseSpaces2) {
			elem = new PhaseSpace(this.hWidth, this.hHeight/spaces, this.VzIntDist, this.zIntDist, this.chirp, this.b, this.totalPulseEnergy*spectroTable[1][i], this.intensityRatio*spectroTable[1][i],
								  this.hDepth, this.hDepthVelocity, this.VxIntDist, this.xIntDist, this.chirpT, this.bT, this.hHeight+(spectroTable[0][i]/1117), this.zC, this.xC);
		}
		count += spaces;
		return phaseSpaces2;
	}
	PhaseSpace modelPhaseSpace(double accuracy){
		auto window = new SimpleWindow(to!int(600), to!int(600)); 
		{// introduce sub-scope;
			auto painter = window.draw(); // begin drawing
			double x, y;
			x = -2*hWidth, y = -2*hHeight;
			while(y < 2*hHeight){
				while(x < 2*hWidth){
					double h = exp((-1*pow(x,2)/(2*pow(hWidth,2)))-(pow(y-chirp*x,2)/(2*pow(VzIntDist,2))))/(2*PI*pow(hWidth*VzIntDist,2));
					
					painter.outlineColor(Color(0, 0, to!int(map(h, 386090.0, 636622.0, 0, 255))));
					painter.drawLine(Point(to!int(x*150+300), to!int(-y*33300+300)), Point(to!int(x*150+300+1), to!int(-y*33300+300+1)));
					x += hWidth/300;
				}
				x = -2*hWidth;
				y += hHeight/300;
			}
		} // end scope, calling `painter`'s, drawing to the screen.
		window.eventLoop(0);// handle events
		return this;
	}
	double getSplitIntensityRatio(double numSections, double sectionNum, double hHeight, double hWidth){
		//Gets intensity % proportionally to 1 (like if its gets .5 its 50% of total intensity)
		//search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space (equal to 1)	
		double ySearchLB = -5.803*hHeight + ((5.803*hHeight*2.0/numSections)*(sectionNum-1));
		double ySearchUB = 5.803*hHeight - (5.803*hHeight*2.0/numSections)*(numSections - sectionNum);
		double xSearchLB = -5.803*hWidth;
		double xSearchUB = 5.803*hWidth;
		double accuracy = ySearchUB-ySearchLB;
		double x = xSearchLB;
		double y = ySearchLB;
		double intensityRatio = 0;
		double negTwohWidthsq = -2*hWidth*hWidth;
		double twoVzIntDistsq = 2*VzIntDist*VzIntDist;
		double twoPIhWidthsqVzIntDistsq = 2*PI*(hWidth*hWidth*VzIntDist*VzIntDist);
		if(numSections==sectionNum){
			ySearchUB += 0.000000000001;
		}
		while(y < ySearchUB-0.000000000001){
			while(x < xSearchUB){
				intensityRatio += 0.1*accuracy*exp((x*x/(negTwohWidthsq))-((y-chirp*x)*(y-chirp*x)/(twoVzIntDistsq)))/(twoPIhWidthsqVzIntDistsq);
				x += 0.1;
			}
			y += accuracy;
			x = xSearchLB;
		}
		return intensityRatio/2000;		
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
		this.chirp = b*pow(VzIntDist/zIntDist,2);
		this.chirpT = bT*pow(VxIntDist/xIntDist,2);
		this.hWidth = sqrt(1/((1/pow(zIntDist,2))-pow(chirp/VzIntDist,2)));
		this.hDepth = sqrt(1/((1/pow(xIntDist,2))-pow(chirpT/VxIntDist,2)));
		return this;
	}
	PhaseSpace RFLens(double changeChirp){
		this.chirp += changeChirp;
		this.zIntDist = sqrt(1/((1/pow(hWidth,2))+pow(chirp/VzIntDist,2)));
		this.b = this.chirp*pow(this.zIntDist/this.VzIntDist,2);
		this.hHeight = sqrt(1/((1/pow(VzIntDist,2))-pow(b/zIntDist,2)));
		return this;
	}
	PhaseSpace magLens(double changechirpT){
		this.chirpT += changechirpT;
		this.xIntDist = sqrt(1/((1/pow(hDepth,2))+pow(chirpT/VxIntDist,2)));
		this.bT = this.chirpT*pow(this.xIntDist/this.VxIntDist,2);
		this.hDepthVelocity = sqrt(1/((1/pow(VxIntDist,2))-pow(bT/xIntDist,2)));
		return this;
	}
	PhaseSpace spectroscopyFunction(){
		this.xC = this.xC + 7172.99042634*this.vZC;
		return this;
	}
	//Conservation Checking - Emittence based
	bool checkAreaConservationLongitudinal(double hDimensionW, double intDistH, double hDimensionH, double intDistW){//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
		double consValue = hDimensionW*intDistH;
		double consValue2 = hDimensionH*intDistW;
		if(consValue > 0.000499 && consValue < 0.000501  && consValue2 > 0.000499 && consValue2 < 0.000501){//randomly decided range to account for data error
			return true;
		}
		else{
			writeln("Area1: ", consValue, "  Area2: ", consValue2);
			return false;
		}
	}
	bool checkAreaConservationTransverse(double hDimensionW, double intDistH, double hDimensionH, double intDistW){//Needs to be reworked, both cons1&2 should describe the same emmittence- be the same value, however height and width are different but the intDist are the same
		double consValue = hDimensionW*intDistH;
		double consValue2 = hDimensionH*intDistW;
		if(consValue > 0.00299 && consValue < 0.00301  && consValue2 > 0.00299 && consValue2 < 0.00301){//randomly decided range to account for data error
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
		writeln("xC: ", xC);
		writeln("totalPulseEnergy: ", totalPulseEnergy);
		writeln("intensityRatio: ", intensityRatio);
		writeln("totalIntensityRatio from PixelSum: ", sumUp());
		writeln("");
		writeln("Longitudinal Emmittence Conserved: ", checkAreaConservationLongitudinal(hWidth, VzIntDist, hHeight, zIntDist));
		writeln("Transverse Emmittence Conserved: ", checkAreaConservationTransverse(hDepth, VxIntDist, hDepthVelocity, xIntDist));
		writeln("");
	}
	PhaseSpace specModeling()
	{ 
		auto window = new SimpleWindow(to!int(400), to!int(400));
		{// introduce sub-scope;
			auto painter = window.draw(); // begin drawing
		
			int count=0;
			double energy, intensity;
			foreach(double pix; pixel)
			{
				intensity = pixel[count];
				//energy = xC/6.421;
				energy = xC;
				painter.drawLine(Point(to!int(energy), 0), Point(to!int(energy), to!int(intensityRatio*electronAmount)));
				count++;
			}
		} // end scope, calling `painter`'s, drawing to the screen.
		window.eventLoop(0);// handle events
		return this;
	}
	void pixelSum(){
		pixel[to!int((xC+5000)/50)] += intensityRatio;//!!Temporary!!
	}
}
double sumUp(){
	double sumA=0;
	for(int i = 0; i<pixel[].length; i++){
		sumA += pixel[i];
		writeln(pixel[i]);
	}
	return sumA;
}

void main(){
	spectroTable = readSpec("hexogon BN-powder-eels.sl0");
	auto test = new Script("test.xml");
	test.run();
	/*PhaseSpace pulse = new PhaseSpace(1.5,0.0006218811,0.000333333,0.804012205,0.00035,2036.27222124,100,1,100,86.6,50,50,0.866,0.866,0,0);
	int i = 0;
	double counterB = 0;
	while(i<100000){
		counterB += pulse.getSplitIntensityRatio(100000, i,0.0006218811,1.5);
		i =i+1;
	}*/
	//writeln(pulse.getSplitIntensityRatio(100, 100,0.0006218811,1.5));
	//writeln(counterB);
	writeln("Total Fragmentated Phase Spaces: ", count);
	writeln("End of Program, enter anything to continue");
	string input = stdin.readln();
}