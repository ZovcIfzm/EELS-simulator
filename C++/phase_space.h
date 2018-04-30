class phase_space{
private:
	double hWidth, hHeight, VzDist, zDist, chirp, b, pulseEnergy, intensityRatio,
	hDepth, hDepthVel, VxDist, xDist, chirpT, bT, vZC, zC, xC;

public:
	phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityRatioC,
	double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC);

phase_space evolution(double time);
};