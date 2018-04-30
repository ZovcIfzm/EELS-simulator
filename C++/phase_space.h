class phase_space{
private:
	double hWidth, hHeight, VzDist, zDist, chirp, b, pulseEnergy, intensityRatio,
	hDepth, hDepthVel, VxDist, xDist, chirpT, bT, vZC, zC, xC;

public:
	phase_space(double hWidthC, double hHeightC, double VzDistC, double zDistC, double chirpC, double bC, double pulseEnergyC, double intensityRatioC,
	double hDepthC, double hDepthVelC, double VxDistC, double xDistC, double chirpTC, double bTC, double vZCC, double zCC, double xCC);

phase_space evolution(double time);
phase_space RFLens(double changeChirp);
phase_space mag_lens(double changeChirpT);
phase_space spectroscopy_function();
bool longitudinal_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
bool transverse_area_conservation(double hDimensionW, double distH, double hDimensionH, double distW);
void print();

};