/*
Title: ionTraj_v3.0.cpp
Project: FY20 MSFC CIF/TIP Direct Thrust Measurements of Electrostatic Sail Tethers
Description: A test particle code that can propogate electrons and ions in an
electric field.
Author: Anthony M. DeStefano
Company: NASA/MSFC/EV44
E-mail: anthony.m.destefano@nasa.gov
Office phone: 256-544-3094
Date last edited: 11/24/2019

Author: Josh Topliss
Company: NASA/MSFC/EV44 Spring Intern
E-mail: joshua.topliss@gmail.com
Cell phone: 256-560-5894
Date last edited: 2/13/2020
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
//#include <time.h>
#include <algorithm> // max(a,b)
//#include <mpi.h>

using namespace std;

#define DEBUG 0 // 0 = no output, 1 = per track, 2 = per time step
#define VERBOSE_INIT 0 // 0 = no output, 1 = all but tracks, 2 = all
#define HELP_INIT 0
#define PRINT_TRACK 0 // 0 = no output, 1 = all track info, 2 = only at specified x-location
#define PRINT_FORCE 0 // 0 = no output, 1 = integrated tracks, 2 = per track

// SI units
const double MASS_ELECTRON      = 9.10938356e-31;   // kg
const double MASS_PROTON        = 1.6726219e-27;    // kg
const double CHARGE_PROTON      = 1.6021766208e-19; // C
const double BOLTZMANN_CONSTANT = 1.38064852e-23;   // J/K
const double EPSILON0           = 8.8541878128e-12; // F/m

// useful functions
inline double sqr(double x) { return x*x; }
inline double mag(double x, double y) { return sqrt(x*x + y*y); }
inline double mag2(double x, double y) { return x*x + y*y; }
/* Used to compute index of a triangular matrix flattened into a 1-D array. */
inline int triIdx(int n, int m) { return ((n * (n+1)) >> 1) + m; }

// unit conversions
inline double eV_to_Joules(double E) { return E*CHARGE_PROTON; }
inline double Joules_to_eV(double E) { return E/CHARGE_PROTON; }
inline double eV_to_Kelvin(double T) { return eV_to_Joules(T)/BOLTZMANN_CONSTANT; }
inline double Kelvin_to_eV(double T) { return Joules_to_eV(T*BOLTZMANN_CONSTANT); }
inline double eV_to_mps(double V, double m) { return sqrt(2.0 * eV_to_Joules(V) / m); }
inline double mps_to_eV(double V, double m) { return Joules_to_eV(0.5 * m * sqr(V)); }
inline double eVpm2_to_density(double d) { return eV_to_Joules(d) * EPSILON0 / sqr(CHARGE_PROTON); }
inline double density_to_eVpm2(double d) { return Joules_to_eV(d) * sqr(CHARGE_PROTON) / EPSILON0; }


// Primative objects ///////////////////////////////////////

struct tether
{
	double potential;
	double radius;
	int numberOfTethers;
	double tetherSeparation;
	void (*EfieldFunction)(tether&, double&, double&, double, double, double);
}; //(tether& t, double& Ex, double& Ey, double x, double y, double r0)

void init_tether(tether& t, double p, double r, int nt, double ts, void (*eff)(tether&, double&, double&, double, double, double))
{
	t.potential = p;
	t.radius = r;
	t.numberOfTethers = nt;
	t.tetherSeparation = ts;
	t.EfieldFunction = eff;

	if(VERBOSE_INIT){
		cout << "Initializing tether paramerters:\n";
		cout << "  potential        = " << p << " V\n";
		cout << "  radius           = " << r << " m\n";
		cout << "  numberOfTethers  = " << nt << " \n";
		cout << "  tetherSeparation = " << ts << " m\n\n";
	}
	if(HELP_INIT){
		cout << "HELP: init_tether\n";
		cout << "  potential        -> The voltage potential of the tether wire. Assumed to be constant.\n";
		cout << "  radius           -> The radius of the tether wire. Assumed to be cylindrical.\n";
		cout << "  numberOfTethers  -> The number of parallel tethers that are in a row with a constant separation.\n";
		cout << "  tetherSeparation -> The separation distance between the tethers.\n";
		cout << "  EfieldFunction   -> The electric field magnitude as a function of position from the tether wire center.\n";
		cout << "                      Used to compute forces on charged particles.\n\n";
	}
}



// E-field of wire in vacuum (from Jackson)
void Efield_WireVacuum(tether& t, double& Ex, double& Ey, double x, double y, double r0)
{
	double r = mag(x, y);
	double Er = -t.potential/(r * log(r0/t.radius));

	Ex = Er * x / r;
	Ey = Er * y / r;
}


// E-field of wire in plasa (from Janhunen & Sandroos 2007)
void Efield_WirePlasma(tether& t, double& Ex, double& Ey, double x, double y, double r0)
{
	double r = mag(x, y);
	double Er = -2.0 * t.potential * sqr(r0) / (pow(r,3) * (1.0 + sqr(r0/r)) * log(1.0 + sqr(r0/t.radius)));

	Ex = Er * x / r;
	Ey = Er * y / r; 
}


// E-field of multiple (20) parallel wires in a plasma, separated by 2 cm
// const double MultiWireSep = 1.0e-3; // m 2e-2
// const int    MultiWireNum = 20;
void Efield_MultiWiresPlasma(tether& t, double& Ex, double& Ey, double x, double y, double r0)
{
	double dy = t.tetherSeparation;
	double ymin = -dy * (t.numberOfTethers - 1.0) / 2.0;
	double Er, xWire, yWire, xcur, ycur, rcur;
	double Exi = 0.0, Eyi = 0.0;
	Ex = 0.0;
	Ey = 0.0;

	for (int i = 0; i < t.numberOfTethers; ++i)
	{
		xWire = 0.0;
		yWire = ymin + dy * i;

		xcur = x - xWire;
		ycur = y - yWire;
		
		rcur = mag(xcur, ycur);

		Efield_WirePlasma(t, Exi, Eyi, xcur, ycur, r0);
		Ex += Exi;
		Ey += Eyi;
	}
}



struct electron
{
	double density;
	double speed;
	double temperature;
	double mass;
	double charge;
	double debyeLength;
};

void init_electron(electron& e, double d, double s, double t, double m, double c)
{
	e.density = d;
	e.speed = s;
	e.temperature = t;
	e.mass = m;
	e.charge = c;
	e.debyeLength = sqrt(EPSILON0 * BOLTZMANN_CONSTANT * t / (d * sqr(c)));

	if(VERBOSE_INIT){
		cout << "Initializing electron paramerters:\n";
		cout << "  density     = " << d << " 1/m^3\n";
		cout << "  speed       = " << s << " m/s\n";
		cout << "  temperature = " << t << " K\n";
		cout << "              = " << Kelvin_to_eV(t) << " eV\n";
		cout << "  mass        = " << m << " kg\n";
		cout << "              = " << m/MASS_ELECTRON << " m_e\n";
		cout << "  charge      = " << c << " C\n";
		cout << "              = " << c/CHARGE_PROTON << " e\n\n";
	}
	if(HELP_INIT){
		cout << "HELP: init_electron\n";
		cout << "  density     -> The electron volume density. Used in computing the Debye length.\n";
		cout << "                   If tracking electrons, the density is used in computing the current flux.\n";
		cout << "  speed       -> The electron flow speed. Assumed to be flowing in the +x direction.\n";
		cout << "                   If tracking electrons, the speed is used in computing the current flux.\n";
		cout << "  temperature -> The electron temperature. Used in computing the Debye length.\n";
		cout << "                   If tracking electrons, the temperature also gives a Maxwellian thermal\n";
		cout << "                   component to the velocity with the speed as the mean.\n";
		cout << "  mass        -> The electron mass. If tracking particles, used in computing the momentum.\n";
		cout << "  charge      -> The electron charge. Used in computing the Debye length.\n";
		cout << "                   If tracking electrons, the charge is used in computing the electric field.\n";
		cout << "  debyeLength -> The Debye length = sqrt(EPSILON0 * BOLTZMANN_CONSTANT * temperature / (density * sqr(charge))).\n\n";
	}
}


struct ion
{
	double density;
	double speed;
	double temperature;
	double mass;
	double charge;
};

void init_ion(ion& i, double d, double s, double t, double m, double c)
{
	i.density = d;
	i.speed = s;
	i.temperature = t;
	i.mass = m;
	i.charge = c;

	if(VERBOSE_INIT){
		cout << "Initializing ion paramerters:\n";
		cout << "  density     = " << d << " 1/m^3\n";
		cout << "  speed       = " << s << " m/s\n";
		cout << "  temperature = " << t << " K\n";
		cout << "              = " << Kelvin_to_eV(t) << " eV\n";
		cout << "  mass        = " << m << " kg\n";
		cout << "              = " << m/MASS_PROTON << " m_p\n";
		cout << "  charge      = " << c << " C\n";
		cout << "              = " << c/CHARGE_PROTON << " e\n\n";
	}
	if(HELP_INIT){
		cout << "HELP: init_ion\n";
		cout << "  density     -> The ion volume density. If tracking ions, the density is used in computing the current flux.\n";
		cout << "  speed       -> The ion flow speed. Assumed to be flowing in the +x direction.\n";
		cout << "                   If tracking ions, the speed is used in computing the current flux.\n";
		cout << "  temperature -> The ion temperature. If tracking ions, the temperature also gives a Maxwellian thermal\n";
		cout << "                   component to the velocity with the speed as the mean.\n";
		cout << "  mass        -> The ion mass. If tracking particles, used in computing the momentum.\n";
		cout << "  charge      -> The ion charge. If tracking electrons, the charge is used in computing the electric field.\n\n";
	}
}

struct codeOptions
{
	double speedFractionEnergyEquivDistance; // used for initial position away from tether, should be less than 1
	double totalEnergyErrorFraction; // if exceeded, will lower trackEnergyErrorPercent, and rerun all tracks
	double trackEnergyErrorFraction; // if exceeded, will lower epsilonErrorPercent
	double totalForceErrorFraction;  // will iterate tracks until this tolerance is met (use Richardson Extrapolation)
	double dpTrackBoundaryFraction;  // used to find boundary track for integrator, dpTrackBoundaryFraction= dp / p0
	double epsilonErrorFractionRK45;     // is epsilon in the RK45 method, percent error
	double epsilonErrorFractionRomberg;  // is epsilon in the Romber Integration method, percent error
	int domainMacroDivisions;       // number of divisions of the domain (better to be an odd #)
	int minRombergDivisions;        // minimum # of Romberg divisions
	int maxRombergDivisions;        // maximum # of Romberg divisions (if reached, should return warning)
};

void init_codeOptions(codeOptions& co, double sfeed, double toeef, double teep, double tfep, double dtbp, double eeprk,
	double eeprm, int dmd, int mrd1, int mrd2)
{
	co.speedFractionEnergyEquivDistance = sfeed;
	co.totalEnergyErrorFraction = toeef;
	co.trackEnergyErrorFraction = teep;
	co.totalForceErrorFraction = tfep;
	co.dpTrackBoundaryFraction = dtbp;
	co.epsilonErrorFractionRK45 = eeprk;
	co.epsilonErrorFractionRomberg = eeprm;
	co.domainMacroDivisions = dmd;
	co.minRombergDivisions = mrd1;
	co.maxRombergDivisions = mrd2;

	if(VERBOSE_INIT){
		cout << "Initializing codeOptions paramerters:\n";
		cout << "  speedFractionEnergyEquivDistance = " << sfeed << " \n";
		cout << "  totalEnergyErrorFraction         = " << toeef << " \n";
		cout << "  trackEnergyErrorFraction         = " << teep << " \n";
		cout << "  totalForceErrorFraction          = " << tfep << " \n";
		cout << "  dpTrackBoundaryFraction          = " << dtbp << " \n";
		cout << "  epsilonErrorFractionRK45         = " << eeprk << " \n";
		cout << "  epsilonErrorFractionRomberg      = " << eeprm << " \n";
		cout << "  domainMacroDivisions             = " << dmd << " \n";
		cout << "  minRombergDivisions              = " << mrd1 << " \n";
		cout << "  maxRombergDivisions              = " << mrd2 << " \n\n";
	}
	// if(HELP_INIT) {
	// 	cout << "HELP: init_codeOptions\n";
	// 	cout << "  speedFractionEnergyEquivDistance -> ";
	// }
}


struct paramList
{
	tether* tetherP;
	electron* electronP;
	ion* ionP;
	codeOptions* codeOptionsP;
};

// need to give track a particle idenitiy
struct trackVars
{
	double x0;
	double y0;
	double u0;
	double v0;

	double xi;
	double yi;
	double ui;
	double vi;

	double dpx;
	double dpy;
	double Fx;
	double Fy;

	double h; // timestep
};

void init_trackVars(trackVars& tv, double x0, double y0, double u0, double v0, double h)
{
	tv.x0 = x0;
	tv.y0 = y0;
	tv.u0 = u0;
	tv.v0 = v0;
	tv.xi = x0;
	tv.yi = y0;
	tv.ui = u0;
	tv.vi = v0;
	tv.dpx = 0.0;
	tv.dpy = 0.0;
	tv.Fx = 0.0;
	tv.Fy = 0.0;
	tv.h = h;

	if(VERBOSE_INIT > 1){
		cout << "Initializing trackVars paramerters:\n";
		cout << "  x0 = " << x0 << " m\n";
		cout << "  y0 = " << y0 << " m\n";
		cout << "  u0 = " << u0 << " m/s\n";
		cout << "  v0 = " << v0 << " m/s\n";
		cout << "  h  = " << h << " s\n\n";
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// ionP->temperature
// if vtFile is defined, assume it is open
void gen_thermalVel(double& vxT, double& vyT, ion* ionP)
{
	//Generate two random variables 
	double U1 = ((double) rand() / (RAND_MAX));
	double U2 = ((double) rand() / (RAND_MAX));
	//cout << U1 << ' ' << U2 << endl;

	//Calc standard deviation based on particle temp and mass.
	double sigma = sqrt((BOLTZMANN_CONSTANT*ionP->temperature)/ionP->mass);

	//Convert to normal distro through Box-Muller Method and adjusted for standard deviation
	double R = sigma * sqrt(-2.0*log(U1));
	double Theta = 2.0*M_PI*U2;
	vxT = R*cos(Theta);
	vyT = R*sin(Theta);
	//cout << vxT << ' ' << vyT << endl;
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// doesn't need to be initialized. Overridden during runtime
struct RK45VarsPosVel
{
	double kx[6];
	double ky[6];
	double ku[6];
	double kv[6];

	double xi1[2];
	double yi1[2];
	double ui1[2];
	double vi1[2];

	double Rx;
	double Ry;
	double Ru;
	double Rv;

	double deltaX;
	double deltaY;
	double deltaU;
	double deltaV;
};


const double RK45Coeff[8][6] =
	{ // k1         , k2          , k3          , k4           , k5      , k6
		{0.         , 0.          , 0.          , 0.           ,       0., 0.    },
		{1./4.      , 0.          , 0.          , 0.           ,       0., 0.    },
		{3./32.     , 9./32.      , 0.          , 0.           ,       0., 0.    },
		{1932./2197., -7200./2197., 7296./2197. , 0.           ,       0., 0.    },
		{439./216.  , -8.         , 3680./513.  , -845./4104.  ,       0., 0.    },
		{-8./27.    , 2.          , -3544./2565., 1859./4104.  , -11./40., 0.    },
		{25./216.   , 0.          , 1408./2565. , 2197./4104.  , -1./5.  , 0.    },
		{16./135.   , 0.          , 6656./12825., 28561./56430., -9./50. , 2./55.}
	};

// functions ///////////////////////////////////////
// see https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
void RK45UpdatePosVel(trackVars& t,
	                  RK45VarsPosVel& r,
	                  const double c[][6],
	                  paramList& p)
{
	int i, j;
	double Extemp, Eytemp, rcur, xcur, ycur, vcur, ucur, maxR;
	double r0 = 2.0 * p.electronP->debyeLength;

	do {
		for (i = 0; i < 6; ++i) // RK45 steps 1-6
		{
			r.kx[i] = t.ui;
			r.ky[i] = t.vi;
			xcur = t.xi;
			ycur = t.yi;
			for (j = 0; j < i; ++j)
			{
				r.kx[i] += c[i][j] * r.ku[j];
				r.ky[i] += c[i][j] * r.kv[j];
				xcur    += c[i][j] * r.kx[j];
				ycur    += c[i][j] * r.ky[j];
			}
			r.kx[i] *= t.h;
			r.ky[i] *= t.h;
			p.tetherP->EfieldFunction(*p.tetherP, Extemp, Eytemp, xcur, ycur, r0);
			r.ku[i] = t.h * (-p.ionP->charge / p.ionP->mass) * Extemp;
			r.kv[i] = t.h * (-p.ionP->charge / p.ionP->mass) * Eytemp;
		}

		// check error and update timestep, t.h
		for (i = 0; i < 2; ++i)
		{
			r.xi1[i] = t.xi;
			r.yi1[i] = t.yi;
			r.ui1[i] = t.ui;
			r.vi1[i] = t.vi;
			for (j = 0; j < 6; ++j)
			{
				r.xi1[i] += c[i+6][j] * r.kx[j];
				r.yi1[i] += c[i+6][j] * r.ky[j];
				r.ui1[i] += c[i+6][j] * r.ku[j];
				r.vi1[i] += c[i+6][j] * r.kv[j];
			}
		}

		r.Rx = fabs(r.xi1[1] - r.xi1[0]) / t.h;
		r.Ry = fabs(r.yi1[1] - r.yi1[0]) / t.h;
		r.Ru = fabs(r.ui1[1] - r.ui1[0]) / t.h;
		r.Rv = fabs(r.vi1[1] - r.vi1[0]) / t.h;

		maxR = max(r.Rx, r.Ry);
		t.h *= (maxR > 0.0 ? 0.84 * pow(p.codeOptionsP->epsilonErrorFractionRK45 / maxR, 0.25) : 1.2);

		if(DEBUG > 1){
			cout << "Rx = " << r.Rx << endl;
			cout << "Ry = " << r.Rx << endl;
			cout << "Ru = " << r.Rx << endl;
			cout << "Rv = " << r.Rx << endl;
			cout << "timestep = " << t.h << endl;
			cout << "maxR = " << maxR << ' ' << p.codeOptionsP->epsilonErrorFractionRK45 << endl << endl;
		}

	} while(maxR > p.codeOptionsP->epsilonErrorFractionRK45);

	// update position and velocity
	t.xi = r.xi1[0];
	t.yi = r.yi1[0];
	t.ui = r.ui1[0];
	t.vi = r.vi1[0];

	if(DEBUG > 1) {
		cout << "xi = " << t.xi << endl;
		cout << "yi = " << t.yi << endl;
		cout << "ui = " << t.ui << endl;
		cout << "vi = " << t.vi << endl << endl;
	}
}

// compute force from finsihed particle track
void forceOnWire(trackVars& t, paramList& p)
{
	t.Fx = p.ionP->mass * (t.ui - t.u0) * (p.ionP->density) * (t.u0); // N/m
	t.Fy = p.ionP->mass * (t.vi - t.v0) * (p.ionP->density) * (t.v0); // N/m
}



// complete one particle track
void particleTrack(trackVars& t, paramList& p)
{
	int i = 0;
	RK45VarsPosVel RK45Vars;
	double E0 = 0.0, Ef = 0.0;
	E0 = 0.5 * p.ionP->mass * mag2(t.u0, t.v0);
	do{
		i = 0;
		do{
			i++;
			RK45UpdatePosVel(t, RK45Vars, RK45Coeff, p);
			// chance to print particle position here
			if(PRINT_TRACK == 1 && !DEBUG && (i % 1) == 0) {
				cout << i << ' ' << t.xi / (p.electronP->debyeLength) << ' ' << t.yi / (p.electronP->debyeLength) << endl;
			}
			if(PRINT_TRACK == 2){
				if((t.yi > 0) && ((t.xi / (p.electronP->debyeLength) > 74.9) &&
					(t.xi / (p.electronP->debyeLength) < 75.1))) {
					cout << t.x0 / (p.electronP->debyeLength) << ' ' << t.y0 / (p.electronP->debyeLength) << ' ';
					cout << t.xi / (p.electronP->debyeLength) << ' ' << t.yi / (p.electronP->debyeLength) << ' ';
					cout << t.u0 << ' ' << t.v0 << ' ';
					cout << t.ui << ' ' << t.vi << endl;
				}
			}

		}while(mag(t.xi, t.yi) < mag(t.x0, t.y0)); // until reach the equipotential point

		Ef = 0.5 * p.ionP->mass * mag2(t.ui, t.vi);

		//////////////////////////////////////////////////////////////////////////////////////////////
		//******************* NEED TO FIX -> include potential energy somehow *********//////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		// check the energy error
		// if(fabs((Ef-E0) / Ef) > p.codeOptionsP->trackEnergyErrorFraction) {
		// 	if(DEBUG){
		// 		cout << "Ef = " << Ef << endl;
		// 		cout << "E0 = " << E0 << endl;
		// 		cout << fabs((Ef-E0) / Ef) << ' ' << p.codeOptionsP->trackEnergyErrorFraction << endl << endl;
		// 	}
		// 	p.codeOptionsP->trackEnergyErrorFraction *= 0.5;
		// } else {
		// 	p.codeOptionsP->trackEnergyErrorFraction *= 1.05;
		// }
	} while(0/*fabs((Ef-E0) / Ef) > p.codeOptionsP->trackEnergyErrorFraction*/); // This doesn't work right now
}


void forceFromParticleTrack(trackVars& t, paramList& p, double x0, double y0)
{
	// generate thermal component of velocity in x and y directions
	double vxT = 0.0, vyT = 0.0;
	gen_thermalVel(vxT, vyT, p.ionP);

	//ofstream vtFile("vtFile.txt", std::ios_base::app);
	cout << vxT << ' ' << vyT << endl;
	//vtFile.close();


	init_trackVars(t,
	/* x0 */       x0,
	/* y0 */       y0,
	/* u0 */       p.ionP->speed + vxT,
	/* v0 */       vyT,
	/* h  */       0.02 * p.electronP->debyeLength);

	particleTrack(t, p);
	forceOnWire(t, p);
}


void rombergIntegrateParticles(paramList& p, double leftBound, double rightBound, double& FxFinal, double& FyFinal, int &n)
{
	double pow4m, hn, err = 10.0 * p.codeOptionsP->epsilonErrorFractionRomberg;
	int m, k;
	vector<double> Rx, Ry;
	trackVars track_i;
	double FxSum = 0.0, FySum = 0.0;
	double boundaryLength = max(fabs(leftBound), fabs(rightBound)); // m

	n = 1;

	// base case, R(0,0), trapezoid rule
	hn = (rightBound - leftBound) / 2.0;

	forceFromParticleTrack(track_i, p, -boundaryLength, leftBound);
	Rx.push_back(hn * track_i.Fx);
	Ry.push_back(hn * track_i.Fy);
	forceFromParticleTrack(track_i, p, -boundaryLength, rightBound);
	Rx[0] += hn * track_i.Fx;
	Ry[0] += hn * track_i.Fy;

	//cout << Rx[0] << endl;

	while(((n < p.codeOptionsP->minRombergDivisions)
			|| (err > p.codeOptionsP->epsilonErrorFractionRomberg))
			&& n < p.codeOptionsP->maxRombergDivisions)
	{
		Rx.push_back(0.0); // extend R array
		Ry.push_back(0.0);

		FxSum = 0.0;
		FySum = 0.0;

		// fill in func eval gaps for next level
		for (k = 1; k <= (1 << (n-1)); k++)
		{
			forceFromParticleTrack(track_i, p,
				                    -boundaryLength,
				                    leftBound + (2.0*k - 1.0)*hn);
			FxSum += track_i.Fx;
			FySum += track_i.Fy;
			// cout << (leftBound + (2.0*k - 1.0)*hn ) / p.electronP->debyeLength << ' ';
			// cout << track_i.Fx << ' ' << FxSum * hn * 1E9 << endl;
		}

		Rx[triIdx(n,0)] = Rx[triIdx(n-1,0)] / 2.0 + hn * FxSum;
		Ry[triIdx(n,0)] = Ry[triIdx(n-1,0)] / 2.0 + hn * FySum;

		pow4m = 1.0;

		// compute n-th row of m's
		for (m = 1; m <= n; m++)
		{
			Rx.push_back(0.0); // extend R array
			Ry.push_back(0.0);

			pow4m *= 4.0;

			Rx[triIdx(n,m)] = Rx[triIdx(n,m-1)]
			                + (Rx[triIdx(n,m-1)] - Rx[triIdx(n-1,m-1)]) / (pow4m - 1.0);
			Ry[triIdx(n,m)] = Ry[triIdx(n,m-1)]
			                + (Ry[triIdx(n,m-1)] - Ry[triIdx(n-1,m-1)]) / (pow4m - 1.0);
		}

		// compute relative error
		err = fabs((Rx[triIdx(n,n)] - Rx[triIdx(n-1,n-1)]) / Rx[triIdx(n,n)]);
		// cout << n << ' ' << err << ' ' << p.codeOptionsP->epsilonErrorFractionRomberg << ' ';
		// cout << Rx[triIdx(n,n)] << ' ' << Rx[triIdx(n-1,n-1)] << endl; 

		n++;
		hn /= 2.0;
	}
	n--;
	FxFinal = Rx[triIdx(n,n)];
	FyFinal = Ry[triIdx(n,n)];
}



void forceOnTetherRomberg(paramList& p)
{
	double lambdaD = p.electronP->debyeLength;
	double boundaryLength = p.tetherP->tetherSeparation * p.tetherP->numberOfTethers + 100.0 * lambdaD; //100; // Debye lengths
	// outputs
	double FxFinal = 0.0, FyFinal = 0.0, Fxi, Fyi;
	double dyDiv = 2.0 * boundaryLength / double(p.codeOptionsP->domainMacroDivisions);
	int nTot = 0, ni, i;


	for (i = 0; i < p.codeOptionsP->domainMacroDivisions; i++)
	{

		rombergIntegrateParticles(p, -boundaryLength + dyDiv*i,
			                          -boundaryLength + dyDiv*(i+1.0), Fxi, Fyi, ni);

		// cout << i << ' ' << (-boundaryLength + dyDiv*i) / p.electronP->debyeLength << ' ';
		// cout << (-boundaryLength + dyDiv*(i+1.0)) / p.electronP->debyeLength << ' ';
		// cout << Fxi << endl;

		FxFinal += Fxi;
		FyFinal += Fyi;
		nTot += pow(2.0, ni);
	}


	if(PRINT_FORCE == 1) {
		cout << nTot << ' ';
		if(p.tetherP->numberOfTethers > 1)
			cout << p.tetherP->tetherSeparation / lambdaD << ' '; // debye lengths
		cout << mps_to_eV(p.ionP->speed, p.ionP->mass) << ' '; // eV, ion energy
		cout << Kelvin_to_eV(p.electronP->temperature) << ' ';  // eV, electron thermal energy
		cout << p.tetherP->potential << ' '; // eV, tether energy
		cout << -1000 * FxFinal / (p.ionP->density * p.ionP->mass * sqr(p.ionP->speed)) << ' '; // mm, stopping distance
		cout << 1000 * lambdaD << ' '; // mm, Debye length
		cout << density_to_eVpm2(p.electronP->density) << ' '; // 1/m^3
		cout << -FxFinal * 1.E9 << endl; // nN/m, net force in x direction
	}

}



void forceOnTether(paramList& p)
{
	int N = 20000;//200000;
	double lambdaD = p.electronP->debyeLength;
	double boundaryLength = p.tetherP->tetherSeparation * p.tetherP->numberOfTethers / lambdaD + 500.0; //100; // Debye lengths
	double px = 0.0, Fx = 0.0;
	double py = 0.0, Fy = 0.0;
	double dy = 2.0 * lambdaD * boundaryLength / double(N-1.0);
	double vxT, vyT;


	for (int i = 0; i < N; ++i)
	{
		// generate thermal component of velocity in x and y directions
		gen_thermalVel(vxT, vyT, p.ionP);

		// initialize track
		trackVars track_i;
		init_trackVars(track_i,
		/* x0 */       -boundaryLength * lambdaD,
		/* y0 */       -boundaryLength * lambdaD + i * dy,
		/* u0 */       p.ionP->speed + vxT,
		/* v0 */       vyT,
		/* h  */       0.02 * lambdaD);


		particleTrack(track_i, p);

		// eval momentum exchange and force (and impact angle)

		track_i.dpx = p.ionP->mass * (track_i.ui - track_i.u0); // kg * m/s
		track_i.dpy = p.ionP->mass * (track_i.vi - track_i.v0); // kg * m/s

		px += track_i.dpx;
		py += track_i.dpy;

		track_i.Fx = p.ionP->mass * (track_i.ui - track_i.u0) * dy * (p.ionP->density) * (track_i.u0); // N/m
		track_i.Fy = p.ionP->mass * (track_i.vi - track_i.v0) * dy * (p.ionP->density) * (track_i.v0); // N/m

		Fx += track_i.Fx;
		Fy += track_i.Fy;

		if(PRINT_FORCE == 2) { // Force as a function of y position
			cout << (-boundaryLength * lambdaD + i * dy) / lambdaD << ' '; // lambdaD
			cout << -track_i.Fx * 1.E9 << ' ';  // nN/m
			cout << -track_i.Fy * 1.E9 << endl; // nN/m
		}

	}


	if(PRINT_FORCE == 1) {
		cout << N << ' ';
		if(p.tetherP->numberOfTethers > 1)
			cout << p.tetherP->tetherSeparation / lambdaD << ' '; // debye lengths
		cout << mps_to_eV(p.ionP->speed, p.ionP->mass) << ' '; // eV, ion energy
		cout << Kelvin_to_eV(p.electronP->temperature) << ' ';  // eV, electron thermal energy
		cout << p.tetherP->potential << ' '; // eV, tether energy
		cout << -1000 * Fx / (p.ionP->density * p.ionP->mass * sqr(p.ionP->speed)) << ' '; // mm, stopping distance
		cout << 1000 * lambdaD << ' '; // mm, Debye length
		cout << p.electronP->density << ' '; // 1/m^3
		cout << -Fx * 1.E9 << endl; // nN/m, net force in x direction

		// cout << "dy = " << dy*1000 << " mm" << endl;
		// cout << "Net Force: Fx = " << Fx * 1.E9 << " nN/m,  Fy = " << Fy * 1.E9 << " nN/m\n";
		// cout << "Net Force: Fx = " << Fx * 1.E3 << " mN/m,  Fy = " << Fy * 1.E3 << " mN/m\n";
		// cout << "Ira comp  = " << 1000 * Fx / (p.ionP->density * p.ionP->mass * sqr(p.ionP->speed)) << " mm\n";
		// // cout << "Net Momentum Exchange: px = " << px * dy << " kg*m/s\n";
	}

	// first find integration bounds, ymax

	// loop over the domain in sections, preferably an odd # (can be 1)

		// for each section, apply Romberg integration 

}



int main(int argc, char const *argv[])
{
	tether eSailTether;
	electron labElectrons;
	ion labArgonIons;
	codeOptions options;
	paramList paramerters;

	paramerters.tetherP      = &eSailTether;
	paramerters.electronP    = &labElectrons;
	paramerters.ionP         = &labArgonIons;
	paramerters.codeOptionsP = &options;


	int N = 50, M = 10, P = 1, i, j, k;
	
	double Te_min = 1e-8;//0.1;  // eV
	double Te_max = 100.0;//100.0; // eV

	double sep_min = 0.15; // debye lengths
	double sep_max = 20.0; // debye lengths

	double ion_min = 0.1;//0.01;   // x potential
	double ion_max = 100.0;//10.0; // x potential

	double potential_min = 1.0;    // V
	double potential_max = 1000.0; // V
	double cur_potential;

	double edensity_min = 1.0e11;
	double edensity_max = 1.0e13;
	double cur_edensity;

	for (k = 0; k < P; ++k)
	{
		for (j = 0; j < M; ++j)
		{

			for (i = 0; i < N; ++i)
			{
				cur_edensity  = edensity_min * pow(edensity_max/edensity_min, j/double(M-1.0));
				//cur_potential = potential_min * pow(potential_max/potential_min, j/double(M-1.0));

				init_electron(labElectrons,
				/* density     */ cur_edensity,//130.0e-2/CHARGE_PROTON/eV_to_mps(50.0, 131.293*MASS_PROTON),//1.00e12, // m^-3
				/* speed       */ 0.0, // m/s
				/* temperature */ eV_to_Kelvin(12.0),//eV_to_Kelvin(Te_min * density_to_eVpm2(cur_edensity) * pow(Te_max/Te_min, i/double(N-1.0))), // eV -> K
				/* mass        */ MASS_ELECTRON, // kg
				/* charge      */ -CHARGE_PROTON); // C

				init_ion(labArgonIons,
			    /* density     */ 1.0e12, //130.0e-2/CHARGE_PROTON/eV_to_mps(50.0, 131.293*MASS_PROTON), //1.8e12,//1.00e12, // m^-3
				/* speed       */ eV_to_mps(1000.0, 131.293*MASS_PROTON),//eV_to_mps(ion_min * cur_potential * pow(ion_max/ion_min, i/double(N-1.0)), 131.293*MASS_PROTON), // eV -> m/s  105
				/* temperature */ eV_to_Kelvin(10.0),//0.0, // K
				/* mass        */ 131.293*MASS_PROTON,//39.948*MASS_PROTON, // kg
				/* charge      */ CHARGE_PROTON); // C

				init_tether(eSailTether,
				/* potential        */ 10000.0,//100.0, // V 
				/* radius           */ 1.0e-3, //7.874e-4 / 2.0, //1.0e-3, // m
				/* numberOfTethers  */ 1,
				/* tetherSeparation */ 3.0 * labElectrons.debyeLength,//(sep_min * pow(sep_max/sep_min, j/double(M-1.0))) * labElectrons.debyeLength,//1e-1, // m
				/* EfieldFunction   */ Efield_WirePlasma);

				
					
				init_codeOptions(options,
				/* speedFractionEnergyEquivDistance */ 1.0e-6, // fraction
				/* totalEnergyErrorFraction         */ 1.0e-3, // fraction
				/* trackEnergyErrorFraction         */ 5.0e-2, // fraction
				/* totalForceErrorFraction          */ 1.0e-5, // fraction
				/* dpTrackBoundaryFraction          */ 1.0e-6, // fraction
				/* epsilonErrorFractionRK45         */ 1.0e-4, // fraction
				/* epsilonErrorFractionRomberg      */ 1.0e-2, // fraction
				/* domainMacroDivisions             */ 1,
				/* minRombergDivisions              */ 8,
				/* maxRombergDivisions              */ 20);

				forceOnTetherRomberg(paramerters);
				//forceOnTether(paramerters);
				
				
		 	}
		 	// if(PRINT_FORCE){
				// 	cout << endl << endl;
				// }
		}
		// if(PRINT_FORCE){
		// 			cout << endl << endl;
		// 		}
	}

	return 0;
}

	// double speedFractionEnergyEquivDistance; // used for initial position away from tether, should be less than 1
	// double totalEnergyErrorFraction; // if exceeded, will lower trackEnergyErrorPercent, and rerun all tracks
	// double trackEnergyErrorFraction; // if exceeded, will lower epsilonErrorPercent
	// double totalForceErrorFraction;  // will iterate tracks until this tolerance is met (use Richardson Extrapolation)
	// double dpTrackBoundaryFraction;  // used to find boundary track for integrator, dpTrackBoundaryFraction= dp / p0
	// double epsilonErrorFractionRK45;     // is epsilon in the RK45 method, percent error
	// double epsilonErrorFractionRomberg;  // is epsilon in the Romber Integration method, percent error
	// int domainMacroDivisions;       // number of divisions of the domain (better to be an odd #)
	// int minRombergDivisions;        // minimum # of Romberg divisions
	// int maxRombergDivisions;   