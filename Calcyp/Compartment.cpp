#include "Compartment.h"
#include <iomanip>

Compartment::Compartment(int Index, int area, double wieldingpoint, double fieldcapacity, int thick , double CCa0, double CSO4)
{
	this->fDepth = (Index +0.5)*thick;
	this->nThetaWeildingPnt = wieldingpoint;
	this->nFieldCapacity = fieldcapacity;
	this->nWhc = fieldcapacity;
	this->nMoist = nThetaWeildingPnt;
	this->C_Ca = CCa0;
	this->C_SO4 = CSO4;
	this->C_CaSO4 = 0;
	this->nthick = thick;
	this->nArea = area;
	this->nLastMoist = 0;
	this->nWetCount = 0;
	this->fTotLeachate = 0;
	this->fAETLoss = 0;
	this->nFloodedCount = 0;
}

Compartment::~Compartment()
{
}




// solubility is a function that calculates the mass - caco3 [g] of CaCO3 and concentration of Ca in solution - cca [mol/L] or [M] in a compartment given the temprature - temp [deg]
// partial pressure of CO2 - pco2 [atm], intial concentration of Ca in solution - cca [mol/L] or [M], current soil moisture - moist [cm], previous soil moisture content - moisti [cm]
// initial mass of CaCO3 - caco3 [g], compartment thickness - thick [cm], and compartment area - area [cm2].
// The calculation is following Marion et al. (1985) and Hirmas et al. (2010)
double Compartment::solubility(double temp)
{
	double I, A, Ksp, ionActivity, MCa, MSo4, MCaSO4,a_Ca,a_So4,a_CaSo4, // a_* is activity = [Molar * ion activity]
		GypOmega, alphaGypsum, a, b, c;
	pair<double, double> Quadsolutions;

	//convert to Molar
	double MoistInLitre = this->nMoist * 0.001F;
	MCa = C_Ca / MoistInLitre;//Molar
	MSo4 = C_SO4 / MoistInLitre;//Molar
	MCaSO4 = C_CaSO4 / MoistInLitre;//Molar

	// Ionic strength
	I = 0.5*(MCa * 4 + MSo4 * 4); // eq. 4. as cca is in [M] or [mol/L], I is also in [M] or [mol/L] 
	A = 0.4918 + (6.6098 * pow(10, -4) * temp) + (5.0231 * pow(10, -6) * pow(temp, 2));
	ionActivity = pow(10, -A * sqrt(I) / (1 + sqrt(I)) - 0.3*I);
	Ksp = pow(10, -4.58);

	//Eq.3 converting to activity
	a_Ca = MCa * ionActivity;
	a_So4 = MSo4 * ionActivity;
	a_CaSo4 = MCaSO4 * ionActivity;
	
	// omega is the difference between current product and equil product
	//Ω > 1 - supersaturated solution, and Ω < 1 - undersaturated solution
	GypOmega =  a_Ca * a_So4/Ksp;
// the concentration we need to add or substract to gypsum
	a = 1;
	b = -(a_Ca + a_So4);
	c = (a_Ca * a_So4) - Ksp;
	Quadsolutions = quadricEquation(a, b, c);
	GypOmega = (round(GypOmega * 10) / 10.0);	
	//cout << GypOmega << endl;
	// percipitation
	if (GypOmega > 1)
	{
		// Check which one of the two solution produce positive value
		alphaGypsum = nonNegetiveX(Quadsolutions, a_Ca, a_So4);		
		alphaGypsum *= -1;
	}
	//dissolution
	else if(GypOmega < 1 && GypOmega > 0)
	{
		// Check which one of the two solution produce positive value
		alphaGypsum = nonNegetiveX(Quadsolutions, a_CaSo4);		
	}
	else {
		alphaGypsum = 0;
	}

	//Rcpp::Rcout << "ag  "<<Quadsolutions.first<< "  " << Quadsolutions.second<< endl;
	a_Ca += alphaGypsum;
	a_So4 += alphaGypsum;
	a_CaSo4 -= alphaGypsum;
	
	//convert to mol
	C_Ca = a_Ca / ionActivity * MoistInLitre;
	C_SO4 = a_So4 / ionActivity * MoistInLitre;
	C_CaSO4 = a_CaSo4/ ionActivity * MoistInLitre;
	setAllToZero();

	//printf_s("Percipitation: omegaga is %f    gypsum is: %f \n", GypOmega, MCaSO4 * MoistInLitre);

	return C_CaSO4;
}

pair<double, double> Compartment::quadricEquation(double a, double b , double c)
{
	pair<double, double> solutions;
	solutions.first = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	solutions.second = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	return solutions;
}
     

void Compartment::setAllToZero()
{
	{
		if (C_Ca < 0) C_Ca = 0;
		if (C_SO4 < 0) C_SO4 = 0;
		if (C_CaSO4 < 0) C_CaSO4 = 0;
	}
}

double Compartment::GetIOnsSum()
{
	return C_CaSO4 + C_SO4 + C_Ca;
}

double Compartment::nonNegetiveX(pair<double, double>& sol, double Ca, double SO4)
{
	if ((sol.first <= Ca) && (sol.first <= SO4)) {
		return(sol.first);
	}
	else {
		return(sol.second);
	}
}

double Compartment::nonNegetiveX(pair<double, double>& sol, double gyp)
{
	if (sol.first <= gyp) {
		return(sol.first);
	}
	else if (sol.second < 0) {
		return(gyp);
	}
	else {
		return(sol.second);
	}
}
