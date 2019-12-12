#include "Compartment.h"


Compartment::Compartment(int Index, int area, float wieldingpoint, float fieldcapacity, float thick , float CCa0, float CSO4)
{
	this->fDepth = (Index +0.5)*thick;
	this->nThetaWeildingPnt = wieldingpoint;
	this->nFieldCapacity = fieldcapacity;
	this->nWhc = fieldcapacity;
	this->nMoist = nThetaWeildingPnt;
	this->nCCa = CCa0;
	this->C_SO4 = CSO4;
	this->C_CaSO4 = 0;
	this->nthick = thick;
	this->nArea = area;
	this->nLastMoist = 0;
	this->nWetCount = 0;
	this->fTotLeachate = 0;
	this->fAETLoss = 0;
}

Compartment::~Compartment()
{
}




// solubility is a function that calculates the mass - caco3 [g] of CaCO3 and concentration of Ca in solution - cca [mol/L] or [M] in a compartment given the temprature - temp [deg]
// partial pressure of CO2 - pco2 [atm], intial concentration of Ca in solution - cca [mol/L] or [M], current soil moisture - moist [cm], previous soil moisture content - moisti [cm]
// initial mass of CaCO3 - caco3 [g], compartment thickness - thick [cm], and compartment area - area [cm2].
// The calculation is following Marion et al. (1985) and Hirmas et al. (2010)
float Compartment::solubility(float temp)
{
	float I, A, k6, ionActivity, EquilConcentrationConstant, MCa, MSo4, MCaSO4,
		CurrentConcentrationProduct, GypOmega, alphaGypsum, a, b, c, limitation;
	pair<float, float> Quadsolutions;

	//convert to Molar
	float MoistInLitre = this->nMoist * 0.01;
	MCa = nCCa / MoistInLitre / 1000;
	MSo4 = C_SO4 / MoistInLitre / 1000;
	MCaSO4 = C_CaSO4 / MoistInLitre / 1000;
	k6 = pow(10, -2.23 - (0.0019*temp));

	// Ionic strength
	I = 0.5*(MCa * 4 + MSo4 * 4); // eq. 14. as cca is in [M] or [mol/L], I is also in [M] or [mol/L] 
	A = 0.4918 + (6.6098 * pow(10, -4) * temp) + (5.0231 * pow(10, -6) * pow(temp, 2));
	ionActivity = pow(10, -A * 4 * (sqrt(I) / (1 + sqrt(I)) - (0.3*I)));

	// check if percipitating
	EquilConcentrationConstant = k6 / (pow(ionActivity, 2));
	CurrentConcentrationProduct = MCa * MSo4;

	// omega is the difference between current product and equil product
	GypOmega = CurrentConcentrationProduct - EquilConcentrationConstant;

	// the concentration we need to add or substract to gypsum
	alphaGypsum = 0;
	a = 1;
	b = (MCa + MSo4);
	c = GypOmega;
	Quadsolutions = quadricEquation(a, b, c);
	alphaGypsum = fmaxf(Quadsolutions.first, Quadsolutions.second);

	// percipitation
	if (GypOmega >= 0)
	{
		limitation = fminf(MCa, MSo4);
		if (limitation < -alphaGypsum)
		{
			alphaGypsum = -limitation;
		}
	}
	//dissolution
	else
	{
		if (MCaSO4 <= alphaGypsum)
		{
			alphaGypsum = MCaSO4;
		}
	}

	//Rcout << "ag  "<<Quadsolutions.first<< "  " << Quadsolutions.second<< endl;
	MCa += alphaGypsum;
	MSo4 += alphaGypsum;
	MCaSO4 -= alphaGypsum;
	
	//convert to mmol
	nCCa = MCa * MoistInLitre*1000;
	C_SO4 = MSo4 * MoistInLitre*1000;
	C_CaSO4 = MCaSO4 * MoistInLitre*1000;
	setAllToZero();

	//printf_s("Percipitation: omegaga is %f    gypsum is: %f \n", GypOmega, MCaSO4 * MoistInLitre);

	return C_CaSO4;
}

pair<float, float> Compartment::quadricEquation(float a, float b , float c)
{
	pair<float, float> solutions;
	solutions.first = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	solutions.second = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	return solutions;
}
     

void Compartment::setAllToZero()
{
	{
		if (nCCa < 0) nCCa = 0;
		if (C_SO4 < 0) C_SO4 = 0;
		if (C_CaSO4 < 0) C_CaSO4 = 0;
	}
}

float Compartment::GetIOnsSum()
{
	return C_CaSO4 + C_SO4 + nCCa;
}
