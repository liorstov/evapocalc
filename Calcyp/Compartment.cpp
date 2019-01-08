#include "Compartment.h"
#include "math.h"

Compartment::Compartment(int Index, int area, float wieldingpoint, float fieldcapacity, float thick , float CCa0, float CSO4)
{
	this->nIndex = Index;
	this->nThetaWeildingPnt = wieldingpoint;
	this->nFieldCapacity = fieldcapacity;
	this->nWhc = fieldcapacity;
	this->nMoist = 0.1*thick;
	this->nCaCO3 = 0;
	this->nCCa = CCa0;
	this->C_SO4 = CSO4;
	this->C_CaSO4 = 0;
	//this->nICO2 = CO2;
	//this->nSolubleCa = this->nICO2;
	this->nTotWhc += this->nWhc*thick;
	this->nTotThetaWP += this->nThetaWeildingPnt*thick;
	this->nInitTotalMoist = this->nTotMoist;
	this->nthick = thick;
	this->nArea = area;
}





Compartment::~Compartment()
{
}

// PCO2 is a function that calculates the partial pressure of CO2 in the soil atmosphere - pco2 [ppm] based on the soil depth [cm].
// We used equations derived from data presented in Table A1 of Marion et al. (2008).

float Compartment::PCO2(int nDayjulian)
{
	float zt, zb, pco2, pco2t, pco2b;

	zt = (nIndex - 1)*nthick;
	zb = nIndex*nthick;

	if ((nDayjulian>68) & (nDayjulian<140)) // spring
	{
		pco2t = (-0.037 / 3 * pow(zt, 3)) + (15.36 / 2 * pow(zt, 2)) + (653.36*zt);
		pco2b = (-0.037 / 3 * pow(zb, 3)) + (15.36 / 2 * pow(zb, 2)) + (653.36*zb);
	}
	else
		if ((nDayjulian>139) & (nDayjulian<228)) // summer
		{
			pco2t = (-0.053 / 3 * pow(zt, 3)) + (16.39 / 2 * pow(zt, 2)) + (698.60*zt);
			pco2b = (-0.053 / 3 * pow(zb, 3)) + (16.39 / 2 * pow(zb, 2)) + (698.60*zb);
		}
		else
			if ((nDayjulian>227) & (nDayjulian<329)) // fall
			{
				pco2t = (0.023 / 3 * pow(zt, 3)) + (1.40 / 2 * pow(zt, 2)) + (664.45*zt);
				pco2b = (0.023 / 3 * pow(zb, 3)) + (1.40 / 2 * pow(zb, 2)) + (664.45*zb);
			}
			else // winter
			{
				pco2t = (-0.042 / 3 * pow(zt, 3)) + (9.33 / 2 * pow(zt, 2)) + (775.62*zt);
				pco2b = (-0.042 / 3 * pow(zb, 3)) + (9.33 / 2 * pow(zb, 2)) + (775.62*zb);
			}

	pco2 = (pco2b - pco2t) / nthick / pow(10, 6); // divide by 10^6 to transfer from [ppm] to [atm]

	return pco2;
}


// solubility is a function that calculates the mass - caco3 [g] of CaCO3 and concentration of Ca in solution - cca [mol/L] or [M] in a compartment given the temprature - temp [deg]
// partial pressure of CO2 - pco2 [atm], intial concentration of Ca in solution - cca [mol/L] or [M], current soil moisture - moist [cm], previous soil moisture content - moisti [cm]
// initial mass of CaCO3 - caco3 [g], compartment thickness - thick [cm], and compartment area - area [cm2].
// The calculation is following Marion et al. (1985) and Hirmas et al. (2010)
float Compartment::solubility(float temp)
{

	{
		float I, A, k6, ionActivity, EquilConcentrationConstant,
			CurrentConcentrationProduct, GypOmega, alphaGypsum, totalConcentrationProduct, a, b, c;

		 k6 = pow(10, -2.23 - (0.0019*temp));

											 // Ionic strength
		I = 0.5*(nCCa*4+C_SO4*4); // eq. 14. as cca is in [M] or [mol/L], I is also in [M] or [mol/L] 
		A = 0.4918 + (6.6098 * pow(10, -4) * temp) + (5.0231 * pow(10, -6) * pow(temp, 2));
		ionActivity = pow(10, -A * 4 * (sqrt(I) / (1 + sqrt(I)) - (0.3*I)));
		
		// check if percipitating
		EquilConcentrationConstant = k6 / (pow(ionActivity, 2));
		CurrentConcentrationProduct = nCCa * C_SO4;

		// omega is the difference between current product and equil product
		GypOmega = CurrentConcentrationProduct - EquilConcentrationConstant;
		totalConcentrationProduct = nCCa * C_SO4;

		// the concentration we need to add or substract to gypsum
		alphaGypsum = 0;
		

		// check if  percipitation occurs
		if (GypOmega >= 0 ) 
		{
			// solving polinum
			a = 1;
			b = -(nCCa + C_SO4);
			c = (CurrentConcentrationProduct - GypOmega);
			alphaGypsum = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
			nCCa -= alphaGypsum;
			C_SO4 -= alphaGypsum;
			C_CaSO4 += alphaGypsum;			
		}
		// gypsum dissolution
		else {
			a = 1;
			b = (nCCa + C_SO4);
			c = (CurrentConcentrationProduct + GypOmega);
			alphaGypsum = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
			nCCa += alphaGypsum;
			C_SO4 += alphaGypsum;
			C_CaSO4 -= alphaGypsum;
		}

		setAllToZero();
		
		return C_CaSO4;
		
	}
	
}

void Compartment::setAllToZero()
{
	{
		if (nCCa < 0) nCCa = 0;
		if (C_SO4 < 0) C_SO4 = 0;
		if (C_CaSO4 < 0) C_CaSO4 = 0;
	}
}
