#include "Compartment.h"
#include "math.h"

Compartment::Compartment(int Index, int area, float wieldingpoint, float fieldcapacity, float thick , float CCa0, float CO2)
{
	this->nIndex = Index;
	this->nThetaWeildingPnt = wieldingpoint;
	this->nWhc = fieldcapacity - wieldingpoint;
	this->nInitMoist = 0.1*thick;
	this->nMoist = this->nInitMoist;
	this->nCaCO3 = 0;
	this->nCCa = CCa0;
	this->nICO2 = CO2;
	this->nSolubleCa = this->nICO2;
	this->nTotMoist += this->nInitMoist;
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
		double k1, k2, k3, k4, I, g1, g2, A, B, C, ah0, ah1, ah, aco2, ahco3, aco3, aca, ccacalc, mcacalc, mcas0, mca0, mcaeq, mcas1, mca1;

		// Equilibrium constants for the carbonate-bicarbonate system as a function of temprature
		// Notation given by MArion et al. (1985)
		k1 = pow(10, -1.14 - (0.0131*temp)); // eq. 5
		k2 = pow(10, -6.54 + (0.0071*temp)); // eq. 7
		k3 = pow(10, -10.59 + (0.0102*temp)); // eq. 9
		k4 = pow(10, -7.96 - (0.0125*temp)); // eq. 11

											 // Ionic strength
		I = 3 * nCCa; // eq. 14. as cca is in [M] or [mol/L], I is also in [M] or [mol/L] 

					 // Activity coefficients for ionic valences of 1 and 2
		g1 = pow(10, -0.505*(sqrt(I) / (1 + sqrt(I)) - (0.3*I))); // eq. 13
		g2 = pow(10, -0.505 * 4 * (sqrt(I) / (1 + sqrt(I)) - (0.3*I))); // eq. 13

																		// Hydrogen ion activity
		A = 2 * k4 / (g2*k3*k2*k1*nCO2); // Left term in eq. 16
		B = k1*k2*nCO2 / g1; // Middle term in eq. 16
		C = 2 * k1*k2*k3*nCO2 / g2; // Right term in eq. 16

									// Newton's method for solving eq. 16
		ah0 = pow(10, -8); // ph=8 intial guess
		ah1 = pow(10, -4); // 
		while (fabs(-1 * log(ah1) + log(ah0))>0.1)
		{
			ah0 = ah1;
			ah1 = ah0 - ((A*pow(ah0, 4) - (B*ah0) - C) / ((4 * A*pow(ah0, 3)) - B));
		}
		ah = ah1;

		// Equalibrium activities of CO2, HCO3, CO3, and Ca
		aco2 = k1*nCO2; // eq. 4
		ahco3 = k2*aco2 / ah; // eq. 6
		aco3 = k3*ahco3 / ah; // eq. 8
		aca = k4 / aco3; // eq. 10

						 // Ca concentration in [mol/L] or [M]
		ccacalc = aca / g2; // eq. 12

							// Determination of final grams of CaCO3 and Ca in solution within a compartment

							// Ca in solution in compartment volume [g]
		mcacalc = ccacalc*40.1*nMoist*nArea / 1000; // molar mass of Ca is 40.1 [g/mol], moist is the water amount in [cm] of the compartment, area is in [cm^2], and dividing by 1000 is to transfer from Liter to cm^3

												  // Ca in the soil phase in compartment volume [g/cm^3]
		mcas0 = nCaCO3 / 100.1*40.1*nthick*nArea; // molar mass of CaCO3 is 100.1 [g/mol]

												// Initial Ca in solution in compartment volume [g]
												// percolating from upper compartment + same compartment previous time step
		mca0 = nCCa*40.1*nMoist / 1000 * nArea; // similar to calculating mcacalc

											 // Ca needed in solution to reach equilibrium [g]		
		mcaeq = mcacalc - mca0;

		// Difference in Ca between Ca in soil phase to Ca in solution [g]
		mcas1 = mcas0 - mcaeq;

		if (mcas1 >= 0) mca1 = mcacalc; // There is enough Ca in solution to reach equilibrium
		else						   // There is not enough Ca in solution to reach equilibrium	
		{
			mca1 = mca0 + mcas0; // Ca in solution equals Ca from solid + Ca from intial solution
			mcas1 = 0;
		}


		// Conversions
		nCaCO3 = mcas1*100.1 / 40.1 / nthick / nArea; // caco3 [g/cm3]
		nCCa = mca1 / 40.1 / nMoist / nArea * 1000; // Ca in solution [mol/L] or [M]

		return nCaCO3;

	}
}