#pragma once
class Compartment
{
public:
	Compartment(int Index, int area, float wieldingpoint = 0, float fieldcapacity = 0, float thick = 0, float CCa0 = 0, float CO2 = 0);
	
	~Compartment();

	float PCO2(int nDay);
	float solubility(float temp);

	int nIndex;
	int nthick;
	int nArea;
	float nThetaWeildingPnt;
	float nWhc;					//water holding capacity
	float	nInitMoist;
	float nMoist;				// moisture level
	float	nCaCO3;
	float nCCa;
	float nICO2;
	float nCO2;
	float nSolubleCa;
	float nTotMoist;			//summing the total moisture in the soil profile. if there are 10 compartments of 10cm each, total moisture is 8.05cm
	float	nTotWhc;			// summing the whc in the soil profile. if there are 10 compartments of 10cm each, total whc is 8.30cm
	float	nTotThetaWP;		 //summing the wp in the soil profile. if there are 10 compartments of 10cm each, total wp is 3.9cm
	float	nInitTotalMoist;

};

