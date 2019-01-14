#pragma once
class Compartment
{
public:
	Compartment(int , int, float, float, float , float , float);
	
	~Compartment();
	
	float solubility(float temp);
	void setAllToZero();
	float GetIOnsSum();

	int nIndex;
	int nthick;
	int nArea;
	float nThetaWeildingPnt;
	float nFieldCapacity;
	float nWhc;					//water holding capacity
	float	nInitMoist;
	float nMoist;				// moisture level
	float	nCaCO3;
	float nCCa;					//soluble Ca ion [M]
	float C_SO4;				//solunle so4 anion [M]
	float C_CaSO4;				//solid gypsum in Molar
	float nICO2;
	float nCO2;
	float nSolubleCa;
	float nTotMoist;			//summing the total moisture in the soil profile. if there are 10 compartments of 10cm each, total moisture is 8.05cm
	float	nTotWhc;			// summing the whc in the soil profile. if there are 10 compartments of 10cm each, total whc is 8.30cm
	float	nTotThetaWP;		 //summing the wp in the soil profile. if there are 10 compartments of 10cm each, total wp is 3.9cm
	float	nInitTotalMoist;
	float CcaLess;
	float CSO4Less;

};

