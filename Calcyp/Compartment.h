#pragma once
#include "math.h"
#include <stdio.h>
#include <map>
#include <Rcpp.h>


using namespace std;
class Compartment
{
public:
	Compartment(int , int, float, float, float , float , float);
	
	~Compartment();
	
	float solubility(float temp);
	pair<float, float> quadricEquation(float, float, float);
	void setAllToZero();
	float GetIOnsSum();

	float fDepth;
	int nthick;
	int nArea;
	float nThetaWeildingPnt;
	float nFieldCapacity;
	float nWhc;					//water holding capacity
	float nInitTotMoist;
	float nMoist;				// moisture level
	float nCCa;					//soluble Ca ion [mmol]
	float C_SO4;				//solunle so4 anion [mmol]
	float C_CaSO4;				//solid gypsum in mmol
	float nLastMoist;
	float nCO2;
	float nSolubleCa;
	float nTotMoist;			//summing the total moisture in the soil profile. if there are 10 compartments of 10cm each, total moisture is 8.05cm
	float CcaLess;
	float CSO4Less;
	int nWetCount;
	float fTotLeachate;
	float fAETLoss;
	int nFloodedCount;

};

