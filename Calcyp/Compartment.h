#pragma once
#include "math.h"
#include <stdio.h>
#include <map>
#include <Rcpp.h>


using namespace std;
class Compartment
{
public:
	Compartment(int , int, double, double, int , double, double);
	
	~Compartment();
	
	double solubility(double temp);
	pair<double, double> quadricEquation(double, double, double);
	void setAllToZero();
	double GetIOnsSum();
	double nonNegetiveX(pair<double, double>&, double, double);
	double nonNegetiveX(pair<double, double>&, double);

	double fDepth;
	int nthick;
	int nArea;
	double nThetaWeildingPnt;
	double nFieldCapacity;
	double nWhc;					//water holding capacity
	double nInitTotMoist;
	double nMoist;				// moisture level
	double C_Ca;					//soluble Ca ion [mmol]
	double C_SO4;				//solunle so4 anion [mmol]
	double C_CaSO4;				//solid gypsum in mmol
	double nLastMoist;
	double nCO2;
	double nSolubleCa;
	double nTotMoist;			//summing the total moisture in the soil profile. if there are 10 compartments of 10cm each, total moisture is 8.05cm
	double CcaLess;
	double CSO4Less;
	int nWetCount;
	double fTotLeachate;
	double fAETLoss;
	int nFloodedCount;

};

