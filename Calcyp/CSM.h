#pragma once
#ifndef _CSM_H
#define _CSM_H

#include <vector>
#include <iostream>
#include <stack>
#include <stdio.h>
#include "Compartment.h"
#include "Month.h"
#include <algorithm>
//#include <Rcpp.h>


//using namespace Rcpp;

using namespace std;




class CSM
{
public:
	CSM();
	vector<double> Calculate(vector<double>, vector<double>,int years, int Depth, int nthick, double WieltingPoint, double InitialCa,
						double initialSO4, double nBulkDensity, int FieldArea, double FieldCapacity, double DustCa, double DustSO4, double AETFactor);
	std::vector<Compartment>* GetCompartments();
	double meqSoil2molar(double, double, double moisture);
	double mol2meqSoil(double, double);
	double moistcm2Litre(double);
	~CSM();
	double GetPrecision(double);
	//Rcpp::List GetResults();
	int nNumOfDays;
	int nDepth;
	int thick;
	int nNumOfCompatments;
	int nArea;
	double nTotalCaDust;
	double nTotalRain;
	double nTotalCaLeachate;
	double nTotalSO4Leachate;
	double nTotalLeachate;
	double nTotalMoist;
	double nInitMoistTotal;
	double nTotalAet;
	double nTotalWP;
	double nTotalWhc;
	void InitMonths();
	double BulkDensity;
	double accumolateDustDays;
	
	//Rcpp::List results;
	void InitCompartments();


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	vector<double> RainArr;
	double AET;

	
	double nDailyDustCa;
	double nDailyDustSO4;
	double wieltingPoint;
	double nFieldCapacity;
	double CCa;
	double CSO4;
	double nCaco3less;
	
	
	double nTemp;
	double nLeachate;

	double TempArr[12] = { 20,21,24.7,29.7,31.4,37.6,39.4,39.4,36.5,32,26.5,27.1};
};

//RCPP_EXPOSED_CLASS(CSM);


#endif