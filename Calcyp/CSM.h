#pragma once
#ifndef _CSM_H
#define _CSM_H

#include <vector>
#include <stack>
#include <stdio.h>
#include "Compartment.h"
#include "Month.h"
#include <Rcpp.h>


using namespace Rcpp;

using namespace std;


typedef enum { WINTER, SPRING, SUMMER } SEASON;


class CSM
{
public:
	CSM();
	Rcpp::List Calculate(Rcpp::DoubleVector, Rcpp::DoubleVector,float years, float Depth, float nthick, float WieltingPoint, float InitialCa,
						float initialSO4, float nBulkDensity, float FieldArea, float FieldCapacity, float DustCa, float DustSO4, float AETFactor);
	std::vector<Compartment>* GetCompartments();
	float meq2mmol(float, float);
	float mmol2meq(float, float);
	~CSM();
	float GetPrecision(float);
	Rcpp::List GetResults();
	int nNumOfDays;
	float nDepth;
	float thick;
	int nNumOfCompatments;
	float nArea;
	float nTotalCaDust;
	float nTotalRain;
	float nTotalCaLeachate;
	float nTotalSO4Leachate;
	float nTotalLeachate;
	float nTotalMoist;
	float nInitMoistTotal;
	float nTotalAet;
	float nTotalWP;
	float nTotalWhc;
	void InitMonths();
	float BulkDensity;
	float accumolateDustDays;
	
	Rcpp::List results;
	void InitCompartments();
	MONTH GetMonth(int nDay);
	float JULIAN(int day);


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	Rcpp::NumericVector RainArr;
	float AET;

	
	float nDailyDustCa;
	float nDailyDustSO4;
	float wieltingPoint;
	float nFieldCapacity;
	float CCa;
	float CSO4;
	float nCaco3less;
	
	
	float nTemp;
	float nLeachate;

	double TempArr[12] = { 20,21,24.7,29.7,31.4,37.6,39.4,39.4,36.5,32,26.5,27.1};
};

RCPP_EXPOSED_CLASS(CSM);


#endif