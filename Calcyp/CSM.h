#pragma once


#include <stdio.h>
#include <vector>
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
	Rcpp::List Calculate(Rcpp::NumericVector , float years, float Depth, float nthick, float nwieltingPoint, float InitialCa, float initialSO4,
		float nBulkDensity, float FieldArea, float FieldCapacity, float DustCa, float DustSO4);
	std::vector<Compartment>* GetCompartments();
	float meq2mmol(float, float);
	float mmol2meq(float, float);
	~CSM();
	float GetPrecision(float);
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
	

	void InitCompartments();
	MONTH GetMonth(int nDay);
	float JULIAN(int day);


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	NumericVector RainArr;
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

	double TempArr[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
};

RCPP_EXPOSED_CLASS(CSM);


