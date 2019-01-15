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
	CSM(void);
	Rcpp::List Calculate(Rcpp::NumericVector , int , int , int , float , float , float , 
		float , float , float , float , float );
	std::vector<Compartment>* GetCompartments();
	void test(int l, float b = 5);
	float meq2mmol(float, float);
	float mmol2meq(float, float);
	~CSM();

	int nNumOfDays;
	int nDepth;
	int thick;
	int nNumOfCompatments;
	int nArea;
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
	

	void InitCompartments();
	MONTH GetMonth(int nDay);
	float JULIAN(int day);


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	NumericVector RainArr;
	float AET;

	
	float nDustCa;
	float nDustSO4;
	float wieltingPoint;
	float nFieldCapacity;
	float CCa;
	float CSO4;
	float nCaco3less;
	
	
	float nTemp;
	float nLeachate;

	float TempArr[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
};

RCPP_EXPOSED_CLASS(CSM);


