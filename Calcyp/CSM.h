#pragma once


#include <stdio.h>
#include <vector>
#include "Compartment.h"
#include "Month.h"

//#include <Rcpp.h>
//using namespace Rcpp;

using namespace std;

// [[Rcpp::plugins(cpp11)]]

typedef enum { WINTER, SPRING, SUMMER } SEASON;


class CSM
{
public:
	CSM(float* rain);
	void Calculate();
	std::vector<Compartment>* GetCompartments();
	~CSM();

	int nNumOfDays;
	int nDepth;
	int thick;
	int nNumOfCompatments;
	int nArea;
	float nTotalCaDust;
	float nTotalCaRain;
	float nTotalCaLeachate;
	float nTotalSO4Leachate;
	float nTotalLeachate;
	float nTotalMoist;
	float nInitMoistTotal;
	float nTotalAet;
	float nTotalWP;
	float nTotalWhc;
	void InitMonths();
	
	

	void InitCompartments();
	MONTH GetMonth(int nDay);
	float JULIAN(int day);


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	float* RainArr;
	float* AET;

	
	float nDust;
	float wieltingPoint;
	float nFieldCapacity;
	float CCa;
	float CSO4;
	float nCaco3less;
	
	
	float nTemp;
	float nLeachate;

	float TempArr[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
};

//RCPP_MODULE(CSM_MODULE) {
//  class_<CSM>( "CSM" )
//  .constructor<float>()
//  
//  .method( "Calculate", &CSM::Calculate,
//  "Docstring for stats")
//  .method("InitMonths", &CSM::InitMonths, "desc")
//
//  ;}

