#pragma once

#include <Rcpp.h>
#include <stdio.h>
#include <vector>
#include "Compartment.h"
#include "Month.h"


using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

typedef enum { WINTER, SPRING, SUMMER } SEASON;


class CSM
{
public:
	CSM(float rain);
	void Calculate();
	std::vector<Compartment>* GetCompartments();
	~CSM();

	int nNumOfDays;
	int nDepth;
	int thick;
	int nArea;
	float nTotalCaDust;
	float nTotalCaRain;
	float nTotalCaLeachate;
	float nTotalLeachate;
	float nTotalMoist;
	float nInitMoistTotal;
	float nTotalAet;
	float nTotalWP;
	float nTotalWhc;
	NumericVector InitMonths();
	
	

	void InitCompartments();
	MONTH GetMonth(int nDay);
	float JULIAN(int day);


	std::vector<Compartment> Compartments;
	std::vector<Month> Months;
	float* RainArr;
	float* AET;
	float* TempArr;

	
	float nDust;
	float wieltingPoint;
	float nFieldCapacity;
	float CCa;
	float nCaco3less;
	
	float nTemp;
	float nLeachate;

	double CO2[10] = { 3.79e-4,6.87e-4,9.76e-4,1.28e-3,1.60e-3, 1.91e-3,2.14e-3,2.27e-3,2.41e-3,2.54e-3  };
	double CO2Mul[12] = { 1,1,1.588,1.588,1.588,1.588,1.588,1.588,1.588,1.588,1,1 };
};

RCPP_MODULE(CSM_MODULE) {
  class_<CSM>( "CSM" )
  .constructor<float>()
  
  .method( "Calculate", &CSM::Calculate,
  "Docstring for stats")
  .method("InitMonths", &CSM::InitMonths, "desc")

  ;}

