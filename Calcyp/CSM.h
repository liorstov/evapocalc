#pragma once
#ifndef _CSM_H
#define _CSM_H

#include <vector>
#include <stack>
#include <stdio.h>
#include "Compartment.h"
#include <iostream>
#include <fstream>
#include <Rcpp.h>


using namespace Rcpp;

using namespace std;




class CSM
{
public:
	CSM();
	Rcpp::List Calculate(Rcpp::DoubleVector, Rcpp::DoubleVector, int years, int Depth, int nthick,
	double WieltingPoint,int FieldArea, double FieldCapacity, double DustSO4, double AETFactor, bool verbose, double dustFlux, double rainCa, double rainSO4,bool withFC);
	std::vector<Compartment>* GetCompartments();
	Rcpp::NumericMatrix output2Matrix(Rcpp::DoubleVector & inputVector, bool verbose);
	double meqSoil2molar(double, double,double);
	double mol2meqSoil(double mol, double);
	double moistcm2Litre(double);

	~CSM();
	double GetPrecision(double);
	//Rcpp::List GetResults();
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
	float nTotalWP;
	float nTotalAet;
	float nTotalWhc;
	float BulkDensity;
	float accumolateDustDays;
	
	Rcpp::List results;
	void InitCompartments();
	void initVector(Rcpp::DoubleVector & inputVector);
	inline bool firstDayInYear(int day) {
		return(!((day + 1) % 365));
	}

//holocene profile have 0.1 FC and plesitocene profiles have 0.19 FC. a linear function updates the fc every year
	inline int updateFieldCapacity(int year) {
		if (year > 10000) {
			int CompDuration = 12500;
			int DustComp = (year-10001)/ (int)CompDuration;
			float yearsInProcess = (year - 10001) % CompDuration;
			//x is yearsinprocess and 0.01 is b and 
			Compartments[DustComp].nFieldCapacity = ((yearsInProcess)*0.09/ (float)CompDuration +0.1)*thick;
			
			updateTotalWHC();

			//Rcout << year <<"   " <<DustComp <<"  "<< Compartments[DustComp].nFieldCapacity/thick<<"   WHC" << nTotalWhc << endl;
			return(DustComp);
		}		
	}

	inline void updateTotalWHC() {
		nTotalWhc = 0.0;
		for (std::vector<Compartment>::iterator it = Compartments.begin(); it != Compartments.end(); ++it) {
			nTotalWhc += (it->nFieldCapacity - it->nThetaWeildingPnt);
		}

	}

	// vector of soil profile compartments
	std::vector<Compartment> Compartments;
	Rcpp::NumericVector RainArr;
	float AET;

	
	float nDailyDustCa;
	float nDailyDustGyp;
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