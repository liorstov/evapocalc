#include "CSM.h"
#include <ctime>


#define MM_Ca = 40.08
#define MM_SO4 = 96.06
#define MM_CaSO4 = 136.14;




//This code calculate the content and depth of pedogenic carbonate with depth
CSM::CSM()
{
	nTotalWhc = 0;
	nTotalCaDust = 0;
	nTotalRain = 0;
	nTotalCaLeachate = 0;
	nTotalSO4Leachate = 0;
	nTotalLeachate = 0;
	nTotalMoist = 0;
	nTotalAet = 0;
	nTotalWP = 0;
	nTemp = 0;
	accumolateDustDays = 0.0F;
	
	nLeachate = 0;

	AET = 0;
	printf("csm Const\n");
}

//flux is gram/cm2/day
//dust concentration is mg/gram; rain concentration is mg/l
Rcpp::List CSM::Calculate(Rcpp::DoubleVector rain, Rcpp::DoubleVector PET, int years, int Depth, int nthick, double WieltingPoint,
	int FieldArea, double FieldCapacity, double DustGyp, double AETFactor, bool verbose, double dustFlux, double rainCa, double rainSO4)
{
	Rcpp::Rcout << "years: " << years << endl <<
		"Depth: " << Depth << endl <<
		"thick: " << nthick << endl <<
		"WieltingPoint: " << WieltingPoint << endl <<
		"FieldArea: " << FieldArea << endl <<
		"FieldCapacity: " << FieldCapacity <<  endl <<
		"DustGyp: "  << DustGyp<< endl << "AETFactor: " << AETFactor << endl << "verbose: " << verbose << endl
		<< "dustFlux: " << dustFlux<< endl<< "rainCa: " << rainCa << endl << "rainSO4: " << rainSO4<< endl;

	double nTemp = 25;
	double nDailyPET = 0;
	double nDailyAET = 0;
	double nDailyRain = 0;
	double DailySO4Rain;
	double DailyCaRain;
	double inputCa = 0;
	double outputCa = 0;
	
	int nWDComp = 0;
	int nRainEvents = 0;
	int nRainVecLength = rain.length();

	nTotalWhc = 0;
	nTotalRain = 0;
	nTotalCaLeachate = 0;
	nTotalSO4Leachate = 0;
	nTotalLeachate = 0;
	nTotalMoist = 0;
	nTotalAet = 0;
	nTotalWP = 0;
	accumolateDustDays = 0.0F;

	nLeachate = 0;

	nNumOfDays= years*365;
	nDepth = Depth;
	thick = nthick;	
	nArea = FieldArea;
	nNumOfCompatments = nDepth / thick;
	wieltingPoint = WieltingPoint*thick ;//in cm
	nFieldCapacity = FieldCapacity*thick;
	RainArr = rain;

	// flux is gram/cm2/day multiple by concentration % convert to mol/cm2/day	
	nDailyDustGyp = dustFlux * DustGyp  / 172.172F; //gypsum molar waight  172.172
	double WHC = nFieldCapacity - wieltingPoint;
	nTotalWhc = (nFieldCapacity - wieltingPoint) * nNumOfCompatments;//in cm
	nTotalMoist = wieltingPoint * nNumOfCompatments;//[cm]
	nTotalWP = nTotalMoist;
	InitCompartments();

	CCa =0; // molar
	CSO4 = 0; // molar
	
	

	// create list with all compartment
	Rcpp::DoubleVector vect = Rcpp::DoubleVector::create(); 
	Rcpp::DoubleVector Gypsum = Rcpp::DoubleVector::create(); 
	Rcpp::DoubleVector GypsumDay = Rcpp::DoubleVector::create(); 
	Rcpp::DoubleVector Ca = Rcpp::DoubleVector::create();
	Rcpp::DoubleVector SO4 = Rcpp::DoubleVector::create();
	Rcpp::DoubleVector moisture = Rcpp::DoubleVector::create();
	Rcpp::DoubleVector CompWash = Rcpp::DoubleVector::create();
	Rcpp::DoubleVector floodComp = DoubleVector::create();
	Rcpp::DoubleVector AETLoss = DoubleVector::create();
	Rcpp::DoubleVector dayRain = DoubleVector::create();
	Rcpp::DoubleVector WD = DoubleVector::create();
	Gypsum.erase(0, Gypsum.length());
	
	Rcout.precision(10);
	/*std::ofstream myfile;
	myfile.open("example.csv");
	myfile << "day,compartment, daily rain [L], dialy AET[L],moistutre [L], input Ca [mol/L], input SO4 [mol/L], input gyp [mol/L], output Ca [mol/L], output SO4 [mol/L], output gyp [mol/L].\n";
*/
	nNumOfDays =   std::min(nNumOfDays, nRainVecLength);
	Rcpp::Rcout << nNumOfDays << endl;
	//main loop over days
	for (int day = 0; ((day < nNumOfDays)); day++)
	{
		nDailyPET = PET[day]/10; 
		nDailyRain = RainArr[day]/10;
		DailySO4Rain = nDailyRain * 0.001*rainSO4 / 1000 / 96.06F;// convert cm3 to littre and multiple with concentration
		DailyCaRain = nDailyRain * 0.001*rainCa / 1000 / 40.078F;// convert cm3 to littre and multiple with concentration

		inputCa += DailyCaRain;
		/*nTotalCaDust += nDust;
		nTotalCaRain += RainArr[day] * CCa*40.0 / 1000.0;*/
		if (nTotalMoist >= (0.546*nTotalWhc)) // according to Marion et al. (1985), for the upper 45% of the total whc the actual evapotranspiration (AET) is the potential evapotranspiration (pet)
			AET = nDailyPET;		// in case of 10 compartments of 10 cm each, if total moistute > 8.465 
														//AET=PETdaily[monthperday[day]];
		else if ((nTotalMoist > nTotalWP*1.001) && (nTotalMoist < 0.546*nTotalWhc)) {										// the lower 55% of the total whc are according to modifeid Thornthwaite-Mather model
			AET = (nTotalMoist / nTotalWhc)*nDailyPET;
		}
		else{// in case the soil is at WP, no evaporation
			AET = 0;
		}
		//if (AET == 0 && nDailyRain == 0) continue;
		//Rcpp::Rcout << nTotalMoist << endl << nTotalWP << endl << nTotalWhc<<endl<<AET<<endl;
		nTotalMoist = 0;
		AET *= AETFactor;
		
		nDailyAET = AET;
		//AET = AET * 10;
		nTotalAet += nDailyAET;
		
		nTotalRain += nDailyRain;
		Compartments[0].nMoist += nDailyRain;   // set the moiture of the 1st compartment to the intial moisture plus the daily rainfall. rainfall is added only the 1st compartment

		//accumolate dust and relaese when its raining 
		if (nDailyRain > 0) {
			Compartments[0].C_Ca +=  DailyCaRain;
			Compartments[0].C_SO4 += DailySO4Rain;
			Compartments[0].C_CaSO4 += accumolateDustDays * nDailyDustGyp;
			accumolateDustDays = 0;
		}
		else
		{
			accumolateDustDays++;
		}

		// This is the second loop that runs through the soil profile
		for (int CurrentComp = 0; CurrentComp < nNumOfCompatments; CurrentComp++)
		{

			//WASHING
			if (Compartments[CurrentComp].nMoist > (Compartments[CurrentComp].nWhc)) 
			{				
				// determines the leachate by substracting the field capacity from the moisture content
				nLeachate = Compartments[CurrentComp].nMoist - (Compartments[CurrentComp].nWhc);
				Compartments[CurrentComp].nFloodedCount++;
				nWDComp = CurrentComp+1;
			}
			// in case the comp is floated without leaching
			else if(Compartments[CurrentComp].nMoist == (Compartments[CurrentComp].nWhc))			{			
				
				Compartments[CurrentComp].nFloodedCount++;
				nLeachate = 0.0;				
			}
			else
			{
				nLeachate = 0.0;
			}
			
			Compartments[CurrentComp].fTotLeachate += nLeachate;

			// subtract the leachate from the moisture of the  compartment
			Compartments[CurrentComp].nMoist -= nLeachate; 


			// EVAPORATING
			//  taking into account the AET for this current compartment, and updating the AET value
			// moist - wieltingPOint is the water available for evaporation
			if (Compartments[CurrentComp].nMoist - (Compartments[CurrentComp].nThetaWeildingPnt) < AET)
			{
				AET -= (Compartments[CurrentComp].nMoist - (Compartments[CurrentComp].nThetaWeildingPnt));
				Compartments[CurrentComp].fAETLoss += Compartments[CurrentComp].nMoist - Compartments[CurrentComp].nThetaWeildingPnt;
				Compartments[CurrentComp].nMoist = Compartments[CurrentComp].nThetaWeildingPnt;
			}
			else
			{
				Compartments[CurrentComp].nMoist -= AET; 
				Compartments[CurrentComp].fAETLoss += AET;
				AET = 0.0;
			}

			
			//CHEMISTRY
			//calculating gypsum concentration and ion available for washing
			// only if moist change since yesterday
			if (Compartments[CurrentComp].nMoist != Compartments[CurrentComp].nLastMoist) {
				Compartments[CurrentComp].nLastMoist = Compartments[CurrentComp].nMoist;
				
		
			}
				//	myfile << day << "," << CurrentComp << "," << nDailyRain * 0.001 << "," << nDailyAET * 0.001 << "," << Compartments[CurrentComp].nMoist * 0.001 << "," << Compartments[CurrentComp].C_Ca / (Compartments[CurrentComp].nMoist * 0.001F) << "," << Compartments[CurrentComp].C_SO4 / (Compartments[CurrentComp].nMoist * 0.001F) << "," << Compartments[CurrentComp].C_CaSO4 / (Compartments[CurrentComp].nMoist * 0.001F) <<",";

			Compartments[CurrentComp].solubility(nTemp);
		//	myfile << Compartments[CurrentComp].C_Ca / (Compartments[CurrentComp].nMoist * 0.001F) << "," << Compartments[CurrentComp].C_SO4 / (Compartments[CurrentComp].nMoist * 0.001F) << "," << Compartments[CurrentComp].C_CaSO4 / (Compartments[CurrentComp].nMoist * 0.001F) << endl;

			
			//LEACHATE
			//  examin if we reached saturation
			if (nLeachate > 0)
			{	
				//wash to next compartment or to leachete
				// ion are in mol
				if (CurrentComp != nNumOfCompatments - 1) {
					//adding the fractional quantitiy of matter
					Compartments[CurrentComp + 1].C_Ca += Compartments[CurrentComp].C_Ca* nLeachate/(WHC+nLeachate);
					Compartments[CurrentComp + 1].C_SO4 += Compartments[CurrentComp].C_SO4* nLeachate / (WHC + nLeachate);
					Compartments[CurrentComp + 1].nMoist += nLeachate;
				}
				else
				{
					//adding the fractional quantitiy of matter
					nTotalCaLeachate += Compartments[CurrentComp].C_Ca* nLeachate / (WHC + nLeachate);
					nTotalSO4Leachate += Compartments[CurrentComp].C_SO4* nLeachate / (WHC + nLeachate);
					nTotalLeachate += nLeachate;
					nWDComp = CurrentComp;
				}
				//leaving  the fractional quantitiy of matter
				Compartments[CurrentComp].C_Ca *= WHC/ (WHC + nLeachate);
				Compartments[CurrentComp].C_SO4 *= WHC / (WHC + nLeachate);
				nLeachate = 0;
			}		

			nTotalMoist += Compartments[CurrentComp].nMoist;
		}

		// Write to output
		
		if (verbose == 1)
		{
			GypsumDay.push_back(day);
			GypsumDay.push_back(nDailyRain);
			GypsumDay.push_back(nDailyAET);
			GypsumDay.push_back(nTotalMoist);
			GypsumDay.push_back((nWDComp + 0.5)*nthick);

			
			WD.push_back(day);		
			WD.push_back(nDailyRain);	
			WD.push_back(nDailyAET);			
			WD.push_back(nTotalMoist);
			WD.push_back((nWDComp + 0.5)*nthick);
			for (std::vector<Compartment>::iterator it = Compartments.begin(); (it - Compartments.begin()) < 15 ; ++it) {
				WD.push_back(it->nMoist);			
				GypsumDay.push_back(it->C_CaSO4);
			}
			
		}
		else if (nDailyRain > 0 )
		{
			nRainEvents++;
			Compartments[nWDComp].nWetCount++;	
			WD.push_back(day);
			WD.push_back(nDailyRain);
			WD.push_back((nWDComp + 0.5)*nthick);
		
		}
		
		nWDComp = 0;
		
		if (fmod(double(day)/ double(nNumOfDays), 0.1 )  == 0) {
			Rcout << double(day) / double(nNumOfDays) * 100 << "%; ";
		}
	}

	for (std::vector<Compartment>::iterator it = Compartments.begin(); it != Compartments.end(); ++it) {
		//convert to meq/100 g soi;; first convert to mol with the moist and then to mmol and then multiply by 100/BDensity = 69
		vect.push_back(it->nWetCount);
		CompWash.push_back(it->fTotLeachate);
		AETLoss.push_back(it->fAETLoss); 
		floodComp.push_back(it->nFloodedCount);
		Gypsum.push_back(mol2meqSoil(it->C_CaSO4,thick));
		Ca.push_back(mol2meqSoil(it->C_Ca, thick));
		SO4.push_back(mol2meqSoil(it->C_SO4, thick));

		moisture.push_back(it->nMoist);

		outputCa += it->C_Ca + it->C_CaSO4;
	}

	outputCa += nTotalCaLeachate;
	Rcpp::Rcout << WD.length() << verbose<< endl;

	//myfile.close();
	
	/*Rcout << "total moist in soil: " << nTotalMoist << endl <<
		"total leachate: " << nLeachate << endl <<
		"total rain: " << nTotalRain << endl <<
		"total AET: " << nTotalAet << endl <<
		"balance: " << nTotalRain - nLeachate - nTotalMoist - nInitMoistTotal - nTotalAet << endl <<
		"rain events: " << nRainEvents <<  endl;*/
	/*Gypsum.attr("dim") = Dimension(20, Gypsum.length() / (20));
	NumericMatrix GypsumMat = transpose(as<NumericMatrix>(Gypsum));*/

	
	results = Rcpp::List::create(_["gypsum"] = Gypsum, _["gypsumDay"] = output2Matrix(GypsumDay, verbose),
		_["Ca"] = Ca,
		_["SO4"] = SO4,
		_["WD"] = output2Matrix(WD, verbose),
		_["moist"] = moisture,
		_["WetZone"] = floodComp,
		_["AET"] = nTotalAet,
		_["AETLoss"] = AETLoss,		
		_["totalRain"] = nTotalRain,
		_["CompWash"] = CompWash,
		_["Leachate"]= nTotalLeachate,
		_["inputCa"]= inputCa,
		_["outputCa"]= outputCa);

	return results;
}

std::vector<Compartment>* CSM::GetCompartments()
{
	return &Compartments;
}

Rcpp::NumericMatrix CSM::output2Matrix(Rcpp::DoubleVector & inputVector, bool verbose) 
{
	Rcpp::NumericMatrix outputMat;
	if (verbose == 1)
	{
		inputVector.attr("dim") = Dimension(20, inputVector.length() / (20));
		outputMat = transpose(as<NumericMatrix>(inputVector));
		Rcpp::colnames(outputMat) = CharacterVector::create("day", "rain[cm]", "AET", "soilMoisture",  "WD","1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15");
	}
	else
	{
		inputVector.attr("dim") = Dimension(3, inputVector.length() / (3));
		outputMat = transpose(as<NumericMatrix>(inputVector));
		Rcpp::colnames(outputMat) = CharacterVector::create("day","rain[cm]", "WD");
	}
	Rcpp::Rcout << outputMat.nrow() << " " << outputMat.ncol() << endl;
	return(outputMat);
}

//input meq/100g return mol/Litre
double CSM::meqSoil2molar(double fMeq, double SoilVolume, double moisture)
{
	//fmeq is meq/100g soil . multiply by bukl denisty[g/cm3] divide by 100[gr soil ] mult by soilVolume [cm3] mult by 0.5[mmol/meq] mult by 0.001 [mol/mmol] divide by moisture
	return(fMeq * 1.44* SoilVolume/100  *0.5* 0.001 / moisture );
}

//input mol and multiply by 2000 [meq/mol] divide by compartment soil wight and multiply to represent 100 g soil
double CSM::mol2meqSoil(double mol, double SoilVolume)
{
	double inMeq = mol * 2000;
	return (( inMeq * 100/(SoilVolume*1.44))  );
}

double CSM::moistcm2Litre(double moist_cm) {
	return(moist_cm*0.001F);
}

 

CSM::~CSM()
{
	//Compartments.~vector();
	//Months.~vector();
	printf ("csm destructor\n");
}

double CSM::GetPrecision(double x)
{
	return(((int)(x*10000.0)) / 10000.0F);
}

//Rcpp::List CSM::GetResults()
//{
//	return results;
//}

void CSM::InitCompartments()
{
	if (Compartments.size() > 0) {
		Compartments.clear();
	}
	Compartment *newCompartment;
	for (int i = 0; i < nNumOfCompatments; i++)
	{
		newCompartment = new Compartment(i, nArea, wieltingPoint, nFieldCapacity, thick, CCa, CSO4);
		Compartments.push_back(*newCompartment);
	}
}


	

	







// [[Rcpp::plugins(cpp11)]]
RCPP_MODULE(CSM_MODULE) {
	using namespace Rcpp;
	class_<CSM>("CSMCLASS")
	.constructor()

	.method("Calculate", &CSM::Calculate,
		"Docstring for stats")
		

	.field("RainArr", &CSM::RainArr, "rain array")
	;
}