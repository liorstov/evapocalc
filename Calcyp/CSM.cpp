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
	nInitMoistTotal = 0;
	nTemp = 0;
	accumolateDustDays = 0.0F;
	
	nLeachate = 0;

	AET = 0;
	printf("csm Const\n");
}


vector<float> CSM::Calculate(vector<float> rain, vector<float> PET, float years, float Depth, float nthick, float WieltingPoint, float InitialCa,
 float initialSO4, float nBulkDensity, float FieldArea, float FieldCapacity, float DustCa, float DustSO4, float AETFactor)
{
	nNumOfDays = years*365;
	nDepth = Depth;
	thick = nthick;	
	nArea = FieldArea;
	nNumOfCompatments = nDepth / thick;
	wieltingPoint = WieltingPoint*thick ;//in cm
	BulkDensity =nBulkDensity; /// gr/cm^3
	nFieldCapacity = FieldCapacity*thick;
	RainArr = rain;
	// carbonate in dust range from 0.5 to 5 [g m-2 yr-1]. 
	//I took 0.51 [g m-2 yr-1] = 5.1*10^-5 [g cm-2 yr-1] = 1.4*10^-7 [g cm-2 day-1]
	// then convert to [mol cm-2 day-1] 
	nDailyDustCa =  (DustCa / (365.0F*10000.0F)) * 0.02495F; //calcium molar waight 0.02495
	nDailyDustSO4 = (DustSO4 / (365.0F*10000.0F)) * 0.0104F; //so4 molar waight 0.0104

	nTotalWhc = (nFieldCapacity - wieltingPoint) * nNumOfCompatments;//in cm
	nTotalMoist = wieltingPoint * nNumOfCompatments;//[cm]
	nInitMoistTotal = nTotalMoist;
	InitCompartments();
	InitMonths();

	CCa =meqSoil2molar(InitialCa, nthick*nArea, moistcm2Litre(nInitMoistTotal)); // molar
	CSO4 = meqSoil2molar(initialSO4, nthick*nArea, moistcm2Litre(nInitMoistTotal)); // molar
	
	int nMonth;
	float nTemp;
	float nDailyPET;
	float nDailyAET;
	float nDailyRain;
	int nWDComp = 0;
	int nRainEvents = 0;
	int nRainVecLength = RainArr.size();
	
	// create list with all compartment
	vector<float> vect;
	vector<float> Gypsum;
	vector<float> CompWash;
	vector<float> floodComp;
	vector<float> AETLoss;
	vector<float> dayRain;
	vector<float> WD;
	//vect.erase(0, vect.length());

	nNumOfDays =   std::min(nNumOfDays, nRainVecLength);

	//main loop over days
	for (int day = 0; ((day < nNumOfDays)); day++)
	{
		nMonth = ((day % 365) / 31);
		nTemp = 25;
		nDailyPET = PET[day]/10; 
		nDailyRain = RainArr[day];
		/*nTotalCaDust += nDust;
		nTotalCaRain += RainArr[day] * CCa*40.0 / 1000.0;*/
		if (nTotalMoist >= (0.546*nTotalWhc)) // according to Marion et al. (1985), for the upper 45% of the total whc the actual evapotranspiration (AET) is the potential evapotranspiration (pet)
			AET = nDailyPET;		// in case of 10 compartments of 10 cm each, if total moistute > 8.465 
														//AET=PETdaily[monthperday[day]];
		else {										// the lower 55% of the total whc are according to modifeid Thornthwaite-Mather model
			AET = (nTotalMoist / nTotalWhc)*nDailyPET;
		}

		nTotalMoist = 0;
		AET *= AETFactor;
		
		Months[nMonth].totalAET += AET;
		nDailyAET = AET;
		//AET = AET * 10;
		nTotalAet += AET;
		
		nTotalRain += nDailyRain;
		Compartments[0].nMoist += nDailyRain;   // set the moiture of the 1st compartment to the intial moisture plus the daily rainfall. rainfall is added only the 1st compartment

		//accumolate dust and relaese when its raining 
		if (nDailyRain > 0) {
			Compartments[0].C_Ca += accumolateDustDays * nDailyDustCa;
			Compartments[0].C_SO4 += accumolateDustDays * nDailyDustSO4;
			accumolateDustDays = 0;
		}
		else
		{
			accumolateDustDays++;
		}

		// This is the second loop that runs through the soil profile
		for (int d = 0; d < nNumOfCompatments; d++)
		{

			//  taking into account the AET for this current compartment, and updating the AET value
			// moist - wieltingPOint is the water available for evaporation
			if (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt) < AET)
			{
				AET = AET - (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt));
				Compartments[d].fAETLoss += Compartments[d].nMoist - Compartments[d].nThetaWeildingPnt;
				Compartments[d].nMoist = Compartments[d].nThetaWeildingPnt;
			}
			else
			{
				Compartments[d].nMoist -= AET; 
				Compartments[d].fAETLoss += AET;
				AET = 0.0;
			}

			
			if (Compartments[d].nMoist > (Compartments[d].nWhc)) 
			{				
				// determines the leachate by substracting the field capacity from the moisture content
				nLeachate = Compartments[d].nMoist - (Compartments[d].nWhc);
				Compartments[d].nFloodedCount++;
				nWDComp = d+1;
			}
			// in case the comp is floated without leaching
			else if(Compartments[d].nMoist == (Compartments[d].nWhc))			{			
				
				Compartments[d].nFloodedCount++;
				nLeachate = 0.0;				
			}
			else
			{
				nLeachate = 0.0;
			}
			
			Compartments[d].fTotLeachate += nLeachate;

			// subtract the leachate from the moisture of the  compartment
			Compartments[d].nMoist -= nLeachate; 

			//calculating gypsum concentration and ion available for washing
			// only if moist change since yesterday
			if (Compartments[d].nMoist != Compartments[d].nLastMoist) {
				Compartments[d].solubility(nTemp);
				Compartments[d].nLastMoist = Compartments[d].nMoist;
			}
			
			//start washing down
			//  examin if we reached saturation
			if (nLeachate > 0)
			{	
				//wash to next compartment or to leachete
				// ion are in mol
				if (d != nNumOfCompatments - 1) {
					Compartments[d + 1].C_Ca += Compartments[d].C_Ca;
					Compartments[d + 1].C_SO4 += Compartments[d].C_SO4;
					Compartments[d + 1].nMoist += nLeachate;
				}
				else
				{
					nTotalCaLeachate += Compartments[d].C_Ca;
					nTotalSO4Leachate += Compartments[d].C_SO4;
					nTotalLeachate += nLeachate;
					nWDComp = d;
				}

				Compartments[d].C_Ca = 0;
				Compartments[d].C_SO4 = 0;
				nLeachate = 0;
			}		

			nTotalMoist += Compartments[d].nMoist;
		}

		
		for (std::vector<Compartment>::iterator it = Compartments.begin(); it != Compartments.end(); ++it) {
			//convert to meq/100 g soi;; first convert to mol with the moist and then to mmol and then multiply by 100/BDensity = 69
			WD.push_back(it->nMoist);

		}

		if (nDailyRain > 0 )
		{
			nRainEvents++;
			Compartments[nWDComp].nWetCount++;

			WD.push_back(day);		
		WD.push_back(nDailyRain);	
		WD.push_back(nDailyAET);
		WD.push_back((nWDComp + 0.5)*nthick);	
		WD.push_back(nTotalMoist);
		}
		nWDComp = 0;

		
		
		
	}

	for (std::vector<Compartment>::iterator it = Compartments.begin(); it != Compartments.end(); ++it) {
		//convert to meq/100 g soi;; first convert to mol with the moist and then to mmol and then multiply by 100/BDensity = 69
		vect.push_back(it->nWetCount);
		CompWash.push_back(it->fTotLeachate);
		AETLoss.push_back(it->fAETLoss); 
		floodComp.push_back(it->nFloodedCount);		
		Gypsum.push_back(molar2meqSoil(it->C_CaSO4, 5, moistcm2Litre(it->nMoist)));

	}
	
	
	/*Rcout << "total moist in soil: " << nTotalMoist << endl <<
		"total leachate: " << nLeachate << endl <<
		"total rain: " << nTotalRain << endl <<
		"total AET: " << nTotalAet << endl <<
		"balance: " << nTotalRain - nLeachate - nTotalMoist - nInitMoistTotal - nTotalAet << endl <<
		"rain events: " << nRainEvents <<  endl;*/
	/*Gypsum.attr("dim") = Dimension(20, Gypsum.length() / (20));
	NumericMatrix GypsumMat = transpose(as<NumericMatrix>(Gypsum));

	WD.attr("dim") = Dimension(5, WD.length() / (5));
	NumericMatrix WDMat = transpose(as<NumericMatrix>(WD));
	colnames(WDMat) = CharacterVector::create("day","rain","AET", "WD","soilMoisture");
	results = Rcpp::List::create(_["gypsum"] = vect,	
		_["WD"] = WDMat,
		_["WetZone"] = floodComp,
		_["AET"] = nTotalAet,
		_["AETLoss"] = AETLoss,
		_["GypsumDay"] = GypsumMat,
		_["totalRain"] = nTotalRain,
		_["CompWash"] = CompWash,
		_["Leachate"]= nTotalLeachate);*/

	return Gypsum;
}

std::vector<Compartment>* CSM::GetCompartments()
{
	return &Compartments;
}



//input meq/100g return mol/Litre
float CSM::meqSoil2molar(float fMeq, float SoilVolume, float moisture)
{
	//fmeq is meq/100g soil . multiply by bukl denisty[g/cm3] divide by 100[gr soil ] mult by soilVolume [cm3] mult by 0.5[mmol/meq] mult by 0.001 [mol/mmol] divide by moisture
	return(fMeq * 1.44* SoilVolume/100  *0.5* 0.001 / moisture );
}

//input mol/Litre mult by moisture [litre] mult by 0.002 [meq/mol] divide by volume return meq/100g soil
float CSM::molar2meqSoil(float molar, float SoilVolume, float moisture)
{
	return ((molar * moisture * 0.002 /(SoilVolume*1.44))*  100  );
}

float CSM::moistcm2Litre(float moist_cm) {
	return(moist_cm*0.001F);
}

 

CSM::~CSM()
{
	//Compartments.~vector();
	//Months.~vector();
	printf ("csm destructor\n");
}

float CSM::GetPrecision(float x)
{
	return(((int)(x*10000.0)) / 10000.0F);
}

//Rcpp::List CSM::GetResults()
//{
//	return results;
//}

void CSM::InitCompartments()
{
	Compartment *newCompartment;
	for (int i = 0; i < nNumOfCompatments; i++)
	{
		newCompartment = new Compartment(i, nArea, wieltingPoint, nFieldCapacity, thick, CCa, CSO4);
		Compartments.push_back(*newCompartment);
	}
}

void CSM::InitMonths()
{
  //NumericVector returnList = NumericVector::create(1,2,3,4,5,6,7,8,9,10,11,12);
  float II = 0.0;
	float a = 0;
	float PET = 0;
	float PETdaily = 0;
	float pan = 0;
	float panDaily = 0;

	// calculate heat index
	for (int i = 0; i < 12; i++)
	{		
		II = II + pow((TempArr[i] / 5.0), 1.514);
	}

	a = (6.75*pow(10.0, -7.0)*pow(II, 3.0)) - (7.71*pow(10.0, -5.0)*pow(II, 2.0)) + (0.01792*II) + 0.49239;

	for (int i = 0; i < 12; i++)
	{
		// calculates PET [cm/month] according to Thornthwaite (1948) and PEV [cm/month] according to equations in figure 3 of Marion et al (1985)
		PET = 1.6*pow((10.0*TempArr[i] / II), a);
		PETdaily = PET / 30;

		
		if (i < 7)
		{
			pan = ((17.07*exp(-0.1309*TempArr[i])) + 1.91)*PET;
			panDaily = ((17.07*exp(-0.1309*TempArr[i])) + 1.91)*PETdaily;
		}
		else
		{
			pan = ((13.67*exp(-0.13*TempArr[i])) + 1.35)*PET;
			panDaily = ((13.67*exp(-0.13*TempArr[i])) + 1.35)*PETdaily;
		}

		Months.push_back(Month(i, TempArr[i], panDaily, PETdaily, PET, pan));
		//returnList[i] = 20;
	}
	
	//return returnList;

	


}

MONTH CSM::GetMonth(int nDay)
{
	MONTH month = Jan;
	nDay = JULIAN(nDay);
	if (nDay <= 31)
		month = Oct;
	else
		if (nDay <= 59)
			month = Nov;
		else
			if (nDay <= 90)
				month =Dec;
			else
				if (nDay <= 120)
					month = Jan;
				else
					if (nDay <= 151)
						month = Feb;
					else
						if (nDay <= 181)
							month = Apr;
						else
							if (nDay <= 212)
								month = May;
							else
								if (nDay <= 243)
									month = Jun;
								else
									if (nDay <= 273)
										month = Jul;
									else
										if (nDay <= 304)
											month = Aug;
										else
											if (nDay <= 334)
												month = Sep;
											else
												if (nDay <= 365)
													month = Oct;
	return month;

}

float CSM::JULIAN(int day)

{
	int julian;

	julian = day - ((long)(day / 365) * 365);
	return julian;
}

// [[Rcpp::plugins(cpp11)]]
//RCPP_MODULE(CSM_MODULE) {
//	using namespace Rcpp;
//	class_<CSM>("CSMCLASS")
//	.constructor()
//
//	.method("Calculate", &CSM::Calculate,
//		"Docstring for stats")
//		
//
//	.method("InitMonths", &CSM::InitMonths, "desc")
//	.method("GetResults", &CSM::GetResults, "GetResults")
//	.field("RainArr", &CSM::RainArr, "rain array")
//	;
//}