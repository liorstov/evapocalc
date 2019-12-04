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

}


Rcpp::List CSM::Calculate(Rcpp::NumericVector rain, Rcpp::NumericVector PET, float years, float Depth, float nthick, float nwieltingPoint, float InitialCa,
 float initialSO4, float nBulkDensity, float FieldArea, float FieldCapacity, float DustCa, float DustSO4, float AETFactor)
{
	
	nNumOfDays = years*365;
	nDepth = Depth;
	thick = nthick;	
	nArea = FieldArea;
	nNumOfCompatments = nDepth / thick;
	wieltingPoint = nwieltingPoint * thick;
	CCa =meq2mmol(InitialCa, nthick*nArea); // mmol
	CSO4 = meq2mmol(initialSO4, nthick*nArea); // mmol
	BulkDensity =nBulkDensity; /// gr/cm^3
	nFieldCapacity = FieldCapacity * nthick;

	Rcpp::IntegerVector WD = Rcpp::IntegerVector::create();
	// carbonate in dust range from 0.5 to 5 [g m-2 yr-1]. 
	//I took 0.51 [g m-2 yr-1] = 5.1*10^-5 [g cm-2 yr-1] = 1.4*10^-7 [g cm-2 day-1]
	// then convert to [mmol cm-2 day-1] 
	nDailyDustCa =  (DustCa / (365.0F*10000.0F)) * 25.0F;
	nDailyDustSO4 = (DustSO4 / (365.0F*10000.0F)) * 10.0F;

	nTotalWhc = (nFieldCapacity - nwieltingPoint) * nNumOfCompatments;
	nTotalMoist = wieltingPoint * nNumOfCompatments;
	nTotalWP = nTotalMoist;
	InitCompartments();
	InitMonths();
	RainArr = rain;
	int nDeepestWetCompartment;
	int nMonth;
	float nTemp;
	float nDailyPET;

	// create list with all compartment
	Rcpp::DoubleVector vect = Rcpp::DoubleVector::create();
	Rcpp::DoubleVector moist = Rcpp::DoubleVector::create();
	vect.erase(0, vect.length());
	
	//main loop over days
	for (int day = 0; day < nNumOfDays; day++)
	{
		nDeepestWetCompartment = 0;
		nMonth = ((day % 365) / 31);
		nTemp = Months[nMonth].nTemp;
		nDailyPET = PET[day];

		/*nTotalCaDust += nDust;
		nTotalCaRain += RainArr[day] * CCa*40.0 / 1000.0;*/
		if (nTotalMoist >= (0.546*nTotalWhc)) // according to Marion et al. (1985), for the upper 45% of the total whc the actual evapotranspiration (AET) is the potential evapotranspiration (pet)
			AET = nDailyPET;		// in case of 10 compartments of 10 cm each, if total moistute > 8.465 
														//AET=PETdaily[monthperday[day]];
		else {										// the lower 55% of the total whc are according to modifeid Thornthwaite-Mather model
			AET = (nTotalMoist / nTotalWhc)*nDailyPET;
		}

		//Rcout << RainArr[day] << "  "<<AET << "  " << nDailyPET<< "  " << nTotalMoist << "   " << nTotalWP<<  endl;
		Months[nMonth].totalAET += AET;
		AET *= AETFactor;
		//AET = AET * 10;
		nTotalMoist = 0;
		nTotalAet += AET;
		
		nTotalRain += RainArr[day];
		Compartments[0].nMoist +=RainArr[day];   // set the moiture of the 1st compartment to the intial moisture plus the daily rainfall. rainfall is added only the 1st compartment

		//accumolate dust and relaese when its raining 
		if (RainArr[day] > 0) {
			Compartments[0].nCCa += accumolateDustDays * nDailyDustCa;
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
				Compartments[d].nMoist = Compartments[d].nThetaWeildingPnt;
			}
			else
			{
				Compartments[d].nMoist -= AET;
				AET = 0.0;
			}

			if (Compartments[d].nMoist > (Compartments[d].nWhc)) {
				// determines the leachate by substracting the field capacity from the moisture content
				nLeachate = Compartments[d].nMoist - (Compartments[d].nWhc);
				nDeepestWetCompartment = d;
			}
			// if there is no excess water, then the leachate is zero
			else  {
				nLeachate = 0.0;				
			}
			

			// subtract the leachate from the moisture of the  compartment
			Compartments[d].nMoist -= nLeachate; 

			//calculating gypsum concentration and ion available for washing
			// only if moist change since yesterday
			if (Compartments[d].nMoist != Compartments[d].nLastMoist) {
				Compartments[d].solubility(nTemp);
				Compartments[d].nLastMoist = Compartments[d].nMoist;
			}


			
			//start washing down
			//Rcout << "moist" << Compartments[d].nMoist << endl;
			//Rcout << "leachet "<< nLeachate << endl;
			//  examin if we reached saturation
			if (nLeachate > 0)
			{	
				//wash to next compartment or to leachete
				if (d != nNumOfCompatments - 1) {
					Compartments[d + 1].nCCa += Compartments[d].nCCa;
					Compartments[d + 1].C_SO4 += Compartments[d].C_SO4;
					Compartments[d + 1].nMoist += nLeachate;
				}
				else
				{
					nTotalCaLeachate += Compartments[d].nCCa;
					nTotalSO4Leachate += Compartments[d].C_SO4;
					nTotalLeachate += nLeachate;
				}

				Compartments[d].nCCa = 0;
				Compartments[d].C_SO4 = 0;
				nLeachate = 0;
			}		

			nTotalMoist += Compartments[d].nMoist;
		}

		WD.push_back(6);

		//Rcout << nDeepestWetCompartment << endl;

		if (day % 365 == 0)
		{
			for (std::vector<Compartment>::iterator it = Compartments.begin(); it != Compartments.end(); ++it) {
				//convert to meq/100 g soi;; first convert to mol with the moist and then to mmol and then multiply by 100/BDensity = 69
				vect.push_back(mmol2meq(it->C_CaSO4, (it->nthick * it->nArea)));
				moist.push_back(it->nTotMoist);
			}
		}

	}
	


	
	//create a list of months
	Rcpp::DoubleVector monthsVector = Rcpp::DoubleVector::create();
	for (std::vector<Month>::iterator it = Months.begin(); it != Months.end(); it++)
	{
		monthsVector.push_back(it->totalAET );
	}

	vect.attr("dim") = Dimension(nNumOfCompatments, years);
	moist.attr("dim") = Dimension(nNumOfCompatments, years);

	Rcpp::List returnList = Rcpp::List::create(_["gypsum"] = vect,
		_["month"] = monthsVector,
		_["WD"] = WD, _["moist"] = moist);

	return returnList;
}

std::vector<Compartment>* CSM::GetCompartments()
{
	return &Compartments;
}



//input meq/100g return mmol/cm3
float CSM::meq2mmol(float fMeq, float SoilVolume)
{
	//fmeq is meq/100g soil . multiply by bukl denisty[g/cm3] mult by 0.5[mmol/meq]
	return(fMeq * 1.44 * 0.01 *0.5 * SoilVolume);
}

//input mmol/cm3 return meq/100g soil
float CSM::mmol2meq(float mmol, float SoilVolume)
{
	return ((mmol * 2 /(SoilVolume*1.44))*  100  );
}

 

CSM::~CSM()
{
	//Compartments.~vector();
	//Months.~vector();
	printf ("csm destructor");
}

float CSM::GetPrecision(float x)
{
	return(((int)(x*10000.0)) / 10000.0F);
}

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
RCPP_MODULE(CSM_MODULE) {
	class_<CSM>("CSMCLASS")
		.constructor()

		.method("Calculate", &CSM::Calculate,
			"Docstring for stats")
		

		.method("InitMonths", &CSM::InitMonths, "desc")
		.field("RainArr", &CSM::RainArr, "rain array")
		;
}