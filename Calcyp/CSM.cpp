#include "CSM.h"
#include <ctime>

#define MM_Ca = 40.08
#define MM_SO4 = 96.06
#define MM_CaSO4 = 136.14;



//This code calculate the content and depth of pedogenic carbonate with depth
CSM::CSM(float * rain)
{
	nNumOfDays = 1 * 365;
	nDepth = 30;
	thick = 5;
	nNumOfCompatments = nDepth / thick;
	wieltingPoint = 0.039;
	CCa = 0.864 * pow(10, -5); // mol
	CSO4 = 1.04 * pow(10, -4); // mol
	BulkDensity = 1.44; /// gr/cm^3
	nArea = 1;
	nFieldCapacity = 0.2;
	nDust = 0.51 / (365.0*10000.0);
	
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
	
	
	nLeachate = 0;
	RainArr = rain;

	AET = 0;

}

void CSM::Calculate()
{
	InitCompartments();
 
  
	InitMonths();
	printf ("Characters: %c %c \n", 'a', 65);
	int nMonth;
	float nTemp;
	//main loop over days
	for (int day = 0; day < nNumOfDays; day++)
	{
		nMonth = ((day % 365) / 31);
		nTemp = Months[nMonth].nTemp;
		 
		/*nTotalCaDust += nDust;
		nTotalCaRain += RainArr[day] * CCa*40.0 / 1000.0;*/
		if (nTotalMoist >= ((0.546*nTotalWhc) + nTotalWP)) // according to Marion et al. (1985), for the upper 45% of the total whc the actual evapotranspiration (AET) is the potential evapotranspiration (pet)
			AET = Months[(int)GetMonth(day)].PanDaily;		// in case of 10 compartments of 10 cm each, if total moistute > 8.465 
														//AET=PETdaily[monthperday[day]];
		else {										// the lower 55% of the total whc are according to modifeid Thornthwaite-Mather model
			AET = (nTotalMoist - nTotalWP) / (0.546*nTotalWhc)*Months[(int)GetMonth(day)].PanDaily;
		}

		nTotalMoist = 0;
		nTotalAet += AET;
		
		nTotalRain += RainArr[day];
		Compartments[0].nMoist +=RainArr[day];   // set the moiture of the 1st compartment to the intial moisture plus the daily rainfall. rainfall is added only the 1st compartment

		// This is the second loop that runs through the soil profile, from the 2nd compartment to the bottom
		for (int d = 0; d < nNumOfCompatments; d++)
		{

			//  taking into account the AET for this current compartment, and updating the AET value
			// moist - wieltingPOint is the water available for evaporation
			if (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt * thick) < AET)
			{
				AET = AET - (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt * thick));
				Compartments[d].nMoist = Compartments[d].nThetaWeildingPnt * thick;
			}
			else
			{
				Compartments[d].nMoist -= AET;
				AET = 0.0;
			}

			if (Compartments[d].nMoist > (Compartments[d].nWhc*thick)) {
				// determines the leachate by substracting the field capacity from the moisture content
				nLeachate = Compartments[d].nMoist - (Compartments[d].nWhc)*thick;
			}
			// if there is no excess water, then the leachate is zero
			else nLeachate = 0.0;

			// subtract the leachate from the moisture of the  compartment
			Compartments[d].nMoist -= nLeachate; 

			//calculating gypsum concentration and ion available for washing
			Compartments[d].solubility(nTemp);

			
			//start washing down
			
			//  examin if we reached saturation
			if (nLeachate > 0)
			{	
				//wash to next compartment or to leachete
				if (d != nNumOfCompatments - 2) {
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

			
			printf("%ld\t%0.3f\t%ld\t%0.3f\t%0.3f\n", day,RainArr[day],d, Compartments[d].nMoist, Compartments[d].C_CaSO4);

			nTotalMoist += Compartments[d].nMoist;
		}

	}

	// aggreagating total moist
	
}

std::vector<Compartment>* CSM::GetCompartments()
{
	return &Compartments;
}

 

CSM::~CSM()
{
	Compartments.~vector();
	Months.~vector();
	printf ("Characters: %c %c \n", 'a', 65);
}

void CSM::InitCompartments()
{
	Compartment *newCompartment;
	for (int i = 0; i < nNumOfCompatments; i++)
	{
		newCompartment = new Compartment(i, nArea, wieltingPoint, nFieldCapacity, thick, CCa, CSO4);
		Compartments.push_back(*newCompartment);
		nTotalWP += wieltingPoint * thick;
		nTotalWhc += nFieldCapacity * thick;
		nTotalMoist += newCompartment->nMoist;
		nInitMoistTotal = nTotalMoist;
	}
}

void CSM::InitMonths()
{
  //NumericVector returnList = NumericVector::create(1,2,3,4,5,6,7,8,9,10,11,12);
  float II = 0.0;
	float ca = 0.0;
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
	MONTH month;
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


//RCPP_MODULE(unif_module) {
//  class_<CSM>( "CSM" )
//  .constructor<int>()
//  .field( "min", &CSM::nNumOfDays )
//  ;
//}