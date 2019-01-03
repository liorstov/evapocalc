#include "CSM.h"
#include <ctime>




//This code calculate the content and depth of pedogenic carbonate with depth
CSM::CSM(float rain)
{
	nNumOfDays = 10000 * 365;
	nDepth = 30;
	thick = 5;
	wieltingPoint = 0.039;
	CCa = 0.122;
	nArea = 1;
	nFieldCapacity = 0.0000113;
	nDust = 0.51 / (365.0*10000.0);
	nCaco3less = 0;
	nTotalWhc = 0;
	nTotalCaDust = 0;
	nTotalCaRain = 0;
	nTotalCaLeachate = 0;
	nTotalLeachate = 0;
	nTotalMoist = 0;
	nTotalAet = 0;
	nTotalWP = 0;
	nTemp = 0;
	TempArr = new float[12];
	//TempArr[] = {10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8};
	nLeachate = 0;
	RainArr = &rain;

	AET = new float[nNumOfDays];

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
		Compartments[0].nCO2 = Compartments[0].nICO2 * CO2Mul[nMonth];
		 
		nTotalCaDust += nDust;
		nTotalCaRain += RainArr[day] * CCa*40.0 / 1000.0;
		if (nTotalMoist >= ((0.546*nTotalWhc) + nTotalWP)) // according to Marion et al. (1985), for the upper 45% of the total whc the actual evapotranspiration (AET) is the potential evapotranspiration (pet)
			AET[day] = Months[(int)GetMonth(day)].PanDaily;		// in case of 10 compartments of 10 cm each, if total moistute > 8.465 
														//AET[day]=PETdaily[monthperday[day]];
		else										// the lower 55% of the total whc are according to modifeid Thornthwaite-Mather model
			AET[day] = (nTotalMoist - nTotalWP) / (0.546*nTotalWhc)*Months[(int)GetMonth(day)].PanDaily;
			//AET[day]=((moisttotal*0.352)-1.36)*pandaily[monthperday[day]];

			nTotalAet += AET[day];

			
		Compartments[0].nMoist = Compartments[0].nInitMoist + RainArr[day];   // set the moiture of the 1st compartment to the intial moisture plus the daily rainfall. rainfall is added only the 1st compartment

		if ((Compartments[0].nMoist - (Compartments[0].nThetaWeildingPnt * thick))<AET[day]) // examine if the daily AET is greater than the available water for evaporation at the 1st compartment
		{
			AET[day] = AET[day] - (Compartments[0].nMoist - (Compartments[0].nThetaWeildingPnt * thick));           // if yes, evaporate all the available moisture and set the moisture to the wilting point
			Compartments[0].nMoist = Compartments[0].nThetaWeildingPnt * thick;
		}
		else
		{
			Compartments[0].nMoist -= AET[day];             // if no, substract the AET from the moisture of the 1st compartment
			AET[day] = 0.0;
		}
		// in both cases, the AET is updated: if the AET is greater than the available water for evaporation at the 1st compartment, than
		// substract the part that was evaporated from the AET. Other wise, set AET to zero. 

		Compartments[0].nCaCO3 += (nDust / thick); // add the daily dust accumulation to the 1st compartment
		

		if (((Compartments[0].nThetaWeildingPnt + Compartments[0].nWhc)*thick)<Compartments[0].nMoist) // examine if the updated moisture of the 1st compartment (after rainfall addition and AET substraction)
													 // is greater than field capacity (=saturation). if yes, starts the solubility calculation
		{
			CCa = Compartments[0].solubility(nTemp); // calculates the carbonate solubility for the 1st compartment. ca equals to the amount of caco3 [g/cm3]
			if (Compartments[0].nCaCO3>CCa) nCaco3less = Compartments[0].nCaCO3 - CCa; // if the caco3 in the 1st compartment is greater than the soluble caco3 (=ca) then caco3less is set to the difference between them
			else nCaco3less = 0.0;                         // if not, caco3less is set to zero. note that caco3less is the dissolution/percipitation amount of carbonate
			Compartments[0].nCaCO3 = Compartments[0].nCaCO3 - nCaco3less;
			Compartments[1].nCaCO3 += nCaco3less;
		}

		// This is the second loop that runs through the soil profile, from the 2nd compartment to the bottom
		for (int d = 1;d < nDepth;d++)
		{
			// set the pco2 values to all compartments
			Compartments[d].nCO2 = Compartments[d].nICO2 * CO2Mul[nMonth];

			if (Compartments[d - 1].nMoist>(Compartments[d - 1].nWhc + Compartments[d - 1].nThetaWeildingPnt)*thick)
				nLeachate = Compartments[d - 1].nMoist - (Compartments[d - 1].nWhc + Compartments[d - 1].nThetaWeildingPnt)*thick; // determines the leachate from the upper compartment by substracting the field capacity from the moisture content
			else nLeachate = 0.0; // if there is no excess water, then the leachate is zero
			Compartments[d - 1].nMoist -= nLeachate; // subtract the leachate from the moisture of the upper compartment
			Compartments[d].nMoist = Compartments[d].nInitMoist + nLeachate;    // add the leachate to the moisture of the current compartment

			if (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt * thick)<AET[day]) // same as for the 1st compartment, taking into account the AET for this current compartment, and updating the AET value
			{
				AET[day] = AET[day] - (Compartments[d].nMoist - (Compartments[d].nThetaWeildingPnt * thick));
				Compartments[d].nMoist = Compartments[d].nThetaWeildingPnt * thick;
			}
			else
			{
				Compartments[d].nMoist -= AET[day];
				AET[day] = 0.0;
			}

			if (((Compartments[d].nThetaWeildingPnt + Compartments[d].nWhc)*thick)<Compartments[d].nMoist) // same as for the 1st compartment, examine if we reached saturation
			{
				CCa = Compartments[d].solubility(nTemp);
				if (Compartments[d].nCaCO3>CCa)
					nCaco3less = Compartments[d].nCaCO3 - CCa;
				else nCaco3less = 0.0;
				Compartments[d].nCaCO3 -= nCaco3less;
				if (d != nDepth-1) Compartments[d + 1].nCaCO3 += nCaco3less;
				else nTotalCaLeachate += nCaco3less;
			}

			if (d == nDepth-1) // examines if we reached the compartment above the lowest one 
			{
				if (AET[day]>0) nTotalAet -= AET[day];

				//   if (moist[d+1]>((whc[d+1]+thetawp[d+1])*cthick)) 
				//   {moist[d+1]=moist[d+1]-(moist[d+1]-(whc[d+1]+thetawp[d+1])*cthick); // calculate the moisture for the last compartment. if it is above field capacity, substracts the excess water from the moisture of the lowest compartment (i.e. the water is leached to the rock)
				//    totalleachate=totalleachate+(moist[d+1]-(whc[d+1]+thetawp[d+1])*cthick);}
				//   else
				//{moist[d+1]=moist[d+1]-0.0;}

				if (Compartments[d].nMoist>(Compartments[d].nWhc + Compartments[d].nThetaWeildingPnt)*thick)
					nLeachate = Compartments[d].nMoist - (Compartments[d].nWhc + Compartments[d].nThetaWeildingPnt)*thick; // determines the leachate from the upper compartment by substracting the field capacity from the moisture content
				else nLeachate = 0.0; // if there is no excess water, then the leachate is zero
				Compartments[d - 1].nMoist -= nLeachate; // subtract the leachate from the moisture of the upper compartment
				nTotalLeachate += nLeachate;


				nTotalMoist = 0.0;
				//pco2[d+1]=PCO2(d+1,julian,cthick);
				for (int w = 1;w < nDepth;w++) // resert the initial soil moisture as the current soil moisture, and calculate the total moisture
				{
					Compartments[w].nInitMoist = Compartments[w].nMoist;
					nTotalMoist += Compartments[w].nInitMoist;
					//fprintf(fp5,"%d\t%0.3f\n",w,moist[w]);
				}

				//fprintf(fp5,"%ld\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",day,rain[day],pandaily[monthperday[day]],moisttotal,AET[day]);

			}

		}

	}
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
	delete[] AET;
}

void CSM::InitCompartments()
{
	Compartment *newCompartment;
	for (int i = 0; i < nDepth; i++)
	{
		newCompartment = new Compartment(i, nArea, wieltingPoint, nFieldCapacity, thick, CCa, CO2[i]);
		Compartments.push_back(*newCompartment);
		nTotalWP += wieltingPoint * thick;
		nTotalWhc += nFieldCapacity * thick;
		nTotalMoist += newCompartment->nInitMoist;
		nInitMoistTotal = nTotalMoist;
	}
}

NumericVector CSM::InitMonths()
{
  NumericVector returnList = NumericVector::create(1,2,3,4,5,6,7,8,9,10,11,12);
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
		returnList[i] = 20;
	}
	
	return returnList;

	


}

MONTH CSM::GetMonth(int nDay)
{
	MONTH month;
	nDay = JULIAN(nDay);
	if (nDay <= 31)
		month = Jan;
	else
		if (nDay <= 59)
			month = Feb;
		else
			if (nDay <= 90)
				month = Mar;
			else
				if (nDay <= 120)
					month = Apr;
				else
					if (nDay <= 151)
						month = May;
					else
						if (nDay <= 181)
							month = Jun;
						else
							if (nDay <= 212)
								month = Jul;
							else
								if (nDay <= 243)
									month = Aug;
								else
									if (nDay <= 273)
										month = Sep;
									else
										if (nDay <= 304)
											month = Oct;
										else
											if (nDay <= 334)
												month = Nov;
											else
												if (nDay <= 365)
													month = Dec;
	return month;

}

float CSM::JULIAN(int day)

{
	int julian;

	julian = day - ((long)(day / 365) * 365);
	return julian;
}


RCPP_MODULE(unif_module) {
  class_<CSM>( "CSM" )
  .constructor<int>()
  .field( "min", &CSM::nNumOfDays )
  ;
}