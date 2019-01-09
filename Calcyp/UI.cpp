#include "UI.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

UI::UI(CSM& _handler, Rain& _rain) : Handler(_handler), rainHandler(_rain)
{
	//Handler = _handler;
}

void UI::Print()
{
	ofstream fOutput;
	fOutput.open("./../DB/output.txt", fstream::out);
	const std::vector<Compartment> &Comp =  *Handler.GetCompartments();

	float totalcaco3 = 0;
	if (fOutput.is_open())
	{
		fOutput << "time (years) = " << Handler.nNumOfDays << endl;
		fOutput << "depth(cm)				caso4(g/cm3) 				 Moist[cm]\n";
		for (int d = 0;d < Handler.nNumOfCompatments;d++)
		{
			fOutput <<  (d)*Handler.thick <<  "\t\t\t\t" << Comp[d].C_CaSO4 << "\t\t\t\t" << Comp[d].nMoist<< endl;
			
		}

		fOutput << "\nCa mass balance (g/cm2):\n";
		fOutput << "Ca accumulation from dust		Ca accumulation from rain		Ca in soil as CaCO3		leachate Ca\n";
		fOutput << 0.4*Handler.nTotalCaDust << "\t\t\t\t\t\t" << Handler.nTotalRain << "\t\t\t\t\t\t" << totalcaco3*0.4 << "\t\t\t\t\t\t" << Handler.nTotalCaLeachate*0.4 << endl;
		fOutput << "balance = ";
		fOutput << (0.4*Handler.nTotalCaDust) + Handler.nTotalRain - (totalcaco3*0.4) - (Handler.nTotalCaLeachate*0.4) << endl << endl;
		fOutput << "water mass balance (cm):\n";
		fOutput <<  "intial soil water		rainfall		final soil water		AET			leachate\n";
		fOutput << Handler.nInitMoistTotal << "\t\t\t\t" << Handler.nTotalRain << "\t\t\t\t" << Handler.nTotalMoist << "\t\t\t\t" << Handler.nTotalAet << "\t\t\t\t" << Handler.nTotalLeachate << endl;
		fOutput << "balance = " <<  Handler.nTotalRain - Handler.nTotalMoist - Handler.nTotalAet - Handler.nTotalLeachate << endl;
	}


	//	fprintf(fp8, "time (years) = ");
	//fprintf(fp8, "%d\n\n", numberofyears);
	//fprintf(fp8, "depth(cm)       pco2(atm)              caco3(g/cm3)		caco3(g/cm2)		caco3(g/m2)\n");
	//for (d = 1;d <= depth;d++)
	//{
	//	fprintf(fp8, "%0.0f-%0.0f\t\t%0.6f\t\t%0.6f\t\t\%0.6f\t\t\%0.1f\n", (d - 1)*cthick, d*cthick, pco2[d], caco3[d], caco3[d] * cthick, caco3[d] * cthick * 10000);
	//	totalcaco3 = totalcaco3 + caco3[d] * cthick;
	//}
	////fprintf(fp8,"\n%0.3f\n",totalcaco3);
	//fprintf(fp8, "\nCa mass balance (g/cm2):\n");
	//fprintf(fp8, "Ca accumulation from dust	Ca accumulation from rain	Ca in soil as CaCO3	leachate Ca\n");
	//fprintf(fp8, "\t%0.5f\t\t\t\t%0.5f\t\t\t\t%0.5f\t\t%0.5f\n", 0.4*totalcadust, totalcarain, totalcaco3*0.4, totalcaleach*0.4);
	//fprintf(fp8, "balance = ");
	//fprintf(fp8, "%0.5f\n\n", (0.4*totalcadust) + totalcarain - (totalcaco3*0.4) - (totalcaleach*0.4));
	//fprintf(fp8, "water mass balance (cm):\n");
	//fprintf(fp8, "intial soil water	rainfall	final soil water	AET		leachate\n");
	//fprintf(fp8, "\t%0.2f\t\t%0.2f\t\t%0.2f\t\t\t%0.2f\t\t%0.2f\n", initialmoisttotal, totalrain, moisttotal, totalAET, totalleachate);
	//fprintf(fp8, "balance = ");
	//fprintf(fp8, "%0.3f\n", initialmoisttotal + totalrain - moisttotal - totalAET - totalleachate);

	fOutput.close();
}

UI::~UI()
{

}
