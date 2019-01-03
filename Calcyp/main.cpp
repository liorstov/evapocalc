#pragma once
#include "CSM.h"
#include "Rain.h"

int main()
{
	int Years = 100;
	float temp[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
	Rain Rainclass(false, Years);
	


	CSM Calc(Years, 5, 20, *Rainclass.rain,*temp, 0.039, 0.122, 0.0000113, 0.51 / (365.0*10000.0));
	Calc.Calculate();

	//UI ui(Calc, Rainclass);
	//ui.Print();

	
	


}