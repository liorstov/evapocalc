#include "CSM.h"
#include "Rain.h"
#include "UI.h"

int main()
{

	double temp[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
	Rain Rainclass(false);
	

	
	CSM * Calc = new CSM();
	
	Calc->Calculate(Rainclass.rain, 1000, 100, 5, float(0.039), 33, 33, float(1.44), 1, float(0.2), 1000, 1000);
	UI ui(*Calc, Rainclass);
	ui.Print();

	//system("pause");
	

}

