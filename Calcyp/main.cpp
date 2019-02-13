#include "CSM.h"
#include "Rain.h"
#include "UI.h"

int main()
{

	float temp[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
	Rain Rainclass(false);
	

	
	CSM * Calc = new CSM();
	
	Calc->Calculate(Rainclass.rain, 1, 50, 5, 0.039, 10, 10, 1.44, 1, 0.2, 0, 0, );
	UI ui(*Calc, Rainclass);
	ui.Print();

	//system("pause");
	

}

