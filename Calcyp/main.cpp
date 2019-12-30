#include "CSM.h"
#include "Rain.h"
#include "UI.h"

int main()
{

	float temp[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
	WG rainPET;
	
	
	
	CSM * Calc = new CSM();
	
	vector<float> res = Calc->Calculate(rainPET.rain, rainPET.PET, 900, 100, 5, 0.02, 0, 0, 1.44, 1, 0.19, 1.5, 1.5,0.6 );

	delete(Calc);
	system("pause");
	

}

