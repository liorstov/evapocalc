#include "CSM.h"
#include "Rain.h"
#include "UI.h"
int main()
{

	double temp[12] = { 10.2,11.8,14.3,18.1,22.6,27.8,30.0,28.9,26.7,20.8,14.6,10.8 };
	WG rainPET;
	
	
	
	CSM * Calc = new CSM();
	
	vector<double> res = Calc->Calculate(rainPET.rain, rainPET.PET, 1000, 80, 5, 0.013, 1, 0.16, 128, 2.4, 0.6, 0, 0.0152 / 365,35.48,20);
	//system("cls");
	rainPET.writeRes(res);

	delete(Calc);
	system("pause");
	

}

