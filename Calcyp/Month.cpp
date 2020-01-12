#include "Month.h"






Month::Month(int _Number, double _temp, double _PanDaily, double _PETdaily, double _PET, double _pan) 
	: nNumber(MONTH(_Number)) , nTemp(_temp), PanDaily(_PanDaily), PETdaily(_PanDaily), PET(_PET), pan(_pan), totalAET(0)
{

}

Month::~Month()
{
}
