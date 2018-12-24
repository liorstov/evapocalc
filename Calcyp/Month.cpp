#include "Month.h"






Month::Month(int _Number, float _temp, float _PanDaily, float _PETdaily, float _PET, float _pan) 
	: nNumber(MONTH(_Number)) , nTemp(_temp), PanDaily(_PanDaily), PETdaily(_PanDaily), PET(_PET), pan(_pan)
{

}

Month::~Month()
{
}
