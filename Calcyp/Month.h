#pragma once
typedef enum { Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec } MONTH;

class Month
{
public:
	Month(int _Number, double _temp, double _PanDaily, double _PETdaily, double _PET, double _pan);
	~Month();

	MONTH nNumber;
	double nTemp;
	double PanDaily;
	double PETdaily;
	double PET;
	double pan;	
	double totalAET;

};

