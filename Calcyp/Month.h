#pragma once
typedef enum { Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec } MONTH;

class Month
{
public:
	Month(int _Number, float _temp, float _PanDaily, float _PETdaily, float _PET, float _pan);
	~Month();

	MONTH nNumber;
	float nTemp;
	float PanDaily;
	float PETdaily;
	float PET;
	float pan;
};

