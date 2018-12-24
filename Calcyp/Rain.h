#pragma once
class Rain
{
public:
	Rain(bool bCreateRain, int _NumberOfDays);
	~Rain();

	float* rain;
	float totalrain;
};

