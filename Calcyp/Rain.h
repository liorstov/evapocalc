#pragma once
#include <vector>

using namespace std;

class Rain
{
public:
	Rain(bool bCreateRain, int _NumberOfDays);

	~Rain();

	float* rain;
	float totalrain;

private:
	vector<pair<float, float> > Interval, Amount;

	int getRandomInterval();
	float getRandomAmount();
};

