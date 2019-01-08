#pragma once
#include <vector>

using namespace std;

class Rain
{
public:
	Rain(bool bCreateRain);

	~Rain();

	float* rain;
	float totalrain;

private:
	vector<pair<float, float> > Interval, Amount;

	int getRandomInterval();
	float getRandomAmount();
};

