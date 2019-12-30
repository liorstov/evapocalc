#pragma once
#include <vector>

using namespace std;

class WG
{
public:
	WG();

	~WG();

	vector<float> rain;
	vector<float> PET;
	float totalrain;

private:
	vector<pair<float, float> > Interval, Amount;

	int getRandomInterval();
	float getRandomAmount();
};

