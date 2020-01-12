#pragma once
#include <vector>

using namespace std;

class WG
{
public:
	WG();

	~WG();

	vector<double> rain;
	vector<double> PET;
	double totalrain;

	void writeRes(vector<double>);
private:
	vector<pair<double, double> > Interval, Amount;

	int getRandomInterval();
	double getRandomAmount();
};

