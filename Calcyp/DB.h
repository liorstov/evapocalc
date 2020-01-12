#pragma once
#include <iostream>
class DB
{
public:
	
	~DB();

private:

	DB();
	static DB instance;

	double RainAmount; //input is probability of getting rainfall greater than x on any day, where x is in mm / d from 0.05 cm / d to 3.05 cm / d
	double	RainFall; //output the rinfall per day
	double	AvrRain; // average annual rainfall
	double	RainEvent; // input is probability of getting rainfall event greater than x, where x is in days from 1 day to 73 days
	double AvrMonthTemp; // input is average monthly temprature in celcius
	double RainAmountProb; //holds the probability for rain amounts
	double	RainIntervalProb; //holds the probability for interval between rain events
	double	Config;

};

