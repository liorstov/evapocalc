#pragma once
#include <iostream>
class DB
{
public:
	
	~DB();

private:

	DB();
	static DB instance;

	float RainAmount; //input is probability of getting rainfall greater than x on any day, where x is in mm / d from 0.05 cm / d to 3.05 cm / d
	float	RainFall; //output the rinfall per day
	float	AvrRain; // average annual rainfall
	float	RainEvent; // input is probability of getting rainfall event greater than x, where x is in days from 1 day to 73 days
	float AvrMonthTemp; // input is average monthly temprature in celcius
	float RainAmountProb; //holds the probability for rain amounts
	float	RainIntervalProb; //holds the probability for interval between rain events
	float	Config;

};

