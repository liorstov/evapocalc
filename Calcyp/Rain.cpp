#include "Rain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

#define RAND rand() % 9999 + 1



Rain::Rain(bool bCreateRain)
{
	std::string tmp;
	vector<string> words;
	char delim = ','; // Ddefine the delimiter to split by
	
	
	rain = new float[1000*365];
	ifstream file("./../DB/rainseriestemp.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	string value;
	int column = 1;
	int i = 0;
	while (std::getline(file, tmp)) {
		rain[i++] = stof(tmp);
	}

	
	file.close();
	//srand(time(NULL));
	//int _NumberOfDays = _NumberOfYears * 365;
	//rain = new float[_NumberOfDays + 1];
	//fstream RainArr;
	//RainArr.open("./../DB/rain.txt", fstream::out | fstream::in | fstream::trunc);

	////rain data already exist
	//if (!true)
	//{
	//	for (size_t i = 0; i < _NumberOfDays; i++)
	//	{
	//		RainArr >> rain[i];
	//	}
	//	RainArr.close();
	//	return;
	//}

	//int i = 0;


	//ifstream QuantTucson, IntervalTucson;
	//QuantTucson.open("./../DB/TucsonAmount.txt");
	//IntervalTucson.open("./../DB/TucsonInterval.txt");

	//if (!QuantTucson.is_open() || !IntervalTucson.is_open())
	//	return;


	//float Quantity[75];
	//float Interval[75];
	//for (size_t i = 0; i < 75; i++)
	//{
	//	QuantTucson >> Quantity[i];

	//}
	//for (size_t i = 0; i < 58; i++)
	//{
	//	IntervalTucson >> Interval[i];
	//}

	//float *sizerain = new float[500000];
	//float *sizeevent = new float[500000];

	//float *rainevent = new float[500000];

	//for (int day = 1; day <= _NumberOfDays; day++) rain[day] = 0.0;

	//// loop that determines when a rain event will occur
	//int n = 1;
	//for (i = 10000; i >= 1; i--)
	//{
	//	if (Interval[n] > 1.0*i / 10000) n++;
	//	sizeevent[i] = n;
	//}

	//// loop that determines how much rainfall will be in each rain day
	//float j = 0.05;
	//for (i = 10000; i >= 1; i--)
	//{
	//	if (Quantity[(int)(j*20.0 + 1.0)] > 1.0*i / 5000) j = j + 0.05;
	//	size =[i] = (j - 0.05) + (0.05*((float)rand() / RAND_MAX));
	//} // for continuous values
	//  //size[i]=j;} // for discrete values

	//int interval = sizeevent[RAND];
	//int day = 0;
	//int q = 1;
	//totalrain = 0;

	//// loop for counting the days and detemine when will be a rain day and the amount of rainfall
	//for (size_t day = 0; day < _NumberOfDays; day++)
	//{
	//	if (day == 0)
	//	{
	//		rain[day] = sizerain[RAND];
	//		rainevent[q] = day;
	//	}

	//	if (day == rainevent[q] + interval)
	//	{
	//		rain[day] = sizerain[RAND];
	//		interval = sizeevent[RAND];
	//		rainevent[q + 1] = day;
	//		q++;
	//	}
	//	totalrain = totalrain + rain[day];
	//	/*if (day == 364999)
	//	{
	//		int s = 2;
	//	}*/
	//	RainArr << rain[day] << std::endl;
	//}
	//RainArr << totalrain << std::endl;

	//delete[] sizerain;
	//delete[] sizeevent;
	//delete[] rainevent;

	//QuantTucson.close();
	//IntervalTucson.close();

}


Rain::~Rain()
{
	 
	delete[] rain;
}

int Rain::getRandomInterval()
{
	return 0;
}

float Rain::getRandomAmount()
{
	return 0.0f;
}
