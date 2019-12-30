#include "Rain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

#define RAND rand() % 9999 + 1



WG::WG()
{

	// Open an existing file 
	ifstream fin("../synthRain.csv", ios::in);

	float I, year, day, rain, ser, PET, K;
	string temp;
	// Read the Data from the file 
	// as String Vector 

	if (fin.is_open()) {
		getline(fin, temp);
		while (fin.good()) {
			getline(fin, temp, ',');
			getline(fin, temp, ',');
			getline(fin, temp, ',');
			getline(fin, temp, ',');

			this->rain.push_back(atof(temp.c_str()));
			getline(fin, temp, ',');
			getline(fin, temp, ',');

			this->PET.push_back(atof(temp.c_str()));
			getline(fin, temp);

		}
	}
		

}


WG::~WG()
{
	 
	
}

int WG::getRandomInterval()
{
	return 0;
}

float WG::getRandomAmount()
{
	return 0.0f;
}
