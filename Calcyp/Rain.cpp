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

	string temp, year;
	// Read the Data from the file 
	// as String Vector 
	year = "0";
	if (fin.is_open()) {
		getline(fin, temp);
		while (fin.good() && stoi(year) <= 2000) {
			getline(fin, temp, ',');
			getline(fin, year, ',');
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

void WG::writeRes(vector<double> res)
{
	ofstream fout("../Res.csv", ios::out);
	if (fout.is_open())
	{
		for (size_t i = 0; i < res.size(); i++)
		{
			fout << res[i] << endl;
		}
		fout.close();
	}
}

int WG::getRandomInterval()
{
	return 0;
}

double WG::getRandomAmount()
{
	return 0.0f;
}
