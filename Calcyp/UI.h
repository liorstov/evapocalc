#pragma once
#include "CSM.h"
#include "Rain.h"
class UI
{
public:
	UI(CSM& _handler, WG& _rain);
	~UI();

	void Print();
	CSM &Handler;
	WG &rainHandler;
};

