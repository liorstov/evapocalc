#pragma once
#include "CSM.h"
#include "Rain.h"
class UI
{
public:
	UI(CSM& _handler, Rain& _rain);
	~UI();

	void Print();
	CSM &Handler;
	Rain &rainHandler;
};

