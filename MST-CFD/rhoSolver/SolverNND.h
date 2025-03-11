#pragma once
#include "../include/CONST.h"
#include "../mesh/Cell.h"


class SolverNnd
{
public:
	void set(Cell*, Cell*);

	VCTDIMU solverx();
	VCTDIMU solvery();

private:

	Cell leftCell;
	Cell rightCell;
	VCTDIMU pMatrix;
};

