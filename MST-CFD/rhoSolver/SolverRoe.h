#pragma once
#include "../include/CONST.h"
#include "../mesh/Cell.h"
#include "../work/FUNCTION.h"
#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;

class solverRoe
{
public:
	void set(VCTDIMU, VCTDIMU);

	VCTDIMU solverx();
	VCTDIMU solvery();
	VCTDIMU solverAll(int);
private:

	VCTDIMU rightQ ;
	VCTDIMU leftQ ;

	NUM roe_rho;
	VCTDIM roeU;
	NUM roe_ht;
	NUM roe_a;
	NUM D;
	NUM entropyRefine(NUM );
	
	
};

