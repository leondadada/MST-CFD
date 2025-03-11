#pragma once
#include "../include/CONST.h"
#include "../mesh/Cell.h"
#include "../work/FUNCTION.h"

#include <Eigen/Dense>
#include "math.h"
#include <algorithm>
class SolverAusm
{
public:
	void set(VCTDIMU, VCTDIMU);

	VCTDIMU solverx();
	VCTDIMU solvery();

	VCTDIMU solverAll(int _iDim);

private:

	VCTDIMU rightQ;
	VCTDIMU leftQ;
	NUM aL;
	NUM aR;
	NUM aFace;//acoustic speed
	NUM machL;
	NUM machR;
	NUM machFace;
	NUM pL;
	NUM pR;
	NUM pFace;
	void calculateMachAndPFace_Ausm();
	void calculateMachAndPFace_AusmPlus();
	VCTDIMU F_caL;
	VCTDIMU F_caR;

	VCTDIMU pSVRDIMU;
	VCTDIMU F_ccFace;//point, new space in rom to use outside  
	
	

};



