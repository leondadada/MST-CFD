#include "SolverAusm.h"

void SolverAusm::set(VCTDIMU _LeftQ, VCTDIMU _RightQ) {
	leftQ = _LeftQ;
	rightQ = _RightQ;
	NUM aLStar = sqrt(2 * getht(leftQ) * (GAMMA - 1) / (GAMMA + 1));
	NUM aRStar = sqrt(2 * getht(rightQ) * (GAMMA - 1) / (GAMMA + 1));
#if DIM ==2
	
	NUM UL = sqrt((leftQ(1)* leftQ(1) + leftQ(2)* leftQ(2)) / leftQ(0) / leftQ(0));
	NUM UR = sqrt((rightQ(1) * rightQ(1) + rightQ(2) * rightQ(2)) / rightQ(0) / rightQ(0));
	
#else
	
	NUM UL = sqrt((leftQ(1)* leftQ(1) + leftQ(2)* leftQ(2) + leftQ(3)* leftQ(3)) / leftQ(0) / leftQ(0));
	NUM UR = sqrt((rightQ(1) * rightQ(1) + rightQ(2) * rightQ(2) + rightQ(3) * rightQ(3)) / rightQ(0) / rightQ(0));
	
#endif

	aL = aLStar * aLStar / std::max(aLStar, UL);
	aR = aRStar * aRStar / std::max(aRStar, UR);
	this->aFace = std::min(aL, aR);
	//may need recalculate aL\aR??
	this->pL = getP(leftQ);
	this->pR = getP(rightQ);
}

VCTDIMU  SolverAusm::solverx() {
	//calculate Ausm flux on face on X
	//calculate aFace simple :
	this->machL = leftQ(1) / leftQ(0) / aFace;
	this->machR = rightQ(1) / rightQ(0) / aFace;

	calculateMachAndPFace_AusmPlus();
	F_ccFace = 0.5 * (machFace * (aL*F_caL + aR*F_caR) - abs(machFace) * (aR * F_caR - aL*F_caL));

	pSVRDIMU = this->F_ccFace;
	pSVRDIMU(1) += this->pFace;
	return this->pSVRDIMU;
}
VCTDIMU  SolverAusm::solvery() {
	//calculate roe flux on face on Y
	this->machL = leftQ(2) / leftQ(0) / aFace;
	this->machR = rightQ(2) / rightQ(0) / aFace;

	calculateMachAndPFace_AusmPlus();
	F_ccFace = 0.5 * (machFace * (aL * F_caL + aR * F_caR) - abs(machFace) * (aR * F_caR - aL * F_caL));

	pSVRDIMU = this->F_ccFace;
	pSVRDIMU(2) += this->pFace;
	return this->pSVRDIMU;
}
VCTDIMU  SolverAusm::solverAll(int _iDim) {
	this->machL = leftQ(_iDim + 1) / leftQ(0) / aFace;
	this->machR = rightQ(_iDim + 1) / rightQ(0) / aFace;

	calculateMachAndPFace_AusmPlus();
	F_ccFace = 0.5 * (machFace * (aL * F_caL + aR * F_caR) - abs(machFace) * (aR * F_caR - aL * F_caL));

	pSVRDIMU = this->F_ccFace;
	pSVRDIMU(_iDim + 1) += this->pFace;
	return this->pSVRDIMU;
}

void SolverAusm::calculateMachAndPFace_Ausm() {

	F_caL = leftQ;
	F_caL(DIMU - 1) += getP(leftQ);
	F_caR = rightQ;
	F_caR(DIMU - 1) += getP(rightQ);

	this->aL = sqrt(GAMMA * getP(leftQ) / leftQ(0));
	this->aR = sqrt(GAMMA * getP(rightQ) / rightQ(0));

	NUM machPlus, machMinus;//mach number:
	if (abs(machL <= 1)) {
		machPlus = 0.25*(machL + 1)*(machL + 1);
	}
	else {
		machPlus = 0.5*(machL + abs(machL));
	}
	if (abs(machR <= 1)) {
		machMinus = -0.25*(machR - 1)*(machR - 1);
	}
	else {
		machMinus = 0.5*(machR - abs(machR));
	}
	NUM pPlus, pMinus;//P Ausm:
	if (abs(machL <= 1)) {
		pPlus = pL*0.25*(machL + 1)*(machL + 1)*(2 - machL);
	}
	else {
		pPlus = pL*0.5*(machL + abs(machL)) / machL;
	}
	if (abs(machR <= 1)) {
		pMinus = pR*0.25*(machR - 1)*(machR - 1)*(2 + machR);
	}
	else {
		pMinus = pR*0.5*(machR - abs(machR)) / machR;
	}

	this->machFace = machMinus + machPlus;
	this->pFace = pMinus + pPlus;
}
void SolverAusm::calculateMachAndPFace_AusmPlus() {
	//AUSM+ ,
	F_caL = leftQ;
	F_caL(DIMU - 1) += getP(leftQ);
	F_caR = rightQ;
	F_caR(DIMU - 1) += getP(rightQ);

	this->aL = sqrt(GAMMA * getP(leftQ) / leftQ(0));
	this->aR = sqrt(GAMMA * getP(rightQ) / rightQ(0));

	NUM machPlus, machMinus;//mach number:
	if (abs(machL <= 1)) {
		machPlus = 0.25 * (machL + 1) * (machL + 1) + 0.125 * (machL * machL - 1) * (machL * machL - 1);
	}
	else {
		machPlus = 0.5*(machL + abs(machL));
	}
	if (abs(machR <= 1)) {
		machMinus = -0.25 * (machR - 1) * (machR - 1) - 0.125 * (machR * machR - 1) * (machR * machR - 1);
	}
	else {
		machMinus = 0.5*(machR - abs(machR));
	}
	NUM pPlus, pMinus;//P Ausm:
	if (abs(machL <= 1)) {
		pPlus = pL * 0.25 * (machL + 1) * (machL + 1) * (2 - machL) + 0.1875 * machL * (machL * machL - 1) * (machL * machL - 1);
	}
	else {
		pPlus = pL * 0.5 * (machL + abs(machL)) / machL;
	}
	if (abs(machR <= 1)) {
		pMinus = pR * 0.25 * (machR - 1) * (machR - 1) * (2 + machR) - 0.1875 * machR * (machR * machR - 1) * (machR * machR - 1);
	}
	else {
		pMinus = pR * 0.5 * (machR - abs(machR)) / machR;
	}
	this->machFace = machMinus + machPlus;
	this->pFace = pMinus + pPlus;
}

