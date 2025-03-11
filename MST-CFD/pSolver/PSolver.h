#pragma once
#include "../mesh/MshBlock.h"
#include "../rhoSolver/SolverAusm.h"
#include "../rhoSolver/SolverRoe.h"
#include "../rhoSolver/SolverNND.h"
#include "../work/FUNCTION.h"
#include "../lusolver/sparseSolverNUM.h"
#include "../lusolver/sparseSolver.h"
#include "../lusolver/sparseSolver.cpp"
#include "../rhoSolver/RhoSolver.h"
#include "../pSolver/PSolver.h"
#include "../data/AllData.h"
#include "../pSolver/SolverSIMPLE.h"
#include "../pSolver/SolverSIMPLE2.h"
#include "../pSolver/SolverCoupled.h"

class PSolver
{
public:
	PSolver(MshBlock*, fstream*, AllData*);
	~PSolver();
	void setDT(NUM);
	void solve();
	VCTDIMU * getOldValue();
	VCTDIMU * getNewValue();
	VCTDIMU * getOldNTimeValue();
	void updateNewToOld();
private:
	fstream* pFlog;
	NUM DT;
	MshBlock* pMesh;
	AllData* pAllData;

	Cell* itBeginCell;
	Face* itBeginFace;
	vector<FacesInf>::iterator itBeginFacesInf;



};

