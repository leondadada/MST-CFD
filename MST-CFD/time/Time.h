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
#include <iomanip>


class Time
{
public:
	Time(MshBlock*, fstream*,AllData*);
	~Time();
	void goNextTimeStep();
	VCTDIMU* getP1NewCells();
	VCTDIMU* getP1OldCells();
	VCTDIMU * getP1OldPUVTs();
	void initialization(VCTDIMU);
private:
	fstream* pFlog;
	NUM DT;
	MshBlock* pMesh;
	AllData* pAllData;
	vector<Cell>::iterator itBeginCell;
	vector<Face>::iterator itBeginFace;
	vector<FacesInf>::iterator itBeginFacesInf;
		
	
};
