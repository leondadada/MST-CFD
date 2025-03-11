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

class RhoSolver
{
public:
	RhoSolver(MshBlock* , fstream* , AllData* );
	~RhoSolver();
	void setDT(NUM);
	void solve();
	VCTDIMU* getOldValue();
	VCTDIMU* getOldNTimeValue();
	void updateNewToOld();
	VCTDIMU* getNewValue();

private: 
	fstream* pFlog;
	NUM DT;
	MshBlock* pMesh;
	AllData* pAllData;
	Cell* itBeginCell;
	Face* itBeginFace;
	vector<FacesInf>::iterator itBeginFacesInf;


	VCTDIMU* p1NewCellQs;
	VCTDIMU* p1OldFaceQs;
	VCTDIMU* p1OldCellQs;


	VCTDIMU* p1OldFaceConvectFluxF;
	VCTDIMU* p1OldFaceConvectFluxG;
	MTDIMU_DIM* pCellGradPUVT;
	MTDIMU_DIM* pFaceGradPUVT;
	MTDIMU_DIM* p1OldFaceConvectFlux;
	MTDIMU_DIM* p1OldFaceViscidFlux;
	MTDIMU_DIM* p1NewCellGradFlux;
	MTDIMU_DIM* p1NewFaceGradFlux;

	VCTDIMU* ptP1OldNTimeCellQs;
	void updateFaceFlux();
	void updateGradFlux();
	void updateFaceFlux2ndOrder();
	void updateViscid();

};

