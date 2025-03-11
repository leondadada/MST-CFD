#pragma once
#include "../include/CONST.h"
#include "../mesh/Cell.h"
#include "../work/FUNCTION.h"
#include "../mesh/MshBlock.h"
#include "../rhoSolver/SolverAusm.h"
#include "../rhoSolver/SolverRoe.h"
#include "../rhoSolver/SolverNND.h"
#include "../lusolver/sparseSolverNUM.h"
#include "../lusolver/sparseSolver.h"
#include "../lusolver/sparseSolver.cpp"
#include "../rhoSolver/RhoSolver.h"
#include "../pSolver/PSolver.h"
#include "../data/AllData.h"
#include <Eigen/Dense>
#include "math.h"
#include <algorithm>
using namespace std;

class SolverSIMPLE {
public:
	SolverSIMPLE(MshBlock* , fstream* , AllData* );
	~SolverSIMPLE();

	void solve();
	void setDT(NUM _DT);

	VCTDIMU * getNewValue();

	

private:
	fstream* pFlog;
	NUM DT;
	VCTDIMU* pFacePUVT;
	MshBlock* pMesh;
	AllData* pAllData;
	Cell* itBeginCell;
	Face* itBeginFace;
	vector<FacesInf>::iterator itBeginFacesInf;
	


	VCTDIMU* pOldCellPUVT;
	VCTDIMU* pNewCellPUVT;
	VCTDIMU* pOldNTimeCellPUVT;
	VCTDIMU* pFacePUVTFUD;
	VCTDIMU* pFacePUVTSC;//middle face uv (centrol 2-order)
	VCTDIMU* pFacePUVTRC;
	VCTDIM* pFaceGradT;
	VCTDIM* pCellGradT;
	MTDIMU_DIM* pCellGradPUVT;
	MTDIMU_DIM* pFaceGradPUVT;
	bool* etaFUDFace0;// 1 means cell0's upwind eta = 1 
	NUM* RHSBond1List;
	vector<bool>* etaCellCFUD;
	vector<NUM>* cofUPCorrect;
	vector<NUM>* cofVPCorrect;

	NUM* uStarD;
	NUM* vStarD;

	void updateFacePUVT();
	void updateFacePUVTRhieChow();
	void updateUStar();
	void updateVStar();
	void updateAllStar();
	void updateP_CorrectCofSIMPLE();//	
	void updateP_CorrectCofSIMPLEC();

	void updatePUV1();
	void updatePUV2();
	void updatePUV3(int _pst);
	void updatePUV3();
	void updatePUV4();
	void updateT();
	void updateGrad();
	void updateGradT();
	void updateT2();
	void updateNewOldPUVT();

	
};