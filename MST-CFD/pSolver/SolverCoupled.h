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

class SolverCoupled {

	public:
		SolverCoupled(MshBlock*, fstream*, AllData*);
		~SolverCoupled();
		void solve();

		void setDT(NUM _DT);
		VCTDIMU * getNewValue();

	private:
		fstream* pFlog;
		NUM DT;
		VCTDIMU* pFacePUVT;
		MshBlock* pMesh;
		AllData* pAllData;
		vector<Cell>::iterator itBeginCell;
		vector<Face>::iterator itBeginFace;
		vector<FacesInf>::iterator itBeginFacesInf;



		VCTDIMU* pOldCellPUVT;
		VCTDIMU* pNewCellPUVT;

		VCTDIMU* pOldNTimeCellPUVT;
		VCTDIMU* pFacePUVTFUD;
		VCTDIMU* pFacePUVTSC;//middle face uv (centrol 2-order)
		VCTDIMU* pFacePUVTRC;
		VCTDIM* pFaceGradT;
		VCTDIM* pCellGradT;
		MTDIM_DIMU* pCellGradPUVT;
		MTDIM_DIMU* pFaceGradPUVT;

		bool* etaFUDFace0;// 1 means cell0's upwind eta = 1 
		NUM* RHSBond1List;
		vector<bool>* etaCellCFUD;

		NUM* uStarD;
		NUM* vStarD;

		void updateFacePUVT();
		void updateUStar();
		void updateVStar();
		void updateFacePUVTRhieChow();
		void solveCoupledUVP();
		void solveCoupledUVP2();
		void updateT();
		void updateGradT();

		void updateGrad();


};