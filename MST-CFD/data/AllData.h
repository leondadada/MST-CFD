#pragma once
#include "../work/FUNCTION.h"
#include "../include/CONST.h"

class AllData {

public:
	~AllData();
	void createAllData(int, int);
	VCTDIMU* getP1NewCellQs();
	VCTDIMU* getP1OldFaceQs();
	VCTDIMU* getP1OldCellQs();
	VCTDIMU* getPtP1OldNTimeCellQs();

	MTDIMU_DIM * getPCellGradPUVT();
	MTDIMU_DIM * getPFaceGradPUVT();

	MTDIMU_DIM* getP1OldFaceConvectFlux();
	MTDIMU_DIM* getP1OldFaceViscidFlux();


	VCTDIMU* getPOldCellPUVT();
	VCTDIMU* getPNewCellPUVT();


	VCTDIMU* getPOldNTimeCellPUVT();
	VCTDIMU* getPFacePUVTFUD() ;
	VCTDIMU* getPFacePUVTSC();
	VCTDIMU* getPFacePUVTRC();

	void pToFlux();
private:
	int numOfCells;
	int numOfFaces;
	VCTDIMU* p1NewCellQs;
	VCTDIMU* p1OldFaceQs;
	VCTDIMU* p1OldCellQs;
	VCTDIMU* ptP1OldNTimeCellQs;
	MTDIMU_DIM* p1OldFaceConvectFlux;
	MTDIMU_DIM* p1OldFaceViscidFlux;


	VCTDIMU* pOldCellPUVT;
	VCTDIMU* pNewCellPUVT;
	MTDIMU_DIM* pCellGradPUVT;
	MTDIMU_DIM* pFaceGradPUVT;// 
	
	VCTDIMU* pOldNTimeCellPUVT;
	VCTDIMU* pFacePUVTFUD;
	VCTDIMU* pFacePUVTSC;//middle face uv (centrol 2-order)
	VCTDIMU* pFacePUVTRC;

};