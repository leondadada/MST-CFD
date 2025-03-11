#include "AllData.h"

void AllData::createAllData(int _numOfCells, int _numOfFaces){
	numOfCells = _numOfCells;
	numOfFaces = _numOfFaces;

	p1NewCellQs = new VCTDIMU[_numOfCells];
	p1OldCellQs = new VCTDIMU[_numOfCells];
	p1OldFaceQs = new VCTDIMU[_numOfFaces];
	ptP1OldNTimeCellQs = new VCTDIMU[_numOfCells];
	
	pCellGradPUVT = new MTDIMU_DIM[_numOfCells];
	pFaceGradPUVT = new MTDIMU_DIM[_numOfFaces];

	p1OldFaceConvectFlux = new MTDIMU_DIM[_numOfFaces];
	p1OldFaceViscidFlux = new MTDIMU_DIM[_numOfFaces];
	
	pOldCellPUVT = new VCTDIMU[_numOfCells];
	pNewCellPUVT = new VCTDIMU[_numOfCells];

	//next will be better: than last

	pOldNTimeCellPUVT = new VCTDIMU[_numOfCells];
	pFacePUVTFUD = new VCTDIMU[_numOfFaces];
	pFacePUVTSC = new VCTDIMU[_numOfFaces]; //middle face uv (centrol 2-order)
	pFacePUVTRC = new VCTDIMU[_numOfFaces];
	
}
AllData::~AllData() {
	delete[] p1OldCellQs;
	delete[] p1NewCellQs;
	delete[] p1OldFaceQs;
	delete[] ptP1OldNTimeCellQs;

	delete[] pCellGradPUVT;
	delete[] pFaceGradPUVT;

	delete[] p1OldFaceConvectFlux;
	delete[] p1OldFaceViscidFlux;


	delete[] pOldCellPUVT;
	delete[] pNewCellPUVT;
	delete[] pOldNTimeCellPUVT;
	delete[] pFacePUVTFUD;
	delete[] pFacePUVTSC;
	delete[] pFacePUVTRC;

}
VCTDIMU* AllData::getP1NewCellQs() {
	return p1NewCellQs;
}
VCTDIMU* AllData::getP1OldFaceQs() {
	return p1OldFaceQs;
}
VCTDIMU* AllData::getP1OldCellQs() {
	return p1OldCellQs;
}
VCTDIMU* AllData::getPtP1OldNTimeCellQs() {
	return ptP1OldNTimeCellQs;
}
MTDIMU_DIM* AllData::getPCellGradPUVT() {
	return pCellGradPUVT;
}
MTDIMU_DIM* AllData::getPFaceGradPUVT() {
	return pFaceGradPUVT;
}
MTDIMU_DIM* AllData::getP1OldFaceConvectFlux() {
	return p1OldFaceConvectFlux;
}
MTDIMU_DIM* AllData::getP1OldFaceViscidFlux() {
	return p1OldFaceViscidFlux;
}
VCTDIMU* AllData::getPOldCellPUVT() {
	return pOldCellPUVT;
}
VCTDIMU* AllData::getPNewCellPUVT() {
	return pNewCellPUVT;
}
VCTDIMU* AllData::getPOldNTimeCellPUVT() {
	return pOldNTimeCellPUVT;
}
VCTDIMU* AllData::getPFacePUVTFUD() {
	return pFacePUVTFUD;
}
VCTDIMU* AllData::getPFacePUVTSC() {
	return pFacePUVTSC;
}
VCTDIMU* AllData::getPFacePUVTRC() {
	return pFacePUVTRC;
}

void AllData::pToFlux(){
	for (int i = 0; i < numOfCells; i++){
		p1OldCellQs[i][0] = getPBasedRho(pOldCellPUVT[i]);
		p1OldCellQs[i][1] = p1OldCellQs[i][0]* pOldCellPUVT[i](1);
		p1OldCellQs[i][2] = p1OldCellQs[i][0]* pOldCellPUVT[i](2);
		p1OldCellQs[i][3] = CV*pOldCellPUVT[i][DIMU-1] + 0.5*p1OldCellQs[i][0] * (pOldCellPUVT[i](1) * pOldCellPUVT[i](1) + pOldCellPUVT[i](2)*pOldCellPUVT[i](2));
	}
}
