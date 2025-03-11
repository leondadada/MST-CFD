#include "PSolver.h"

PSolver::PSolver(MshBlock* _mesh, fstream* _pFlog, AllData* _pAllData) {
	pMesh = _mesh;
	pFlog = _pFlog;
	pAllData = _pAllData;

	itBeginCell = pMesh->getBeginItCellsList();
	itBeginFace = pMesh->getBeginItFacesList();
	itBeginFacesInf = pMesh->getBeginItFacesInfList();

}
PSolver::~PSolver(){

}
void PSolver::setDT(NUM _DT) {
	this->DT = _DT;
}

void PSolver::solve(){

	PSOLVER ppSolver(pMesh, pFlog, pAllData);
	ppSolver.setDT(DT);
	ppSolver.solve();
	
}

VCTDIMU* PSolver::getOldValue() {
	return pAllData->getPOldNTimeCellPUVT();
}
VCTDIMU* PSolver::getNewValue(){
	return pAllData->getPNewCellPUVT();
}
VCTDIMU* PSolver::getOldNTimeValue() {
	return pAllData->getPOldNTimeCellPUVT();
}
void PSolver::updateNewToOld() {
	auto pOldCellPUVT = pAllData->getPOldCellPUVT();
	auto pNewCellPUVT = pAllData->getPNewCellPUVT();
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pOldCellPUVT[iCell] = pNewCellPUVT[iCell];
	}
}


