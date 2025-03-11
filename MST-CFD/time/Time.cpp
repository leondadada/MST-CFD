#pragma once
#include "./Time.h"

Time::Time(MshBlock* _mesh, fstream* _pFlog,AllData* _pAllData) {
	pMesh = _mesh;
	pFlog = _pFlog;
	pAllData = _pAllData;

}
Time::~Time() {

}
void Time::initialization(VCTDIMU _iniQ) {

	for (int i = 0; i < pMesh->getNumOfCells(); i++) {
		pAllData->getP1NewCellQs()[i] = _iniQ;
		pAllData->getP1OldCellQs()[i] = _iniQ;
		pAllData->getPtP1OldNTimeCellQs()[i] = _iniQ;

		pAllData->getPCellGradPUVT()[i] = MTDIMU_DIM::Zero();

		pAllData->getPNewCellPUVT()[i] << getP(_iniQ), _iniQ[1] / _iniQ[0], _iniQ[2] / _iniQ[0], getT(_iniQ);
		pAllData->getPOldCellPUVT()[i] << getP(_iniQ), _iniQ[1] / _iniQ[0], _iniQ[2] / _iniQ[0], getT(_iniQ);
		
		pAllData->getPOldNTimeCellPUVT()[i] << getP(_iniQ), _iniQ[1] / _iniQ[0], _iniQ[2] / _iniQ[0], getT(_iniQ);
		if (pMesh->getBeginItCellsList()[i].getCenter()[0] > 0.5 ){
			auto iniQ = _iniQ;
			iniQ[0] *= 0.125;
			iniQ[DIMU - 1] *= 0.1;
			pAllData->getP1NewCellQs()[i] = iniQ;
			pAllData->getP1OldCellQs()[i] = iniQ;
			pAllData->getPtP1OldNTimeCellQs()[i] = iniQ;
			pAllData->getPNewCellPUVT()[i] << getP(iniQ), iniQ[1] / iniQ[0], iniQ[2] / iniQ[0], getT(iniQ);
			pAllData->getPOldCellPUVT()[i] << getP(iniQ), iniQ[1] / iniQ[0], iniQ[2] / iniQ[0], getT(iniQ);
			//cout << pAllData->getP1OldCellQs()[i].transpose() << endl;
		}

	}
	for (int i = 0; i < pMesh->getNumOfFaces(); i++) {

		pAllData->getP1OldFaceQs()[i] = _iniQ;

		pAllData->getPFaceGradPUVT()[i] = MTDIMU_DIM::Zero();

		pAllData->getP1OldFaceViscidFlux()[i] = MTDIMU_DIM::Zero();
		pAllData->getP1OldFaceConvectFlux()[i] = MTDIMU_DIM::Zero();

		pAllData->getPFacePUVTFUD()[i] << getP(_iniQ), _iniQ[1] / _iniQ[0], _iniQ[2] / _iniQ[0], getT(_iniQ);
		pAllData->getPFacePUVTSC()[i] << getP(_iniQ), _iniQ[1] / _iniQ[0], _iniQ[2] / _iniQ[0], getT(_iniQ);
	}


}
void Time::goNextTimeStep(){
	//choose which solver to use 
	//RhoSolver thisSolver(pMesh, pFlog, pAllData);
#if RHO_P == 0 //rho 
	RhoSolver thisSolver(pMesh, pFlog, pAllData);
#else
	PSolver thisSolver(pMesh, pFlog, pAllData);
#endif
	thisSolver.setDT(1./STEP_TIME);

#if FLAGPSUEDO == 0 //no psuedo time 
		thisSolver.solve();
		VCTDIMU* oldValue = thisSolver.getOldValue();
		VCTDIMU* newValue = thisSolver.getNewValue();
		
		VCTDIMU residual = VCTDIMU::Zero();
		for (int i = 0; i < pMesh->getNumOfCells(); i++) {
			for (int iR = 0; iR < DIMU; iR++) {
				if (residual[iR] < abs(newValue[i][iR] - oldValue[i][iR]) / oldValue[i][iR]) {
					residual[iR] = abs(newValue[i][iR] - oldValue[i][iR]) / oldValue[i][iR];
				}
			}
		}
		cout << "Time - residual: " << residual.transpose() << endl;
		(*pFlog) << fixed << setprecision(accuracyFile) << setw(accuracyFile) << /*"Presidual: " <<*/ residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << " " << endl;
		
		thisSolver.updateNewToOld();
	
#else   //with psuedo time
		VCTDIMU* oldValue = thisSolver.getOldValue();
		VCTDIMU* newValue = thisSolver.getNewValue();
		VCTDIMU* oldNTimeValue = thisSolver.getOldNTimeValue();
		for (int i = 0; i != pMesh->getNumOfCells(); i++) {
			pAllData->getPtP1OldNTimeCellQs()[i] = oldValue[i];
		}
		for (int pt = 0; pt < PT_STEP; pt++) {

			thisSolver.solve();

			VCTDIMU residual = VCTDIMU::Zero();
			for (int i = 0; i < pMesh->getNumOfCells(); i++) {
				for (int iR = 0; iR < DIMU; iR++) {
					if (residual[iR] < abs(newValue[i][iR] - oldValue[i][iR]) / oldValue[i][iR]) {
						residual[iR] = abs(newValue[i][iR] - oldValue[i][iR]) / oldValue[i][iR];
					}
				}
			}

			//cout << "PSresidual: " << residual.transpose() << endl;
			//(*pFlog) << "PSresidual: " << residual.transpose() << endl;
			thisSolver.updateNewToOld();
		}
#endif	
}

VCTDIMU * Time::getP1NewCells(){
	return pAllData->getP1NewCellQs();
}
VCTDIMU * Time::getP1OldCells(){
	return pAllData->getP1OldCellQs();
}
VCTDIMU * Time::getP1OldPUVTs() {
	return pAllData->getPOldCellPUVT();
}


