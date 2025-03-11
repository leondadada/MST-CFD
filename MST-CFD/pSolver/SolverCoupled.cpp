#include "SolverCoupled.h"

SolverCoupled::SolverCoupled(MshBlock* _mesh, fstream* _pFlog, AllData* _pAllData) {
	pMesh = _mesh;
	pFlog = _pFlog;
	pAllData = _pAllData;

	itBeginCell = pMesh->getBeginItCellsList();
	itBeginFace = pMesh->getBeginItFacesList();
	itBeginFacesInf = pMesh->getBeginItFacesInfList();


	pOldCellPUVT = _pAllData->getPOldCellPUVT();
	pNewCellPUVT = _pAllData->getPNewCellPUVT();

	pOldNTimeCellPUVT = _pAllData->getPOldNTimeCellPUVT();
	pFacePUVTFUD = _pAllData->getPFacePUVTFUD();//middle face rho* uv (centrol 2-order)
	pFacePUVTSC = _pAllData->getPFacePUVTSC();
	pFacePUVTRC = _pAllData->getPFacePUVTRC();
	pCellGradPUVT = _pAllData->getPCellGradPUVT();
	pFaceGradPUVT = _pAllData->getPFaceGradPUVT();

	pFacePUVT = pFacePUVTRC;
	etaFUDFace0 = new bool[pMesh->getNumOfFaces()];
	etaCellCFUD = new vector<bool>[pMesh->getNumOfCells()];
	RHSBond1List = new NUM[pMesh->getNumOfBond1Cells()];
	uStarD = new NUM[pMesh->getNumOfCells()];
	vStarD = new NUM[pMesh->getNumOfCells()];

	for (int i = 0; i < pMesh->getNumOfCells(); i++) {
		pOldNTimeCellPUVT[i] = pOldCellPUVT[i];
	}

}
SolverCoupled::~SolverCoupled() {
	delete[] etaFUDFace0;
	delete[] etaCellCFUD;
	delete[] RHSBond1List;
	delete[] uStarD;
	delete[] vStarD;
}
void SolverCoupled::solve(){
	for (int pst = 0; pst < COUPLED_INTERVAL; pst++) {

		updateFacePUVT();//get face uf vf so on ,FUD 1st-order & SC 2nd-order
		updateFacePUVTRhieChow();
		solveCoupledUVP();
		updateT();


		VCTDIMU* oldValue = pOldCellPUVT;
		VCTDIMU* newValue = pNewCellPUVT;

		VCTDIMU residual = VCTDIMU::Zero();
		for (int i = 0; i < pMesh->getNumOfCells(); i++) {
			for (int iR = 0; iR < DIMU; iR++) {
				if (residual[iR] < abs(newValue[i][iR] - oldValue[i][iR]) / (oldValue[i][iR] + EOR)) {
					residual[iR] = abs(newValue[i][iR] - oldValue[i][iR]) / (oldValue[i][iR] + EOR);
				}
			}
			pOldCellPUVT[i] = pNewCellPUVT[i];
		}
		cout << "Coupled " << pst << "- residual: " << residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << endl;
		(*pFlog) << /*"Presidual: " <<*/ residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << " " << endl;
		if (residual.norm() < EOR2) {
			break;
		}

	}	
}
void SolverCoupled::solveCoupledUVP() {
	//my IILU solver 
	SparseSolverNUM coupledSolverUVP;
	coupledSolverUVP.reloadUVP(pMesh->getNumOfCells()*3,pNewCellPUVT);
	NUM enlargeOnP = 1e4; // to make matrix diagonally dominant 

	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T * pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T * pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];

		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
			MTDIM1 cofFX = MTDIM1::Zero();
			MTDIM1 cofFY = MTDIM1::Zero();
			cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
			cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
			cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFX(0, 2) = 1;
			cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
			cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
			cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFY(1, 2) = 1;
			cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			MTDIM1 cofF = (1 - etafC) * (directCfFace[0] * cofFX + directCfFace[1] * cofFY);
			cofC += etafC / (1 - etafC) * cofF;
			coupledSolverUVP.setELE(cofF(0, 0), 3 * iCell, 3 * idFCell);
			coupledSolverUVP.setELE(cofF(1, 1), 3 * iCell + 1, 3 * idFCell + 1);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 2), 3 * iCell + 2, 3 * idFCell + 2);
			coupledSolverUVP.setELE(cofF(0, 2), 3 * iCell, 3 * idFCell + 2);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 0), 3 * iCell + 2, 3 * idFCell);
			coupledSolverUVP.setELE(cofF(1, 2), 3 * iCell + 1, 3 * idFCell + 2);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 1), 3 * iCell + 2, 3 * idFCell + 1);
			RHSf += getPBasedRho(pFacePUVT[idNbCFace]) * (directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
		}
		coupledSolverUVP.setELE( cofC(0, 0),3 * iCell, 3 * iCell);
		coupledSolverUVP.setELE( cofC(1, 1),3 * iCell + 1, 3 * iCell + 1);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 2),3 * iCell + 2, 3 * iCell + 2);
		coupledSolverUVP.setELE( cofC(0, 2),3 * iCell, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 0),3 * iCell + 2, 3 * iCell);
		coupledSolverUVP.setELE( cofC(1, 2),3 * iCell + 1, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 1),3 * iCell + 2, 3 * iCell + 1);


		coupledSolverUVP.setRHSb(V_T * (getPBasedRho(pOldCellPUVT[iCell]) * pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][1]), 3 * iCell);
		coupledSolverUVP.setRHSb(V_T * (getPBasedRho(pOldCellPUVT[iCell]) * pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][2]), 3 * iCell + 1);
		coupledSolverUVP.setRHSb((enlargeOnP * V_T * getPBasedRho(pOldCellPUVT[iCell]) + RHSf), 3 * iCell + 2);

	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T * pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T * pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];

		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
			MTDIM1 cofFX = MTDIM1::Zero();
			MTDIM1 cofFY = MTDIM1::Zero();
			cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
			cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
			cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFX(0, 2) = 1;
			cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
			cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
			cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFY(1, 2) = 1;
			cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			MTDIM1 cofF = (1 - etafC) * (directCfFace[0] * cofFX + directCfFace[1] * cofFY);
			cofC += etafC / (1 - etafC) * cofF;
			coupledSolverUVP.setELE(cofF(0, 0), 3 * iCell, 3 * idFCell);
			coupledSolverUVP.setELE(cofF(1, 1), 3 * iCell + 1, 3 * idFCell + 1);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 2), 3 * iCell + 2, 3 * idFCell + 2);
			coupledSolverUVP.setELE(cofF(0, 2), 3 * iCell, 3 * idFCell + 2);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 0), 3 * iCell + 2, 3 * idFCell);
			coupledSolverUVP.setELE(cofF(1, 2), 3 * iCell + 1, 3 * idFCell + 2);
			coupledSolverUVP.setELE(enlargeOnP * cofF(2, 1), 3 * iCell + 2, 3 * idFCell + 1);
			RHSf += getPBasedRho(pFacePUVT[idNbCFace]) * (directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
		}
		coupledSolverUVP.setELE(cofC(0, 0), 3 * iCell, 3 * iCell);
		coupledSolverUVP.setELE(cofC(1, 1), 3 * iCell + 1, 3 * iCell + 1);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 2), 3 * iCell + 2, 3 * iCell + 2);
		coupledSolverUVP.setELE(cofC(0, 2), 3 * iCell, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 0), 3 * iCell + 2, 3 * iCell);
		coupledSolverUVP.setELE(cofC(1, 2), 3 * iCell + 1, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP * cofC(2, 1), 3 * iCell + 2, 3 * iCell + 1);


		coupledSolverUVP.setRHSb(V_T * (getPBasedRho(pOldCellPUVT[iCell]) * pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][1]), 3 * iCell);
		coupledSolverUVP.setRHSb(V_T * (getPBasedRho(pOldCellPUVT[iCell]) * pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][2]), 3 * iCell + 1);
		coupledSolverUVP.setRHSb(enlargeOnP* (V_T * getPBasedRho(pOldCellPUVT[iCell]) + RHSf), 3 * iCell + 2);


	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		VCTDIM1 RHSBound = VCTDIM1::Zero();
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T * getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T * pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T * pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];

		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			if (itBeginFace[idNbCFace].getType() == 2) {

				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
				MTDIM1 cofFX = MTDIM1::Zero();
				MTDIM1 cofFY = MTDIM1::Zero();
				cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
				cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
				cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFX(0, 2) = 1;
				cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
				cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
				cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFY(1, 2) = 1;
				cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				MTDIM1 cofF = (1 - etafC) * (directCfFace[0] * cofFX + directCfFace[1] * cofFY);
				cofC += etafC / (1 - etafC) * cofF;
				coupledSolverUVP.setELE(cofF(0, 0), 3 * iCell, 3 * idFCell);
				coupledSolverUVP.setELE(cofF(1, 1), 3 * iCell + 1, 3 * idFCell + 1);
				coupledSolverUVP.setELE(enlargeOnP* cofF(2, 2), 3 * iCell + 2, 3 * idFCell + 2);
				coupledSolverUVP.setELE(cofF(0, 2), 3 * iCell, 3 * idFCell + 2);
				coupledSolverUVP.setELE(enlargeOnP* cofF(2, 0), 3 * iCell + 2, 3 * idFCell);
				coupledSolverUVP.setELE(cofF(1, 2), 3 * iCell + 1, 3 * idFCell + 2);
				coupledSolverUVP.setELE(enlargeOnP* cofF(2, 1), 3 * iCell + 2, 3 * idFCell + 1);
				RHSf += getPBasedRho(pFacePUVT[idNbCFace]) * (directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);

			}
			else {

				MTDIM1 cofFX = MTDIM1::Zero();
				MTDIM1 cofFY = MTDIM1::Zero();
				VCTDIM1 UVP;
				UVP << pFacePUVT[idNbCFace][1], pFacePUVT[idNbCFace][2], pFacePUVT[idNbCFace][0];
				cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
				cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][1];
				cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFX(0, 2) = 1;
				cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
				cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace]) * pFacePUVT[idNbCFace][2];
				cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFY(1, 2) = 1;
				cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];
				RHSf += getPBasedRho(pFacePUVT[idNbCFace]) * (directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
				RHSBound -= (directCfFace[0] * cofFX + directCfFace[1] * cofFY) * UVP;

			}
		}
		coupledSolverUVP.setELE(cofC(0, 0), 3 * iCell, 3 * iCell);
		coupledSolverUVP.setELE(cofC(1, 1), 3 * iCell + 1, 3 * iCell + 1);
		coupledSolverUVP.setELE(enlargeOnP* cofC(2, 2), 3 * iCell + 2, 3 * iCell + 2);
		coupledSolverUVP.setELE(cofC(0, 2), 3 * iCell, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP* cofC(2, 0), 3 * iCell + 2, 3 * iCell);
		coupledSolverUVP.setELE(cofC(1, 2), 3 * iCell + 1, 3 * iCell + 2);
		coupledSolverUVP.setELE(enlargeOnP* cofC(2, 1), 3 * iCell + 2, 3 * iCell + 1);


		coupledSolverUVP.setRHSb(V_T* (getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][1]), 3 * iCell);
		coupledSolverUVP.setRHSb(V_T* (getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell]) * pNewCellPUVT[iCell][2]), 3 * iCell + 1);
		coupledSolverUVP.setRHSb(enlargeOnP*( V_T* getPBasedRho(pOldCellPUVT[iCell]) + RHSf), 3 * iCell + 2);

	}
	coupledSolverUVP.solveILUSGS();
		
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pNewCellPUVT[iCell][0] = coupledSolverUVP.getPNewX()[iCell * 3 + 2];//P
		pNewCellPUVT[iCell][1] = coupledSolverUVP.getPNewX()[iCell * 3];//U
		pNewCellPUVT[iCell][2] = coupledSolverUVP.getPNewX()[iCell * 3 + 1];//V

	}

}
void SolverCoupled::solveCoupledUVP2() {
	//Eigen solver 
	vector<Triplet<NUM>> cof;
	VectorXd RHSb(pMesh->getNumOfCells()*3);

	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T*pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T*pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		

		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
			MTDIM1 cofFX = MTDIM1::Zero();
			MTDIM1 cofFY = MTDIM1::Zero();
			cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
			cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
			cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFX(0, 2) = 1;
			cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
			cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
			cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFY(1, 2) = 1;
			cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];
			
			MTDIM1 cofF = (1- etafC)*(directCfFace[0] * cofFX + directCfFace[1] * cofFY);
			cofC += etafC/ (1 - etafC)*cofF;
			cof.emplace_back(3 * iCell, 3 * idFCell, cofF(0, 0));
			cof.emplace_back( 3 * iCell + 1, 3 * idFCell + 1, cofF(1, 1));
			cof.emplace_back( 3 * iCell + 2, 3 * idFCell + 2,cofF(2, 2));
			cof.emplace_back( 3 * iCell, 3 * idFCell + 2,cofF(0, 2));
			cof.emplace_back( 3 * iCell + 2, 3 * idFCell, cofF(2, 0));
			cof.emplace_back( 3 * iCell + 1, 3 * idFCell + 2, cofF(1, 2));
			cof.emplace_back(3 * iCell + 2, 3 * idFCell + 1, cofF(2, 1));
			RHSf += getPBasedRho(pFacePUVT[idNbCFace])*(directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
		}
		cof.emplace_back(3 * iCell, 3 * iCell, cofC(0, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 1, cofC(1, 1));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 2, cofC(2, 2));
		cof.emplace_back(3 * iCell, 3 * iCell + 2, cofC(0, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell, cofC(2, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 2, cofC(1, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 1, cofC(2, 1));


		RHSb(3 * iCell) = (V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][1]));
		RHSb(3 * iCell + 1) = (V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][2]));
		RHSb(3 * iCell + 2) = (V_T*getPBasedRho(pOldCellPUVT[iCell])+RHSf);


		
	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T*pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T*pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];


		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
			MTDIM1 cofFX = MTDIM1::Zero();
			MTDIM1 cofFY = MTDIM1::Zero();
			cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
			cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
			cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFX(0, 2) = 1;
			cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
			cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
			cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
			cofFY(1, 2) = 1;
			cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

			MTDIM1 cofF = (1 - etafC)*(directCfFace[0] * cofFX + directCfFace[1] * cofFY);
			cofC += etafC / (1 - etafC)*cofF;
			cof.emplace_back(3 * iCell, 3 * idFCell, cofF(0, 0));
			cof.emplace_back(3 * iCell + 1, 3 * idFCell + 1, cofF(1, 1));
			cof.emplace_back(3 * iCell + 2, 3 * idFCell + 2, cofF(2, 2));
			cof.emplace_back(3 * iCell, 3 * idFCell + 2, cofF(0, 2));
			cof.emplace_back(3 * iCell + 2, 3 * idFCell, cofF(2, 0));
			cof.emplace_back(3 * iCell + 1, 3 * idFCell + 2, cofF(1, 2));
			cof.emplace_back(3 * iCell + 2, 3 * idFCell + 1, cofF(2, 1));
			RHSf += getPBasedRho(pFacePUVT[idNbCFace])*(directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
		}
		cof.emplace_back(3 * iCell, 3 * iCell, cofC(0, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 1, cofC(1, 1));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 2, cofC(2, 2));
		cof.emplace_back(3 * iCell, 3 * iCell + 2, cofC(0, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell, cofC(2, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 2, cofC(1, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 1, cofC(2, 1));


		RHSb(3 * iCell) = (V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][1]));
		RHSb(3 * iCell + 1) = (V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][2]));
		RHSb(3 * iCell + 2) = (V_T*getPBasedRho(pOldCellPUVT[iCell]) + RHSf);

		
	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		MTDIM1 cofC = MTDIM1::Zero();
		NUM RHSf = 0;
		VCTDIM1 RHSBound = VCTDIM1::Zero();
		NUM V_T = pMesh->getBeginItCellsList()[iCell].getVolume() / DT;
		cofC(0, 0) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(1, 1) = V_T*getPBasedRho(pNewCellPUVT[iCell]);
		cofC(0, 2) = V_T*pNewCellPUVT[iCell][1] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(1, 2) = V_T*pNewCellPUVT[iCell][2] / RGAS / pNewCellPUVT[iCell][DIMU - 1];
		cofC(2, 2) = V_T / RGAS / pNewCellPUVT[iCell][DIMU - 1];

		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			if (itBeginFace[idNbCFace].getType() == 2) {
								
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();
				MTDIM1 cofFX = MTDIM1::Zero();
				MTDIM1 cofFY = MTDIM1::Zero();
				cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
				cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
				cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFX(0, 2) = 1;
				cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
				cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
				cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFY(1, 2) = 1;
				cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				MTDIM1 cofF = (1 - etafC)*(directCfFace[0] * cofFX + directCfFace[1] * cofFY);
				cofC += etafC / (1 - etafC)*cofF;
				cof.emplace_back(3 * iCell, 3 * idFCell, cofF(0, 0));
				cof.emplace_back(3 * iCell + 1, 3 * idFCell + 1, cofF(1, 1));
				cof.emplace_back(3 * iCell + 2, 3 * idFCell + 2, cofF(2, 2));
				cof.emplace_back(3 * iCell, 3 * idFCell + 2, cofF(0, 2));
				cof.emplace_back(3 * iCell + 2, 3 * idFCell, cofF(2, 0));
				cof.emplace_back(3 * iCell + 1, 3 * idFCell + 2, cofF(1, 2));
				cof.emplace_back(3 * iCell + 2, 3 * idFCell + 1, cofF(2, 1));
				RHSf += getPBasedRho(pFacePUVT[idNbCFace])*(directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);

			}
			else {

				MTDIM1 cofFX = MTDIM1::Zero();
				MTDIM1 cofFY = MTDIM1::Zero();
				VCTDIM1 UVP;
				UVP << pFacePUVT[idNbCFace][1], pFacePUVT[idNbCFace][2], pFacePUVT[idNbCFace][0];
				cofFX(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
				cofFX(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][1];
				cofFX(2, 0) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFX(0, 2) = 1;
				cofFX(2, 2) = pFacePUVT[idNbCFace][1] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];

				cofFY(0, 0) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
				cofFY(1, 1) = getPBasedRho(pFacePUVT[idNbCFace])*pFacePUVT[idNbCFace][2];
				cofFY(2, 1) = getPBasedRho(pFacePUVT[idNbCFace]);
				cofFY(1, 2) = 1;
				cofFY(2, 2) = pFacePUVT[idNbCFace][2] / RGAS / pFacePUVT[idNbCFace][DIMU - 1];
				RHSf += getPBasedRho(pFacePUVT[idNbCFace])*(directCfFace[0] * pFacePUVT[idNbCFace][1] + directCfFace[1] * pFacePUVT[idNbCFace][2]);
				RHSBound -= (directCfFace[0] * cofFX + directCfFace[1] * cofFY)*UVP;
				
			}
		}
		cof.emplace_back(3 * iCell, 3 * iCell, cofC(0, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 1, cofC(1, 1));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 2, cofC(2, 2));
		cof.emplace_back(3 * iCell, 3 * iCell + 2, cofC(0, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell, cofC(2, 0));
		cof.emplace_back(3 * iCell + 1, 3 * iCell + 2, cofC(1, 2));
		cof.emplace_back(3 * iCell + 2, 3 * iCell + 1, cofC(2, 1));
		
		RHSb(3 * iCell) = V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][1] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][1] + RHSBound[0]);
		RHSb(3 * iCell + 1) = V_T*(getPBasedRho(pOldCellPUVT[iCell])* pOldCellPUVT[iCell][2] + getPBasedRho(pNewCellPUVT[iCell])* pNewCellPUVT[iCell][2] + RHSBound[1]);
		RHSb(3 * iCell + 2) = V_T*getPBasedRho(pOldCellPUVT[iCell]) + RHSf + RHSBound[2];
		
	}

	SparseMatrix<NUM> A(3 * pMesh->getNumOfCells(), 3 * pMesh->getNumOfCells());
	A.setFromTriplets(cof.begin(), cof.end());
	SPARSESOLVER<SparseMatrix<NUM>> PUVSolver(A);//how to solve ?
	VectorXd x = PUVSolver.solve(RHSb);
	//cout << x << endl;
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pNewCellPUVT[iCell][0] = x[iCell * 3 + 2];//P
		pNewCellPUVT[iCell][1] = x[iCell * 3];//U
		pNewCellPUVT[iCell][2] = x[iCell * 3 + 1];//V

	}
	
}
void SolverCoupled::updateFacePUVT() {//first order upwind & second order 
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				auto itFace = itBeginFace[i];
				VCTDIM faceUV;
				pFacePUVTSC[i] = pNewCellPUVT[(*itFace.getBeginItPNbCells())->getId()] * itFace.getEta0() + (1 - itFace.getEta0())* pNewCellPUVT[itFace.getBeginItPNbCells()[1]->getId()];
				faceUV[0] = pFacePUVTSC[i][1];
				faceUV[1] = pFacePUVTSC[i][2];

				if (itFace.getDirectAndCells()*itFace.getDirect().transpose()*faceUV >= 0) {
					pFacePUVTFUD[i] = pNewCellPUVT[itFace.getBeginItPNbCells()[0]->getId()];
					etaFUDFace0[i] = true;
				}
				else {
					pFacePUVTFUD[i] = pNewCellPUVT[itFace.getBeginItPNbCells()[1]->getId()];
					etaFUDFace0[i] = false;
				}
			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				VCTDIMU inletQ;
				inletQ << inletrho, inletrho*inletu, inletrho*inletv, inletE;
				pFacePUVTFUD[i] << getP(inletQ), inletu, inletv, inletT;
				pFacePUVTSC[i] << getP(inletQ), inletu, inletv, inletT;

				auto itFace = itBeginFace[i];
				if (itFace.getDirectAndCells()*(itFace.getDirect()(0)*pFacePUVTSC[i][1] + itFace.getDirect()(1)*pFacePUVTSC[i][2]) >= 0) {
					etaFUDFace0[i] = true;
				}
				else {
					etaFUDFace0[i] = false;
				}
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {

				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellPUVT = pNewCellPUVT[pCellBoundary->getId()];
				VCTDIMU outCellPUVT = inCellPUVT;
				NUM dX = itBeginFace[i].getBeginItPNbNodes()[1]->getCenter()(0) - itBeginFace[i].getBeginItPNbNodes()[0]->getCenter()(0);
				NUM dY = itBeginFace[i].getBeginItPNbNodes()[1]->getCenter()(1) - itBeginFace[i].getBeginItPNbNodes()[0]->getCenter()(1);

				if (FLAG_WALLSMOOTH == 1) {// if wall smooth is avaliable
					outCellPUVT(2) = (2 * dX * dY * inCellPUVT(1) - (dX * dX - dY * dY) * inCellPUVT(2)) / (dX * dX + dY * dY);
					outCellPUVT(1) = (2 * dX * dY * inCellPUVT(2) - (dY * dY - dX * dX) * inCellPUVT(1)) / (dX * dX + dY * dY);
				}
				else {
					outCellPUVT(2) = 0;
					outCellPUVT(1) = 0;
				}

				int idInnerCell = (*(itBeginFace[i]).getBeginItPNbCells())->getId();
				pFacePUVTSC[i] = (outCellPUVT + inCellPUVT) / 2;
				pFacePUVTFUD[i] = (outCellPUVT + inCellPUVT) / 2;

				auto itFace = (itBeginFace)[i];
				if (itFace.getDirectAndCells()*(itFace.getDirect()(0)*pFacePUVTSC[i][1] + itFace.getDirect()(1)*pFacePUVTSC[i][2]) >= 0) {
					etaFUDFace0[i] = true;
				}
				else {
					etaFUDFace0[i] = false;
				}
				VCTDIM uvFace;
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				int idInnerCell = (*(itBeginFace[i]).getBeginItPNbCells())->getId();
				pFacePUVTSC[i] = pNewCellPUVT[idInnerCell];
				pFacePUVTFUD[i] = pNewCellPUVT[idInnerCell];

				auto itFace = (itBeginFace)[i];
				if (itFace.getDirectAndCells()*(itFace.getDirect()(0)*pFacePUVTSC[i][1] + itFace.getDirect()(1)*pFacePUVTSC[i][2]) >= 0) {
					etaFUDFace0[i] = true;
				}
				else {
					etaFUDFace0[i] = false;
				}
			}
			continue;
		}
		}

	}//all face (interior & boundary) update
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto it1NbDirect = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
			int idNbFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			if (it1NbDirect[iNbFace][0] * pFacePUVTSC[idNbFace][1] + it1NbDirect[iNbFace][1] * pFacePUVTSC[idNbFace][2] >0) {
				etaCellCFUD[iCell].push_back(true);
			}
			else {
				etaCellCFUD[iCell].push_back(false);
			}
		}
	}
}
void SolverCoupled::updateUStar() {
	// solve U component
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * pOldNTimeCellPUVT[iCell][1] * (itBeginCell[iCell]).getVolume() / DT;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
		}
		uStarD[iCell] = cofC;

	}// system("pause");
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() * pOldNTimeCellPUVT[iCell][1] / DT;
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM faceUV;
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
		}
		uStarD[iCell] = cofC;

	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() * pOldNTimeCellPUVT[iCell][1] / DT;
		NUM RHSBoundary = 0;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			if (idFace < pMesh->getNumOfIntFaces()) {

				VCTDIM faceUV;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			}
			
		}		
		uStarD[iCell] = cofC;
		
	}

}
void SolverCoupled::updateVStar() {
	// solve V component
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
		
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
		}
		vStarD[iCell] = cofC;

	}// system("pause");
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
		}
		vStarD[iCell] = cofC;


	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
		
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			if (idFace < pMesh->getNumOfIntFaces()) {
				VCTDIM faceUV;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			}
			
		}		
		vStarD[iCell] = cofC;
	}

}
void SolverCoupled::updateFacePUVTRhieChow() {//first order upwind & second order 
	updateVStar();
	updateUStar();
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				auto itFace = itBeginFace[i];
				pFacePUVTRC[i] = pFacePUVTSC[i];
				pFacePUVTRC[i][1] = pFacePUVTSC[i][1] + itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[0]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[0]).getId()](0, 0) + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[1]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[1]).getId()](0, 0) - (itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[0]).getId()] + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[1]).getId()]) * pFaceGradPUVT[i](0, 0);
				pFacePUVTRC[i][2] = pFacePUVTSC[i][2] + itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / vStarD[(*itFace.getBeginItPNbCells()[0]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[0]).getId()](1, 0) + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / vStarD[(*itFace.getBeginItPNbCells()[1]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[1]).getId()](1, 0) - (itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / vStarD[(*itFace.getBeginItPNbCells()[0]).getId()] + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / vStarD[(*itFace.getBeginItPNbCells()[1]).getId()]) * pFaceGradPUVT[i](1, 0);

			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				VCTDIMU inletQ;
				inletQ << inletrho, inletrho* inletu, inletrho* inletv, inletE;
				pFacePUVTRC[i] << getP(inletQ), inletu, inletv, inletT;
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {

				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellPUVT = pNewCellPUVT[pCellBoundary->getId()];
				VCTDIMU outCellPUVT = inCellPUVT;
				NUM dX = itBeginFace[i].getBeginItPNbNodes()[1]->getCenter()(0) - itBeginFace[i].getBeginItPNbNodes()[0]->getCenter()(0);
				NUM dY = itBeginFace[i].getBeginItPNbNodes()[1]->getCenter()(1) - itBeginFace[i].getBeginItPNbNodes()[0]->getCenter()(1);

				if (FLAG_WALLSMOOTH == 1) {// if wall smooth is avaliable
					outCellPUVT(2) = (2 * dX * dY * inCellPUVT(1) - (dX * dX - dY * dY) * inCellPUVT(2)) / (dX * dX + dY * dY);
					outCellPUVT(1) = (2 * dX * dY * inCellPUVT(2) - (dY * dY - dX * dX) * inCellPUVT(1)) / (dX * dX + dY * dY);
				}
				else {
					outCellPUVT(2) = 0;
					outCellPUVT(1) = 0;
				}
				int idInnerCell = (*(itBeginFace[i]).getBeginItPNbCells())->getId();
				pFacePUVTRC[i] = (outCellPUVT + inCellPUVT) / 2;
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				int idInnerCell = (*(itBeginFace[i]).getBeginItPNbCells())->getId();
				pFacePUVTRC[i] = pNewCellPUVT[idInnerCell];
			}
			continue;
		}
		}

	}//all face (interior & boundary) update
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto it1NbDirect = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
			int idNbFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			if (it1NbDirect[iNbFace][0] * pFacePUVTSC[idNbFace][1] + it1NbDirect[iNbFace][1] * pFacePUVTSC[idNbFace][2] > 0) {
				etaCellCFUD[iCell].push_back(true);
			}
			else {
				etaCellCFUD[iCell].push_back(false);
			}
		}
	}
}

void SolverCoupled::updateT() {


	if (FLAG_PBASED_EXPLICIT != 0) {
		//explicit method on energy equation 
#pragma omp parallel for num_threads(NUM_CPU_THREADS)		
		for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
			NUM RHS = 0;
			for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();

				RHS += -getPBasedRho(pFacePUVT[idFace]) * pFacePUVT[idFace][DIMU - 1] * (directCfFace[0] * pFacePUVT[idFace][1] + directCfFace[1] * pFacePUVT[idFace][2]);
				RHS += 1 / CP * pFacePUVT[idFace][0] * (directCfFace[0] + directCfFace[1]);
			}
			pNewCellPUVT[iCell][DIMU - 1] = (RHS / itBeginCell[iCell].getVolume() * DT + 1 / CP * (pNewCellPUVT[iCell][0] - pOldNTimeCellPUVT[iCell][0]) + getPBasedRho(pOldNTimeCellPUVT[iCell]) * pOldNTimeCellPUVT[iCell][DIMU - 1]) / getPBasedRho(pNewCellPUVT[iCell]);

		}


	}
	else {
		//implicit method on energy equation 


	}

}
void SolverCoupled::updateGrad() {
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto itDF = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		MTDIM_DIMU tempMT = MTDIM_DIMU::Zero();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {

			for (int d = 0; d < DIM; d++) {
				tempMT.row(d) += pFacePUVTSC[iCell] * (itDF[iNbFace])(d);
			}
		}
		pCellGradPUVT[iCell] = tempMT / itBeginCell[iCell].getVolume();
	}
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		pFaceGradPUVT[i] = itBeginFace[i].getEta0() * pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + (1 - itBeginFace[i].getEta0()) * pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[1]->getId()];
	}//simple interpolation on inter faces
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		pFaceGradPUVT[i] = pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary faces

}

void SolverCoupled::setDT(NUM _DT) {
	this->DT = _DT;
}
VCTDIMU* SolverCoupled::getNewValue() {
	return pAllData->getPNewCellPUVT();
}


