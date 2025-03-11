#include "SolverSIMPLE.h"

SolverSIMPLE::SolverSIMPLE(MshBlock* _mesh, fstream* _pFlog, AllData* _pAllData) {
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


	pFacePUVT = pFacePUVTSC;
	etaFUDFace0 = new bool[pMesh->getNumOfFaces()];
	etaCellCFUD = new vector<bool>[pMesh->getNumOfCells()];
	cofUPCorrect = new vector<NUM>[pMesh->getNumOfCells()];
	cofVPCorrect = new vector<NUM>[pMesh->getNumOfCells()];
	RHSBond1List = new NUM[pMesh->getNumOfBond1Cells()];
	uStarD = new NUM[pMesh->getNumOfCells()];
	vStarD = new NUM[pMesh->getNumOfCells()];


	for (int i = 0; i < pMesh->getNumOfCells(); i++) {
		pOldNTimeCellPUVT[i] = pOldCellPUVT[i];
	}
	
}
SolverSIMPLE::~SolverSIMPLE() {
	delete[] etaFUDFace0; 
	delete[] etaCellCFUD;
	delete[] cofUPCorrect;
	delete[] cofVPCorrect;
	delete[] RHSBond1List;
	delete[] uStarD;
	delete[] vStarD;

}
void SolverSIMPLE::solve(){
	
	for (int pst = 0; pst < SIMPLE_INTERVAL; pst++)	{

		updateFacePUVT();
		updateGrad();
		updateFacePUVTRhieChow();//get face uf vf so on ,FUD 1st-order & SC 2nd-order
		updateUStar();
		updateVStar();
		updateFacePUVT();
		updatePUV2();
		
		updateT();
		
			
		
		VCTDIMU* oldValue = pOldNTimeCellPUVT;
		VCTDIMU* newValue = pNewCellPUVT;

		VCTDIMU residual = VCTDIMU::Zero();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
		for (int i = 0; i < pMesh->getNumOfCells(); i++) {
			for (int iR = 0; iR < DIMU; iR++) {
				if (residual[iR] < abs(newValue[i][iR] - oldValue[i][iR]) / (oldValue[i][iR] + EOR)) {					
					residual[iR] = abs(newValue[i][iR] - oldValue[i][iR]) / (oldValue[i][iR] + EOR);
				}
			}
			pOldCellPUVT[i] = pNewCellPUVT[i];
		}
		cout << "SIMPLE "<< pst <<" - residual: " << residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << endl;
		(*pFlog) << /*"Presidual: " <<*/ residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << " " << endl;
	}
}

void SolverSIMPLE::updateFacePUVT() {//first order upwind & second order 
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				Face itFace = itBeginFace[i];
				VCTDIM faceUV;
				pFacePUVTSC[i] = pNewCellPUVT[(*itFace.getBeginItPNbCells())->getId()] * itFace.getEta0() + (1 - itFace.getEta0())* pNewCellPUVT[itFace.getBeginItPNbCells()[1]->getId()];
				for (int  iDim = 0; iDim < DIM; iDim++){
					faceUV[iDim] = pFacePUVTSC[iDim + 1][1];
				}
				if (itFace.getDirectAndCells()*itFace.getDirect().transpose()*faceUV >= 0) {
					pFacePUVTFUD[i] = pNewCellPUVT[itFace.getBeginItPNbCells()[0]->getId()];
					etaFUDFace0[i] = true;
				}
				else {// upwind
					pFacePUVTFUD[i] = pNewCellPUVT[itFace.getBeginItPNbCells()[1]->getId()];
					etaFUDFace0[i] = false;
				}
			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				Face itFace = itBeginFace[i];
				VCTDIMU inletQ;
				inletQ << inletrho, inletrho*inletu, inletrho*inletv, inletE;
				pFacePUVTSC[i] << getP(inletQ), inletu, inletv, inletT;
				pFacePUVTFUD[i] = pFacePUVTSC[i];
				
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
				
				Cell* pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
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

				Face itFace = (itBeginFace)[i];
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
				Face itFace = itBeginFace[i];
				VCTDIMU inletQ;
				inletQ << inirho, inirho*iniu, inirho*iniv, iniE;

				pFacePUVTSC[i] = pNewCellPUVT[idInnerCell];
				pFacePUVTFUD[i] = pNewCellPUVT[idInnerCell];
				pFacePUVTSC[i][0] = getP(inletQ);
				pFacePUVTFUD[i][0] = getP(inletQ);



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
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++){
		auto it1NbDirect = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++){
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
void SolverSIMPLE::updateFacePUVTRhieChow() {//first order upwind & second order 
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				Face itFace = itBeginFace[i];
				pFacePUVTRC[i] = pFacePUVTSC[i];

				for (int iDim = 0; iDim < DIM; iDim++){
					pFacePUVTRC[i][iDim+1] = pFacePUVTSC[i][iDim+1] + itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[0]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[0]).getId()](iDim, 0) + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[1]).getId()] * pCellGradPUVT[(*itFace.getBeginItPNbCells()[1]).getId()](iDim, 0) - (itFace.getEta0() * itFace.getBeginItPNbCells()[0]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[0]).getId()] + (1 - itFace.getEta0()) * itFace.getBeginItPNbCells()[1]->getVolume() / uStarD[(*itFace.getBeginItPNbCells()[1]).getId()]) * pFaceGradPUVT[i](iDim, 0);
				}

			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVTRC[i] = pFacePUVTSC[i];
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVTRC[i] = pFacePUVTSC[i];
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVTRC[i] = pFacePUVTSC[i];
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

void SolverSIMPLE::updateUStar(){
	// solve U component
	SparseSolverNUM uStarSparse(pMesh->getNumOfCells(), 1, pNewCellPUVT);
	//SparseSolverNUM uStarSparse(pMesh->getNumOfCells(), 1, pNewCellPUVT);//   , 1, pNewCellPUVT
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell])*pOldNTimeCellPUVT[iCell][1] * (itBeginCell[iCell]).getVolume() / DT;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV ;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace] );
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			cofF  = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
			uStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
			RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace] [0];
			
		}
		uStarSparse.setD(cofC, iCell);
		uStarD[iCell] = cofC;
		uStarSparse.setRHSb((RHSP + RHSTime), iCell);
		
		
	}// system("pause");
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell])*(itBeginCell[iCell]).getVolume()*pOldNTimeCellPUVT[iCell][1] / DT;
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM faceUV;
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
			uStarSparse.setELE(cofF , iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
			RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][0];
		}
		uStarSparse.setD(cofC , iCell);
		uStarD[iCell] = cofC ;
		uStarSparse.setRHSb((RHSP + RHSTime) , iCell);

	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF  = 0;
		NUM RHSP  = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell])*(itBeginCell[iCell]).getVolume()*pOldNTimeCellPUVT[iCell][1] / DT;
		NUM RHSBoundary = 0;
		
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			if (idFace < pMesh->getNumOfIntFaces()){
				
				VCTDIM faceUV;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
				cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);

				uStarSparse.addELE(cofF , iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][0];
			}
			else {
				VCTDIM faceUV;
				//cout << itBeginFace[idFace].getType() << endl;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				RHSBoundary += -getPBasedRho(pFacePUVT[idFace]) * pFacePUVT[idFace][1] * faceUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];

				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][0];
			}
		}
			uStarSparse.setD(cofC, iCell);
		uStarD[iCell] = cofC;
		uStarSparse.setRHSb((RHSP + RHSTime + RHSBoundary), iCell);
		

	}

	uStarSparse.solveGS();
	auto newPX = uStarSparse.getPNewX();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++){
		pNewCellPUVT[iCell][1] = newPX[iCell];
	}


}
void SolverSIMPLE::updateVStar(){
	// solve V component
	

	SparseSolverNUM vStarSparse(pMesh->getNumOfCells(), 2, pNewCellPUVT);
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume()*pOldNTimeCellPUVT[iCell][2] / DT;
		
		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
			vStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
			RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
		}
		vStarSparse.setD(cofC, iCell);
		vStarD[iCell] = cofC ;

		vStarSparse.setRHSb((RHSP + RHSTime), iCell);


	}// system("pause");
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume()*pOldNTimeCellPUVT[iCell][2] / DT;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			VCTDIM faceUV;
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			faceUV[0] = pFacePUVT[idFace][1];
			faceUV[1] = pFacePUVT[idFace][2];
			NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
			cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
			cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
			vStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
			RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
		}
		vStarSparse.setD(cofC, iCell);
		vStarD[iCell] = cofC;
		vStarSparse.setRHSb((RHSP + RHSTime), iCell);


	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofC = getPBasedRho(pNewCellPUVT[iCell])*(itBeginCell[iCell]).getVolume() / DT;
		NUM cofF = 0;
		NUM RHSP = 0;
		NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell])*(itBeginCell[iCell]).getVolume()*pOldNTimeCellPUVT[iCell][2] / DT;
		NUM RHSBoundary = 0;

		for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
			int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			if (idFace < pMesh->getNumOfIntFaces()) {
				VCTDIM faceUV;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace])* faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
				cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
				vStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
			}
			else {
				VCTDIM faceUV;
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				RHSBoundary += -getPBasedRho(pFacePUVT[idFace]) * pFacePUVT[idFace][2] * faceUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
			}
		}

		vStarSparse.setD(cofC, iCell);
		vStarD[iCell] = cofC;
		vStarSparse.setRHSb((RHSP + RHSTime + RHSBoundary), iCell);
	}

	vStarSparse.solveGS();
	auto newPX = vStarSparse.getPNewX();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pNewCellPUVT[iCell][2] = newPX[iCell];
	}
}
void SolverSIMPLE::updateAllStar() {
	// solve All velocity component

	SparseSolverNUM vAllStarSparse(pMesh->getNumOfCells(), 2, pNewCellPUVT);

	for (int iDimStar = 0; iDimStar < DIM; iDimStar++){


		for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
			int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
			NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
			NUM cofF = 0;
			NUM RHSP = 0;
			NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() * pOldNTimeCellPUVT[iCell][2] / DT;

			for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
				VCTDIM faceUV;
				int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
				cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
				vAllStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
			}
			vAllStarSparse.setD(cofC, iCell);

			vAllStarSparse.setRHSb((RHSP + RHSTime), iCell);


		}// system("pause");
		for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
			int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
			// cout << iCell << endl;
			NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
			NUM cofF = 0;
			NUM RHSP = 0;
			NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() * pOldNTimeCellPUVT[iCell][2] / DT;

			for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
				VCTDIM faceUV;
				int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
				faceUV[0] = pFacePUVT[idFace][1];
				faceUV[1] = pFacePUVT[idFace][2];
				NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
				cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
				cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
				vAllStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
				RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
			}
			vAllStarSparse.setD(cofC, iCell);
			vAllStarSparse.setRHSb((RHSP + RHSTime), iCell);


		}
		for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {
			int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
			NUM cofC = getPBasedRho(pNewCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() / DT;
			NUM cofF = 0;
			NUM RHSP = 0;
			NUM RHSTime = getPBasedRho(pOldNTimeCellPUVT[iCell]) * (itBeginCell[iCell]).getVolume() * pOldNTimeCellPUVT[iCell][2] / DT;
			NUM RHSBoundary = 0;

			for (int iNbFace = 0; iNbFace < (itBeginCell[iCell]).getNumOfNbFaces(); iNbFace++) {
				int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
				if (idFace < pMesh->getNumOfIntFaces()) {
					VCTDIM faceUV;
					faceUV[0] = pFacePUVT[idFace][1];
					faceUV[1] = pFacePUVT[idFace][2];
					NUM  tempIntergral = getPBasedRho(pFacePUVT[idFace]) * faceUV.transpose() * ((itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace]);
					cofC += tempIntergral * itBeginCell[iCell].getBeginItEta0C()[iNbFace];// etaCellCFUD[iCell][iNbFace];
					cofF = tempIntergral * (1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]); // (1 - etaCellCFUD[iCell][iNbFace]);
					vAllStarSparse.setELE(cofF, iCell, (itBeginCell[iCell]).getBeginItPNbCells()[iNbFace]->getId());
					RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
				}
				else {
					VCTDIM faceUV;
					faceUV[0] = pFacePUVT[idFace][1];
					faceUV[1] = pFacePUVT[idFace][2];
					RHSBoundary += -getPBasedRho(pFacePUVT[idFace]) * pFacePUVT[idFace][2] * faceUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
					RHSP += -pFacePUVT[idFace][0] * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace][1];
				}
			}
			vAllStarSparse.setD(cofC, iCell);
			vAllStarSparse.setRHSb((RHSP + RHSTime + RHSBoundary), iCell);
		}

		vAllStarSparse.solveGS();
		auto newPX = vAllStarSparse.getPNewX();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
		for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
			pNewCellPUVT[iCell][2] = newPX[iCell];
		}
	}

}

void SolverSIMPLE::updatePUV1(){
	//update p' u' v'  
	// maybe a new method  :: lihangbo's
	// wrong 
	SparseSolverNUM pCorrectSparse(pMesh->getNumOfCells());
	cout << pMesh->getNumOfCells() << endl;
	cout << pCorrectSparse.getDim() << endl;
	
	// get cof of p' implicit ,with continuity equation
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		vector<NUM > cofPFList;
		for (int iNbCell = 0; iNbCell < numOfCFace; iNbCell++){
			cofPFList.push_back(0);
		}
		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0] / DT;

		NUM RHS = -itBeginCell[iCell].getVolume()* (getPBasedRho(pNewCellPUVT[iCell])-getPBasedRho(pOldNTimeCellPUVT[iCell]))/DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace])*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			cofPC += etafC*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			cofPFList[iNbFace] += ((1 - etafC)*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]));
			
			//uc' vc' part (C and F)
			for (int iCNbFace = 0; iCNbFace < numOfCFace; iCNbFace++) {
				cout << cofUPCorrect[iCell][iNbFace] << endl;
				cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][iNbFace] * directCfFace[0] + cofVPCorrect[iCell][iNbFace] * directCfFace[1]);

			}
			cofPC += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][numOfCFace] * directCfFace[0] + cofVPCorrect[iCell][numOfCFace] * directCfFace[1]);
			
			//uF' vF' part 
			int numOfFFace = itBeginCell[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()].getNumOfNbCells();
			for (int iFNbFace = 0; iFNbFace < numOfFFace; iFNbFace++) {
				int idFfFace = itBeginCell[idFCell].getBeginItPNbFaces()[iFNbFace]->getId();
				int idFfCell = itBeginCell[idFCell].getBeginItPNbCells()[iFNbFace]->getId();
				VCTDIM directFfFace = itBeginCell[idFCell].getBeginItDirectOfNbFaces()[iFNbFace];
				pCorrectSparse.addELE(getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][iFNbFace] * directCfFace[0] + cofVPCorrect[idFCell][iFNbFace] * directCfFace[1]),iCell, idFfCell);
				
			}
			cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][numOfFFace] * directCfFace[0] + cofVPCorrect[idFCell][numOfFFace] * directCfFace[1]);
			pCorrectSparse.addELE(cofPFList[iNbFace], iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	//system("pause");
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		vector<NUM > cofPFList;
		for (int iNbCell = 0; iNbCell < numOfCFace; iNbCell++) {
			cofPFList.push_back(0);
		}
		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0] / DT;

		NUM RHS = -itBeginCell[iCell].getVolume()*(getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();


			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace])*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			cofPC += etafC*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			cofPFList[iNbFace] += ((1 - etafC)*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]));

			//uc' vc' part (C and F)
			for (int iCNbFace = 0; iCNbFace < numOfCFace; iCNbFace++) {
				cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][iNbFace] * directCfFace[0] + cofVPCorrect[iCell][iNbFace] * directCfFace[1]);
			}
			cofPC += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][numOfCFace] * directCfFace[0] + cofVPCorrect[iCell][numOfCFace] * directCfFace[1]);

			//uF' vF' part 
			int numOfFFace = itBeginCell[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()].getNumOfNbCells();
			for (int iFNbFace = 0; iFNbFace < numOfFFace; iFNbFace++) {
				int idFfFace = itBeginCell[idFCell].getBeginItPNbFaces()[iFNbFace]->getId();
				int idFfCell = itBeginCell[idFCell].getBeginItPNbCells()[iFNbFace]->getId();
				VCTDIM directFfFace = itBeginCell[idFCell].getBeginItDirectOfNbFaces()[iFNbFace];
				pCorrectSparse.addELE(getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][iFNbFace] * directCfFace[0] + cofVPCorrect[idFCell][iFNbFace] * directCfFace[1]), iCell, idFfCell);
				

			}
			cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][numOfFFace] * directCfFace[0] + cofVPCorrect[idFCell][numOfFFace] * directCfFace[1]);
			pCorrectSparse.addELE(cofPFList[iNbFace], iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);


	}
	
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		vector<NUM > cofPFList;
		for (int iNbCell = 0; iNbCell < numOfCFace; iNbCell++) {
			cofPFList.push_back(0);
		}
		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0] / DT;

		NUM RHS = -itBeginCell[iCell].getVolume()*(getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();


			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace])*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			cofPC += etafC*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			cofPFList[iNbFace] += ((1 - etafC)*getPBasedRho(pNewCellPUVT[iCell])/pNewCellPUVT[iCell][0]*(pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]));
			if (itBeginFace[idNbCFace].getType() != 2){
				continue;
			}

			//uc' vc' part (C and F)
			for (int iCNbFace = 0; iCNbFace < numOfCFace; iCNbFace++) {
				cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][iNbFace] * directCfFace[0] + cofVPCorrect[iCell][iNbFace] * directCfFace[1]);
			}
			cofPC += getPBasedRho(pFacePUVT[idNbCFace])*etafC*(cofUPCorrect[iCell][numOfCFace] * directCfFace[0] + cofVPCorrect[iCell][numOfCFace] * directCfFace[1]);

			//uF' vF' part 
			int numOfFFace = itBeginCell[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()].getNumOfNbCells();
			for (int iFNbFace = 0; iFNbFace < numOfFFace; iFNbFace++) {
				int idFfFace = itBeginCell[idFCell].getBeginItPNbFaces()[iFNbFace]->getId();
				int idFfCell = itBeginCell[idFCell].getBeginItPNbCells()[iFNbFace]->getId();
				VCTDIM directFfFace = itBeginCell[idFCell].getBeginItDirectOfNbFaces()[iFNbFace];
				pCorrectSparse.addELE(getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][iFNbFace] * directCfFace[0] + cofVPCorrect[idFCell][iFNbFace] * directCfFace[1]), iCell, idFfCell);


			}
			cofPFList[iNbFace] += getPBasedRho(pFacePUVT[idNbCFace])*(1 - etafC)*(cofUPCorrect[idFCell][numOfFFace] * directCfFace[0] + cofVPCorrect[idFCell][numOfFFace] * directCfFace[1]);
			pCorrectSparse.addELE(cofPFList[iNbFace], iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}





	pCorrectSparse.solveGS();

	auto newPX = pCorrectSparse.getPNewX();

	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++){
		pNewCellPUVT[iCell][0] += newPX[iCell];
		for (int iNbCell = 0; iNbCell < itBeginCell[iCell].getNumOfNbCells();iNbCell++) {
			pNewCellPUVT[iCell][1] += cofUPCorrect[iCell][iNbCell] * newPX[itBeginCell[iCell].getBeginItPNbCells()[iNbCell]->getId()];
			pNewCellPUVT[iCell][2] += cofVPCorrect[iCell][iNbCell] * newPX[itBeginCell[iCell].getBeginItPNbCells()[iNbCell]->getId()];
		}
		pNewCellPUVT[iCell][1] += cofUPCorrect[iCell][itBeginCell[iCell].getNumOfNbCells()];
		pNewCellPUVT[iCell][2] += cofVPCorrect[iCell][itBeginCell[iCell].getNumOfNbCells()];
	}
	//system("pause");
	

}
void SolverSIMPLE::updatePUV2() {

	// china bowei's paper:updatePUV2 
	//update p' u' v'  
	SparseSolverNUM pCorrectSparse(pMesh->getNumOfCells());

	// get cof of p' implicit ,with continuity equation
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();
		
		NUM cofPC = itBeginCell[iCell].getVolume() / DT / (CP - CV) / pNewCellPUVT[iCell][DIMU - 1];

		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {
			
			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / uStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[0] / directCfFace.norm();//could be better 
			NUM cofVCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / vStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[1] / directCfFace.norm();//could be better 



			cofPC += etafC * getPBasedRho(pFacePUVT[idNbCFace])/pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			NUM cofPF = (1- etafC) * getPBasedRho(pFacePUVT[idNbCFace])/pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) *  (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}	
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = itBeginCell[iCell].getVolume() / DT / (CP - CV) / pNewCellPUVT[iCell][DIMU - 1];

		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / uStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[0] / directCfFace.norm();//could be better 
			NUM cofVCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / vStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[1] / directCfFace.norm();//could be better 

			cofPC += etafC * getPBasedRho(pFacePUVT[idNbCFace]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			NUM cofPF = (1 - etafC) * getPBasedRho(pFacePUVT[idNbCFace]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = itBeginCell[iCell].getVolume() / DT / (CP - CV) / pNewCellPUVT[iCell][DIMU - 1];

		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			
			if (itBeginFace[idNbCFace].getType() == 2){
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				//rho'( p') part //rho'( p') part 
				NUM cofUCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / uStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[0] / directCfFace.norm();//could be better 
				NUM cofVCorrectOnP = (etafC * itBeginCell[iCell].getVolume() / vStarD[iCell] + (1 - etafC) * itBeginCell[idFCell].getVolume() / uStarD[idFCell]) / (itBeginCell[iCell].getCenter() - itBeginCell[idFCell].getCenter()).norm() * directCfFace[1] / directCfFace.norm();//could be better 

				cofPC += etafC * getPBasedRho(pFacePUVT[idNbCFace]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				NUM cofPF = (1 - etafC) * getPBasedRho(pFacePUVT[idNbCFace]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

				cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
				cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

				pCorrectSparse.addELE(cofPF, iCell, idFCell);
			}
			else {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			}
			
		}
		//cout << endl;

		pCorrectSparse.setD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}

	pCorrectSparse.solveGS();
	auto newPX = pCorrectSparse.getPNewX();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pNewCellPUVT[iCell][0] += SIMPLE_ALPHA_P * newPX[iCell];

		
	}
	//system("pause");
}

void SolverSIMPLE::updatePUV3(int _pst) {

	// normal rhie-chow interpolation method //  interpolation method (not rhie-chow)
	// update p' u' v'
	SparseSolverNUM pCorrectSparse(pMesh->getNumOfCells());
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (directCfFace[0]) ;//could be better 
			NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (directCfFace[1]);//could be better 

			cofPC += etafC * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			NUM cofPF = (1 - etafC) * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (directCfFace[0]);//could be better 
			NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (directCfFace[1]);//could be better 

			cofPC += etafC * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			NUM cofPF = (1 - etafC) * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);

	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();

			if (itBeginFace[idNbCFace].getType() == 2) {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				//rho'( p') part //rho'( p') part 
				NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (directCfFace[0]);//could be better 
				NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (directCfFace[1]);//could be better 

				cofPC += etafC * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				NUM cofPF = (1 - etafC) * getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

				cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
				cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

				pCorrectSparse.addELE(cofPF, iCell, idFCell);
			}
			else {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			}

		}
		
		
		pCorrectSparse.setD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}

	pCorrectSparse.solveGS();
	auto newPX = pCorrectSparse.getPNewX();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		pNewCellPUVT[iCell][0] += SIMPLE_ALPHA_P * newPX[iCell];

		if (_pst == SIMPLE_INTERVAL) {
			for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etaCF = itBeginCell[iCell].getBeginItEta0C()[iNbFace];
				NUM tempCorrectU = SIMPLE_ALPHA_U * ((1 - etaCF) * newPX[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()] + etaCF * newPX[iCell]) / uStarD[iCell];
				NUM tempCorrectV = SIMPLE_ALPHA_U * ((1 - etaCF) * newPX[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()] + etaCF * newPX[iCell]) / vStarD[iCell];

				pNewCellPUVT[iCell][1] += -directCfFace[0] * tempCorrectU;
				pNewCellPUVT[iCell][2] += -directCfFace[1] * tempCorrectV;
			} // corrcect on u v,can be ignored
		}
	}
	//system("pause");
}
void SolverSIMPLE::updatePUV4() {

	// normal rhie-chow interpolation method //  interpolation method (not rhie-chow)
	// upwind on p' rho' FUD
	// update p' u' v'
	SparseSolverNUM pCorrectSparse(pMesh->getNumOfCells());

	// get cof of p' implicit ,with continuity equation
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][0]);//could be better 
			NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][1]);//could be better 

			NUM cofPF = 0;
			if (etaCellCFUD[iCell][iNbFace]){
				cofPC += getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			}
			else {
				cofPF -= getPBasedRho(pNewCellPUVT[idFCell]) / pNewCellPUVT[idFCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			}

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
			int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

			//explicit RHS part calculate 
			RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			//rho'( p') part //rho'( p') part 
			NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][0]);//could be better 
			NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][1]);//could be better 

			NUM cofPF = 0;
			if (etaCellCFUD[iCell][iNbFace]) {
				cofPC += getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			}
			else {
				cofPF -= getPBasedRho(pNewCellPUVT[idFCell]) / pNewCellPUVT[idFCell][0]  * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
			}

			cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
			cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

			pCorrectSparse.addELE(cofPF, iCell, idFCell);
		}
		pCorrectSparse.addD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		int numOfCFace = itBeginCell[iCell].getNumOfNbFaces();

		NUM cofPC = getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0]  / DT * itBeginCell[iCell].getVolume();
		NUM RHS = -itBeginCell[iCell].getVolume() * (getPBasedRho(pNewCellPUVT[iCell]) - getPBasedRho(pOldNTimeCellPUVT[iCell])) / DT;
		for (int iNbFace = 0; iNbFace < numOfCFace; iNbFace++) {

			int idNbCFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();

			if (itBeginFace[idNbCFace].getType() == 2) {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				NUM etafC = itBeginCell[iCell].getBeginItEta0C()[iNbFace];//etaCellCFUD[iCell][iNbFace]; 
				int idFCell = itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId();

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				//rho'( p') part //rho'( p') part 
				NUM cofUCorrectOnP = (etafC / uStarD[iCell] + (1 - etafC) / uStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][0]);//could be better 
				NUM cofVCorrectOnP = (etafC / vStarD[iCell] + (1 - etafC) / vStarD[idFCell]) * (itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace][1]);//could be better 

				NUM cofPF = 0;
				if (etaCellCFUD[iCell][iNbFace]) {
					cofPC += getPBasedRho(pNewCellPUVT[iCell]) / pNewCellPUVT[iCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				}
				else {
					cofPF -= getPBasedRho(pNewCellPUVT[idFCell]) / pNewCellPUVT[idFCell][0] * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);
				}

				cofPC += getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);
				cofPF += -getPBasedRho(pFacePUVT[idNbCFace]) * (cofUCorrectOnP * directCfFace[0] + cofVCorrectOnP * directCfFace[1]);

				pCorrectSparse.addELE(cofPF, iCell, idFCell);
			}
			else {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];

				//explicit RHS part calculate 
				RHS += -getPBasedRho(pFacePUVT[idNbCFace]) * (pFacePUVT[idNbCFace][1] * directCfFace[0] + pFacePUVT[idNbCFace][2] * directCfFace[1]);

			
			}

		}
		//cout << endl;

		pCorrectSparse.setD(cofPC, iCell);
		pCorrectSparse.setRHSb(RHS, iCell);
	}
	pCorrectSparse.solveGS();


#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		// cout << pCorrectSparse.getPNewX()[iCell] << endl;
		pNewCellPUVT[iCell][0] += SIMPLE_ALPHA_P * pCorrectSparse.getPNewX()[iCell];
		
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
			NUM etaCF = itBeginCell[iCell].getBeginItEta0C()[iNbFace];
			NUM tempCorrectU = SIMPLE_ALPHA_U * ((1 - etaCF) * pCorrectSparse.getPNewX()[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()] + etaCF * pCorrectSparse.getPNewX()[iCell]) / uStarD[iCell];
			NUM tempCorrectV = SIMPLE_ALPHA_U * ((1 - etaCF) * pCorrectSparse.getPNewX()[itBeginCell[iCell].getBeginItPNbCells()[iNbFace]->getId()] + etaCF * pCorrectSparse.getPNewX()[iCell]) / vStarD[iCell];

			pNewCellPUVT[iCell][1] += -directCfFace[0] * tempCorrectU;
			pNewCellPUVT[iCell][2] += -directCfFace[1] * tempCorrectV;
		}


	}
}


void SolverSIMPLE::setDT(NUM _DT) {
	this->DT = _DT;
}

void SolverSIMPLE::updateGrad() {
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto itDF = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		MTDIMU_DIM tempMT = MTDIMU_DIM::Zero();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
			
			for (int iDim = 0; iDim < DIM; iDim++) {
				tempMT.col(iDim) += pFacePUVTSC[iCell] * (itDF[iNbFace])(iDim);
			}
		}
		pCellGradPUVT[iCell] = tempMT / itBeginCell[iCell].getVolume();
	}
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		pFaceGradPUVT[i] = itBeginFace[i].getEta0()*pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + (1 - itBeginFace[i].getEta0())*pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[1]->getId()];
	}//simple interpolation on inter faces
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		pFaceGradPUVT[i] = pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary faces  

}





void SolverSIMPLE::updateP_CorrectCofSIMPLE() {


	//U component 
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofUCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			VCTDIM pFaceRhoUV;
			int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
			pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
			cofUCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
		}
		cofUCOfF += getPBasedRho(pNewCellPUVT[iCell])*itBeginCell[iCell].getVolume() / DT;
		NUM cofUPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			cofUPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
			cofUPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0) / cofUCOfF);
		}
		cofUPCorrect[iCell].push_back(cofUPC / cofUCOfF);


	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofUCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			VCTDIM pFaceRhoUV;
			int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
			pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
			cofUCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
		}
		cofUCOfF += getPBasedRho(pNewCellPUVT[iCell]) / DT;
		NUM cofUPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			cofUPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
			cofUPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0) / cofUCOfF);
		}
		cofUPCorrect[iCell].push_back(cofUPC / cofUCOfF);


	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofUCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			if (itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getType() == 2) {
				VCTDIM pFaceRhoUV;
				int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
				pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
				pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
				cofUCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
			}
		}
		cofUCOfF += getPBasedRho(pNewCellPUVT[iCell]) / DT;
		NUM cofUPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			if (itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getType() == 2) {
				cofUPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
				cofUPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0) / cofUCOfF);
			}
			else {
				cofUPCorrect[iCell].push_back(0);
			}
		}
		cofUPCorrect[iCell].push_back(cofUPC / cofUCOfF);


	}



	//V component 
	for (int iInnerCell = 0; iInnerCell < pMesh->getNumOfInnerCells(); iInnerCell++) {//inner cell momentum
		int iCell = pMesh->getBeginItInnerIdList()[iInnerCell];
		NUM cofVCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			VCTDIM pFaceRhoUV;
			int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
			pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
			cofVCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
		}
		cofVCOfF += getPBasedRho(pNewCellPUVT[iCell])*itBeginCell[iCell].getVolume() / DT;
		NUM cofVPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			cofVPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
			cofVPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](1) / cofVCOfF);
		}
		cofVPCorrect[iCell].push_back(cofVPC / cofVCOfF);


	}
	for (int iBond2Cell = 0; iBond2Cell < pMesh->getNumOfBond2Cells(); iBond2Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond2IdList()[iBond2Cell];
		NUM cofVCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			VCTDIM pFaceRhoUV;
			int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
			pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
			pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
			cofVCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
		}
		cofVCOfF += getPBasedRho(pNewCellPUVT[iCell]) / DT;
		NUM cofVPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			cofVPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
			cofVPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](1) / cofVCOfF);
		}
		cofVPCorrect[iCell].push_back(cofVPC / cofVCOfF);


	}
	for (int iBond1Cell = 0; iBond1Cell < pMesh->getNumOfBond1Cells(); iBond1Cell++) {//inner cell momentum
		int iCell = pMesh->getBeginItBond1IdList()[iBond1Cell];
		NUM cofVCOfF = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			if (itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getType() == 2) {
				VCTDIM pFaceRhoUV;
				int idNbFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();
				pFaceRhoUV[0] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][1];
				pFaceRhoUV[1] = getPBasedRho(pFacePUVT[idNbFace])*pFacePUVT[idNbFace][2];
				cofVCOfF += itBeginCell[iCell].getBeginItEta0C()[iNbFace] * pFaceRhoUV.transpose() * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace];
			}
		}
		cofVCOfF += getPBasedRho(pNewCellPUVT[iCell]) / DT;
		NUM cofVPC = 0;
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
			if (itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getType() == 2) {
				cofVPC += -(itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](0);
				cofVPCorrect[iCell].push_back(-(1 - itBeginCell[iCell].getBeginItEta0C()[iNbFace]) * (itBeginCell[iCell]).getBeginItDirectOfNbFaces()[iNbFace](1) / cofVCOfF);
			}
			else {
				cofVPCorrect[iCell].push_back(0);
			}
		}
		cofVPCorrect[iCell].push_back(cofVPC / cofVCOfF);



	}






}
void SolverSIMPLE::updateP_CorrectCofSIMPLEC() {

}

void SolverSIMPLE::updateNewOldPUVT() {

}
VCTDIMU* SolverSIMPLE::getNewValue() {
	return pAllData->getPNewCellPUVT();
}
void SolverSIMPLE::updateT() {

	if (FLAG_PBASED_EXPLICIT != 0) {
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
		//explicit method on energy equation 
		for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
			NUM RHS = 0;
			for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbCells(); iNbFace++) {
				VCTDIM directCfFace = itBeginCell[iCell].getBeginItDirectOfNbFaces()[iNbFace];
				int idFace = (itBeginCell[iCell]).getBeginItPNbFaces()[iNbFace]->getId();

				RHS += -getPBasedRho(pFacePUVT[idFace]) * pFacePUVT[idFace][DIMU - 1] * (directCfFace[0] * pFacePUVT[idFace][1] + directCfFace[1] * pFacePUVT[idFace][2]);
				RHS += 1 / CP * pFacePUVT[idFace][0] * (directCfFace[0] + directCfFace[1]);
				// RHS += TEMPK / CP * pFaceGradT[idFace].transpose() * directCfFace;
			}
			pNewCellPUVT[iCell][DIMU - 1] = (RHS / itBeginCell[iCell].getVolume() * DT + 1 / CP * (pNewCellPUVT[iCell][0] - pOldNTimeCellPUVT[iCell][0]) + getPBasedRho(pOldNTimeCellPUVT[iCell]) * pOldNTimeCellPUVT[iCell][DIMU - 1]) / getPBasedRho(pNewCellPUVT[iCell]);

			//if (itBeginCell[iCell].getCenter()[1] > 0.75 && itBeginCell[iCell].getCenter()[1] < 1.25 && itBeginCell[iCell].getCenter()[0] > 0.55 && itBeginCell[iCell].getCenter()[0] < 0.7) {
			//	cout <<"xy: "<< itBeginCell[iCell].getCenter() [0]<<" "<< itBeginCell[iCell].getCenter()[1] <<" " << pCorrectSparse.getPNewX()[iCell] << endl;
			//}

		}


	}
	else {
		//implicit method on energy equation 


	}


}

