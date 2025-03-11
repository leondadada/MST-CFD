#include "RhoSolver.h"

RhoSolver::RhoSolver(MshBlock* _mesh, fstream* _pFlog, AllData* _pAllData){
	pMesh = _mesh;
	pFlog = _pFlog;
	pAllData = _pAllData;
	itBeginCell = pMesh->getBeginItCellsList();
	itBeginFace = pMesh->getBeginItFacesList();
	itBeginFacesInf = pMesh->getBeginItFacesInfList();

	p1NewCellQs = _pAllData->getP1NewCellQs();
	p1OldFaceQs = _pAllData->getP1OldFaceQs();
	p1OldCellQs = _pAllData->getP1OldCellQs();
	ptP1OldNTimeCellQs = _pAllData->getPtP1OldNTimeCellQs();

	pCellGradPUVT = _pAllData->getPCellGradPUVT();
	pFaceGradPUVT = _pAllData->getPFaceGradPUVT();


	p1OldFaceConvectFlux = _pAllData->getP1OldFaceConvectFlux();
	p1OldFaceViscidFlux = _pAllData->getP1OldFaceViscidFlux();

	p1NewFaceGradFlux = new MTDIMU_DIM[pMesh->getNumOfFaces()];
	p1NewCellGradFlux = new MTDIMU_DIM[pMesh->getNumOfCells()];

}

RhoSolver::~RhoSolver(){
	delete[] p1NewFaceGradFlux;
	delete[] p1NewCellGradFlux;
}

void RhoSolver::setDT(NUM _DT){
	this->DT = _DT;
}

void RhoSolver::solve(){
	
#if ACCURACY == 1
	updateFaceFlux();
#elif ACCURACY == 2
	updateFaceFlux2ndOrder();
#endif

#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int i = 0; i < pMesh->getNumOfCells(); i++) {
		VCTDIMU integralConvectFlux = VCTDIMU::Zero();

		auto itDF = itBeginCell[i].getBeginItDirectOfNbFaces();
		auto itNbFace = itBeginCell[i].getBeginItPNbFaces();
		for (int iNbFace = 0; iNbFace < itBeginCell[i].getNumOfNbFaces();iNbFace++) {
			for (int iDim = 0; iDim < DIM ; iDim++){
				integralConvectFlux += (itDF[iNbFace])(iDim) * p1OldFaceConvectFlux[(itNbFace[iNbFace])->getId()].col(iDim);

			}
		}
#if FLAGPSUEDO 
		{//with psuedo time 
			NUM DPT = 1. / (NUM)PT_STEP_TIME;
			p1NewCellQs[i] = DPT*DT / (DPT + DT)*(ptP1OldNTimeCellQs[i] / DT + p1OldCellQs[i] / DPT - 1. / itBeginCell[i].getVolume() * (integralConvectFlux));
		}
#else
		{
			p1NewCellQs[i] = p1OldCellQs[i] - DT / (itBeginCell[i].getVolume()) *(integralConvectFlux);
		}
#endif // FLAGPSUEDO == 0
		
	}

#if FLAGVISCID 
	{//with viscid 
		updateViscid();
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
		for (int i = 0; i < pMesh->getNumOfCells(); i++) {
			VCTDIMU integralViscidFlux = VCTDIMU::Zero();
			auto itDF = itBeginCell[i].getBeginItDirectOfNbFaces();
			auto itNbFace = itBeginCell[i].getBeginItPNbFaces();
			for (int iNbFace = 0; iNbFace < itBeginCell[i].getNumOfNbFaces(); ) {
				for (int iDim = 0; iDim < DIM; iDim++) {
					integralViscidFlux += (itDF[iNbFace])(0)* p1OldFaceViscidFlux[(itNbFace[iNbFace])->getId()].col(iDim);
				}
			}
			p1NewCellQs[i] += DT / (itBeginCell[i].getVolume()) *(integralViscidFlux);

		}
	}
#endif // FLAGPSUEDO == 0

}
void RhoSolver::updateFaceFlux() {

	//solverRoe  solverSingle;//choose solver!!!!!!!!!!!
	//SolverNnd  solverSingle;//choose solver!!!!!!!!!!!
	//SolverAusm solverSingle;//choose solver!!!!!!!!!!!

	int x = pMesh->getNumOfFacesInfs();
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()], p1OldCellQs[itFace.getBeginItPNbCells()[1]->getId()]);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[1]->getId()], p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()]);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}
				

			}
			continue;
		}
		case 10: {//inlet face update 
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				VCTDIMU inletQ;
				inletQ << inletrho, inletrho*inletu, inletrho*inletv, inletE;
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto itFace = (itBeginFace)[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()], inletQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(inletQ, p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()]);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}

			}
			continue;
		}
		case 3: {//wall face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()];
				VCTDIMU outCellQ = inCellQ;

#if FLAGVISCID == 0
				// if wall smooth is avaliable
				VCTDIM directFace = itBeginFace[i].getDirect() / itBeginFace[i].getDirect().norm();
				VCTDIM uvwCellQ, uvwFaceQ;
				for (int iDim = 0; iDim < DIM; iDim++) {
					uvwCellQ[iDim] = inCellQ[iDim + 1];
				}
				uvwFaceQ = uvwCellQ - 2 * directFace*(directFace.transpose()*uvwCellQ);
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = uvwFaceQ[iDim];
				}
				
#else  //viscid wall :
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = -inCellQ[iDim + 1];
				}// if wall  uvw = 0

#endif
				
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(inCellQ, outCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(outCellQ, inCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}
			}
			continue;

		}
		case 7: {//symmetry face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()];
				VCTDIMU outCellQ = inCellQ;
				// if wall smooth is avaliable
				VCTDIM directFace = itBeginFace[i].getDirect() / itBeginFace[i].getDirect().norm();
				VCTDIM uvwCellQ, uvwFaceQ;
				for (int iDim = 0; iDim < DIM; iDim++) {
					uvwCellQ[iDim] = inCellQ[iDim + 1];
				}
				uvwFaceQ = uvwCellQ - 2 * directFace * (directFace.transpose() * uvwCellQ);
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = uvwFaceQ[iDim];
				}
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(inCellQ, outCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(outCellQ, inCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}
			}
			continue;

		}
		case 5: {//outlet face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()];
				VCTDIMU outCellQ = inCellQ;

				for (int iDim = 0; iDim < DIM; iDim++) {
					solverSingleTemp.set(outCellQ, inCellQ);
					p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
				}
			}
			continue;
		}
		}

	}//all face (interior & boundary) update
	
}
void RhoSolver::updateFaceFlux2ndOrder() {// 20230213 2nd order try
	updateGradFlux();
	//solverRoe  solverSingle;//choose solver!!!!!!!!!!!
	//SolverNnd  solverSingle;//choose solver!!!!!!!!!!!
	//SolverAusm solverSingle;//choose solver!!!!!!!!!!!

	int x = pMesh->getNumOfFacesInfs();
	for (int inf = 0; inf < pMesh->getNumOfFacesInfs(); inf++) {
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[0]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[0]->getId()].getCenter()), p1OldCellQs[itFace.getBeginItPNbCells()[1]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[1]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[1]->getId()].getCenter()));
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[1]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[1]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[1]->getId()].getCenter()), p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[0]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[0]->getId()].getCenter()));
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}				

			}
			continue;
		}
		case 10: {//inlet face update 
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				VCTDIMU inletQ;
				inletQ << inletrho, inletrho* inletu, inletrho* inletv, inletE;
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto itFace = (itBeginFace)[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[0]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[0]->getId()].getCenter()), inletQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(inletQ, p1OldCellQs[itFace.getBeginItPNbCells()[0]->getId()] + p1NewCellGradFlux[itFace.getBeginItPNbCells()[0]->getId()] * (itFace.getCenter() - itBeginCell[itFace.getBeginItPNbCells()[0]->getId()].getCenter()));
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}				

			}
			continue;
		}
		case 3: {//wall face update 
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()] + p1NewCellGradFlux[pCellBoundary->getId()] * (itBeginFace[i].getCenter() - itBeginCell[pCellBoundary->getId()].getCenter());
				VCTDIMU outCellQ = inCellQ;
#if FLAGVISCID == 0
				// if wall smooth is avaliable
				VCTDIM directFace = itBeginFace[i].getDirect() / itBeginFace[i].getDirect().norm();
				VCTDIM uvwCellQ, uvwFaceQ;
				for (int iDim = 0; iDim < DIM; iDim++) {
					uvwCellQ[iDim] = inCellQ[iDim + 1];
				}
				uvwFaceQ = uvwCellQ - 2 * directFace * (directFace.transpose() * uvwCellQ);
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = uvwFaceQ[iDim];
				}
#else  //viscid wall :
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = -inCellQ[iDim + 1];
				}// if wall  uvw = 0
#endif
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(inCellQ, outCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(outCellQ, inCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}
			}
			continue;
		}
		case 7: {//symmetry face update!!!!
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()];
				VCTDIMU outCellQ = inCellQ + p1NewCellGradFlux[pCellBoundary->getId()] * (itBeginFace[i].getCenter() - itBeginCell[pCellBoundary->getId()].getCenter());
				// if wall smooth is avaliable
				VCTDIM directFace = itBeginFace[i].getDirect() / itBeginFace[i].getDirect().norm();
				VCTDIM uvwCellQ, uvwFaceQ;
				for (int iDim = 0; iDim < DIM; iDim++) {
					uvwCellQ[iDim] = inCellQ[iDim + 1];
				}
				uvwFaceQ = uvwCellQ - 2 * directFace * (directFace.transpose() * uvwCellQ);
				for (int iDim = 0; iDim < DIM; iDim++) {
					outCellQ[iDim + 1] = uvwFaceQ[iDim];
				}
				auto itFace = (pMesh->getBeginItFacesList())[i];
				for (int iDim = 0; iDim < DIM; iDim++) {
					if (itFace.getFlagLeftRight()[iDim]) {
						solverSingleTemp.set(inCellQ, outCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
					else {
						solverSingleTemp.set(outCellQ, inCellQ);
						p1OldFaceConvectFlux[i].col(iDim) = solverSingleTemp.solverAll(iDim);
					}
				}
			}
			continue;
		}
		case 5: {//outlet face update
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				RHOSOLVER solverSingleTemp;//choose solver!!!!!!!!!!!
				auto pCellBoundary = *(itBeginFace[i]).getBeginItPNbCells();
				VCTDIMU inCellQ = p1OldCellQs[pCellBoundary->getId()] ; 
				VCTDIMU outCellQ = inCellQ + p1NewCellGradFlux[pCellBoundary->getId()] * (itBeginFace[i].getCenter() - itBeginCell[pCellBoundary->getId()].getCenter());

				solverSingleTemp.set(p1OldCellQs[pCellBoundary->getId()], outCellQ);
				p1OldFaceConvectFluxF[i] = solverSingleTemp.solverx();
				p1OldFaceConvectFluxG[i] = solverSingleTemp.solvery();
			}
			continue;
		}
		}

	}//all face (interior & boundary) update

}

void RhoSolver::updateViscid() {
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		p1OldFaceQs[i] = itBeginFace[i].getEta0()*p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + (1 - itBeginFace[i].getEta0())*p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[1]->getId()];
	}//simple interpolation   
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		p1OldFaceQs[i] = p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfCells(); i++) {
		auto itDF = itBeginCell[i].getBeginItDirectOfNbFaces();
		for (int iNbFace = 0; iNbFace < itBeginCell[i].getNumOfNbFaces();iNbFace++ ) {
			for (int iDim = 0; iDim < DIM; iDim++)	{
				pCellGradPUVT[i]( 0, iDim) += getP(p1OldFaceQs[i])*(itDF[iNbFace])(iDim);
				for (int iDim2 = 0; iDim2 < DIM; iDim2++) {
					pCellGradPUVT[i]( iDim2+1, iDim) += (p1OldFaceQs[i](iDim2) / p1OldFaceQs[i](0))*(itDF[iNbFace])(iDim);
				}
				pCellGradPUVT[i]( DIMU - 1, iDim) += getT(p1OldFaceQs[i])*(itDF[iNbFace])(iDim);
			}
		}
		pCellGradPUVT[i] /= itBeginCell[i].getVolume();
	}
#pragma omp parallel for num_threads(NUM_CPU_THREADS)     
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		pFaceGradPUVT[i] = itBeginFace[i].getEta0()*pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + (1 - itBeginFace[i].getEta0())*pCellGradPUVT[itBeginFace[i].getBeginItPNbCells()[1]->getId()];
	}//simple interpolation on inter faces
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		pFaceGradPUVT[i] = pFaceGradPUVT[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary faces
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfFaces(); i++) {
		
		//MTDIM_DIM tau;
		p1OldFaceViscidFlux[i] = MTDIMU_DIM::Zero();
		for (int iDim1 = 0; iDim1 < DIM; iDim1++){
			for (int iDim2 = 0; iDim2 < DIM; iDim2++) {
				p1OldFaceViscidFlux[i](iDim1+1, iDim2) = pFaceGradPUVT[i](iDim2 + 1, iDim1 + 1)+ pFaceGradPUVT[i](iDim1 + 1, iDim2 + 1);
			}
		}
		p1OldFaceViscidFlux[i] *= VISCIDMU;
		NUM bLambda = 0;
		for (int iDim = 0; iDim < DIM; iDim++) {
			bLambda += VISCIDLAMBDA*pFaceGradPUVT[i](iDim + 1, iDim + 1);
		}
		for (int iDim = 0; iDim < DIM; iDim++) {
			p1OldFaceViscidFlux[i](iDim+1, iDim) += bLambda;
		}
		
		for (int iDimFlux = 0; iDimFlux < DIM; iDimFlux++) { // loop of every flux
			for (int iDim2 = 0; iDim2 < DIM; iDim2++) {// loop of tau in vis flux
				p1OldFaceViscidFlux[i](DIMU - 1, iDimFlux) += p1OldFaceQs[i](iDim2 + 1) / p1OldFaceQs[i](0)*p1OldFaceViscidFlux[i](iDim2 + 1, iDimFlux);
			}
			p1OldFaceViscidFlux[i](DIMU - 1, iDimFlux) += TEMPK*pFaceGradPUVT[i](DIMU - 1);
		}
	}//just interpolation 

}
void RhoSolver::updateGradFlux() {

#if GRAD_INTERVAL == 1 // GRAD_INTERVAL == 1
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		p1OldFaceQs[i] = itBeginFace[i].getEta0() * p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + (1 - itBeginFace[i].getEta0()) * p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[1]->getId()];
	}//simple interpolation
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		p1OldFaceQs[i] = p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto itDF = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		MTDIMU_DIM tempMT = MTDIMU_DIM::Zero();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
			int idFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			for (int d = 0; d < DIM; d++) {
				tempMT.col(d) += p1OldFaceQs[idFace] * (itDF[iNbFace])(d);
			}
		}
		p1NewCellGradFlux[iCell] = tempMT / itBeginCell[iCell].getVolume();
	}


#else
	VCTDIMU* p1OldFaceQMiddle = new VCTDIMU[pMesh->getNumOfFaces()];
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int i = 0; i < pMesh->getNumOfIntFaces(); i++) {
		p1OldFaceQMiddle[i] = 0.5 * (p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()] + p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[1]->getId()]);
	}//simple interpolation   
#pragma omp parallel for num_threads(NUM_CPU_THREADS)
	for (int i = pMesh->getNumOfIntFaces() - 1; i < pMesh->getNumOfFaces(); i++) {
		p1OldFaceQMiddle[i] = p1OldCellQs[itBeginFace[i].getBeginItPNbCells()[0]->getId()];
	}//simple interpolation on doundary
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		auto itDF = itBeginCell[iCell].getBeginItDirectOfNbFaces();
		MTDIMU_DIM tempMT = MTDIMU_DIM::Zero();
		for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
			int idFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
			for (int d = 0; d < DIM; d++) {
				tempMT.col(d) += p1OldFaceQMiddle[idFace] * (itDF[iNbFace])(d);
			}
		}
		p1NewCellGradFlux[iCell] = tempMT / itBeginCell[iCell].getVolume();
	}
	// below is the interation of gradient calculation :
	for (int iGradloop = 0; iGradloop < GRAD_INTERVAL; iGradloop++) {
		// below is the interation of gradient calculation :
		for (int iFaceInt = 0; iFaceInt < pMesh->getNumOfIntFaces(); iFaceInt++) {
			p1OldFaceQs[iFaceInt] = p1OldFaceQMiddle[iFaceInt] + 0.5 * (p1NewCellGradFlux[itBeginFace[iFaceInt].getBeginItPNbCells()[0]->getId()] + p1NewCellGradFlux[itBeginFace[iFaceInt].getBeginItPNbCells()[1]->getId()]) * itBeginFace[iFaceInt].getCenterMiddleDef();
		}
#pragma omp parallel for num_threads(NUM_CPU_THREADS) 
		for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
			auto itDF = itBeginCell[iCell].getBeginItDirectOfNbFaces();
			MTDIMU_DIM tempMT = MTDIMU_DIM::Zero();
			for (int iNbFace = 0; iNbFace < itBeginCell[iCell].getNumOfNbFaces(); iNbFace++) {
				int idFace = itBeginCell[iCell].getBeginItPNbFaces()[iNbFace]->getId();
				for (int d = 0; d < DIM; d++) {
					tempMT.col(d) += p1OldFaceQs[idFace] * (itDF[iNbFace])(d);
				}
			}
			p1NewCellGradFlux[iCell] = tempMT / itBeginCell[iCell].getVolume();
		}


	}

#endif 


}

VCTDIMU* RhoSolver::getOldValue() {
	return pAllData->getP1OldCellQs();
}
VCTDIMU * RhoSolver::getNewValue() {
	return pAllData->getP1NewCellQs();
}
VCTDIMU* RhoSolver::getOldNTimeValue() {
	return pAllData->getPtP1OldNTimeCellQs();
}
void RhoSolver::updateNewToOld() {
	for (int iCell = 0; iCell < pMesh->getNumOfCells(); iCell++) {
		p1OldCellQs[iCell] = p1NewCellQs[iCell];
	}
}


