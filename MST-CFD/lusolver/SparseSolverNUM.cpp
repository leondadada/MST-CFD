#pragma once
#include "SparseSolverNUM.h"


SparseSolverNUM::SparseSolverNUM() {

}
SparseSolverNUM::SparseSolverNUM(int _dim){
	dim = _dim;
	this->giveNew();
	for (int i = 0; i < dim; i++) {
		X[i] = 0;
		D[i] = 0;
	}
}
SparseSolverNUM::SparseSolverNUM(int _dim, int _flagX, VCTDIMU* _pOldX) {
	dim = _dim;
	giveNew();

	for (int i = 0; i < dim; i++) {
		X[i] = _pOldX[i](_flagX);
		D[i] = 0;
	}
}
SparseSolverNUM::SparseSolverNUM(int _dim, NUM* _pOldX){
	dim = _dim;
	giveNew();

	for (int i = 0; i < dim; i++) {
		X[i] = _pOldX[i];
		D[i] = 0;
	}
}

void SparseSolverNUM::reload(int _dim, int _flagX, VCTDIMU* _pOldX) {
	dim = _dim;
	deleteAll();
	giveNew();

	for (int i = 0; i < dim; i++) {
		X[i] = _pOldX[i](_flagX);
	}
}
void SparseSolverNUM::reload(int _dim, NUM* _pOldX) {
	dim = _dim;
	deleteAll();
	giveNew();

	for (int i = 0; i < dim; i++) {
		X[i] = _pOldX[i];
	}
}
void SparseSolverNUM::reload(int _dim) {
	dim = _dim;
	deleteAll();
	giveNew();
}
void SparseSolverNUM::reloadUVP(int _dim,VCTDIMU* _PUVT) {
	dim = _dim;
	deleteAll();
	giveNew();

	for (int i = 0; i < round(dim / 3); i++) {
		X[3 * i] = _PUVT[i][1];
		X[3 * i + 1] = _PUVT[i][2];
		X[3 * i + 2] = _PUVT[i][0];
	}
}
void SparseSolverNUM::reloadRHS(int _dim, int _flagX, VCTDIMU* _pOldX) {
	dim = _dim;

	delete[] RHS;
	delete[] RHS1;
	delete[] X;
	delete[] X1;
	delete[] RHSb;
	delete[] RHSLDUx;
	delete[] RHSUx;

	RHS = new NUM[dim];
	RHSb = new NUM[dim];
	RHSLDUx = new NUM[dim];
	RHSUx = new NUM[dim];
	RHS1 = new NUM[dim];
	X = new NUM[dim];
	X1 = new NUM[dim];

	for (int i = 0; i < dim; i++) {
		X[i] = _pOldX[i](_flagX);
	}
}
SparseSolverNUM::~SparseSolverNUM(){
	this->deleteAll();	
}
int SparseSolverNUM::getDim(){
	return dim;
}
void SparseSolverNUM::addD(NUM _mt, int i){
	D[i] += _mt;
}
void SparseSolverNUM::setD(NUM _mt, int i){
	D[i] = _mt;
}
NUM SparseSolverNUM::getD(int _pst) {
	return this->D[_pst];
}
void SparseSolverNUM::setL(NUM _mt,  int _row, int _column) {
	columnL[_column].push_back(_mt);
	idColumnDL[_column].push_back(_row);
}
void SparseSolverNUM::setU(NUM _mt, int _row, int _column) {
	columnU[_column].push_back(_mt);
	idColumnDU[_column].push_back(_row);
}
void SparseSolverNUM::setELE(NUM _mt, int _row, int _column) {
	if (_row == _column) {
		setD(_mt,_row);
	}
	else {
		if (_row > _column) {
			setL(_mt, _row,_column);
		}
		else {
			setU(_mt, _row, _column);
		}
	}
}
void SparseSolverNUM::addELE(NUM _mt, int _row, int _column) {
	if (_row == _column) {
		addD(_mt, _row);
	}
	else {
		if (_row > _column) {
			setL(_mt, _row, _column);
		}
		else {
			setU(_mt, _row, _column);
		}
	}
}
void SparseSolverNUM::setRHSb(NUM _vt, int _column) {
	RHSb[_column] =(_vt);
}
void SparseSolverNUM::solveILUSGS(){
	//cout << "interval: ";
	//for (int i = 0; i < dim; i++){
	//	cout << i << ": " << D[i] << endl;
	//}
	//system("pause");


	for (int interval = 0; interval < LU_INTERVAL;interval++) {
		 //cout<< interval ;
		for (int i = 0; i < dim; i++) {
			RHSUx[i] = 0.;
			RHSLDUx[i] = 0.;
		}
		for (int i = 0; i < dim; i++) {
			if (idColumnDU[i].size() == 0) { continue; }
			for (int iCU = 0; iCU < idColumnDU[i].size(); iCU++) {
				RHSUx[idColumnDU[i][iCU]] += columnU[i][iCU] * X[i];
			}
		}
		for (int i = 0; i < dim; i++) {
			RHSUx[i] /= D[i];
		}
		for (int i = 0; i < dim; i++) {
			if (idColumnDL[i].size() == 0) { continue; }
			for (int iCL = 0; iCL < idColumnDL[i].size(); iCL++) {
				RHSLDUx[idColumnDL[i][iCL]] += columnL[i][iCL] * RHSUx[i];
			}
		}
		for (int i = 0; i < dim; i++) {
			RHS[i] = RHSb[i] + RHSLDUx[i];
		}
		

		for (int i = 0; i < dim; i++) {			
			if (idColumnDL[i].size() == 0) { continue; }
			for (int iCL = 0; iCL < idColumnDL[i].size();iCL++) {
				RHS[(idColumnDL[i])[iCL]]  -= (columnL[i][iCL])*(1./D[i]) *RHS[i];
			}
		}
		for (int i = 0; i < dim; i++){
			X1[i] = (1. / D[i])*RHS[i];//!!!
			RHS1[i] = D[i] * X1[i];
		}
		for (int i = dim - 1; i != -1; i--) {
			if (idColumnDU[i].size() == 0) { continue; }
			for (int iCU = 0; iCU < idColumnDU[i].size(); iCU++) {
				RHS1[idColumnDU[i][iCU]] -= (columnU[i][iCU])*(1. / D[i]) * RHS1[i];
			}
		}
		NUM  residualNew = 0;
		for (int i = 0; i < dim; i++) {
			if (residualNew < abs(X[i] - (1. / D[i]) * RHS1[i]) / X[i]) {
				residualNew = abs(X[i] - (1. / D[i]) * RHS1[i]) / X[i];
			}
			X[i] = (1. / D[i]) * RHS1[i];
			
			// cout << "RHS["<<i<<"]"<<RHSb[i] << endl;
			//cout << "X["<<i<<"]"<<X[i] << endl;
		}
		cout <<"LU-residual-IILU "<<interval<<" : "<< residualNew << endl;
		if (residualNew  > EOR* EOR && residualNew < EOR2) {
			break;
		}
	}
	
	//cout << endl;

}

void SparseSolverNUM::solveGS(){
	for (int interval = 0; interval < LU_INTERVAL; interval++) {
		for (int i = 0; i < dim; i++) {
			RHSUx[i] = 0.;
		}
		for (int i = 0; i < dim; i++) {
			if (idColumnDU[i].size() == 0) { continue; }
			for (int iCU = 0; iCU < idColumnDU[i].size(); iCU++) {
				RHSUx[idColumnDU[i][iCU]] += columnU[i][iCU] * X[i];
			}
		}
		for (int i = 0; i < dim; i++) {
			RHS[i] = RHSb[i] - RHSUx[i];
		}


		for (int i = 0; i < dim; i++) {
			if (idColumnDL[i].size() == 0) { continue; }
			for (int iCL = 0; iCL < idColumnDL[i].size(); iCL++) {
				RHS[(idColumnDL[i])[iCL]] -= (columnL[i][iCL])*(1. / D[i]) *RHS[i];
			}
		}
		NUM  residualNew = 0;
		for (int i = 0; i < dim; i++) {
			if (residualNew < abs(X[i] - (1. / D[i]) * RHS[i]) / X[i]) {
				residualNew = abs(X[i] - (1. / D[i]) * RHS[i]) / X[i];
			}
			X[i] = (1. / D[i]) *RHS[i] ;

			// cout << "RHS["<<i<<"]"<<RHSb[i] << endl;
			//cout << "X["<<i<<"]"<<X[i] << endl;
		}
		cout << "LU-residual-GS " << interval << " : " << residualNew << endl;
		if (residualNew  > EOR* EOR && residualNew < 1e-7) {
			break;
		}
	}

}

NUM* SparseSolverNUM::getPNewX(){
	return X;
}
void SparseSolverNUM::giveNew(){
	D = new NUM[dim];
	RHS = new NUM[dim];
	RHSb = new NUM[dim];
	RHSLDUx = new NUM[dim];
	RHSUx = new NUM[dim];
	RHS1 = new NUM[dim];
	X = new NUM[dim];
	X1 = new NUM[dim];

	columnL = new vector<NUM>[dim];
	columnU = new vector<NUM>[dim];

	idColumnDL = new vector<int>[dim];
	idColumnDU = new vector<int>[dim];
}
void SparseSolverNUM::deleteAll(){
	delete[] D;
	delete[] RHS;
	delete[] RHS1;
	delete[] X;
	delete[] X1;
	delete[] RHSb;
	delete[] RHSLDUx;
	delete[] RHSUx;
	delete[] idColumnDL;
	delete[] idColumnDU;
	delete[] columnL;
	delete[] columnU;
}




