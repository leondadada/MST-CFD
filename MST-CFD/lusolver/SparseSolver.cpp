#pragma once
#include "SparseSolver.h"


template <typename TMP_MT, typename TMP_VCT>
SparseSolver<TMP_MT, TMP_VCT>::SparseSolver(int _dim, TMP_VCT* _pOldX){
	dim = _dim;
	pOldX = _pOldX;

	D = new TMP_MT[dim];
	RHS = new TMP_VCT[dim];
	RHSb = new TMP_VCT[dim];
	RHSLDUx = new TMP_VCT[dim];
	RHSUx = new TMP_VCT[dim];
	RHS1 = new TMP_VCT[dim];
	X = new TMP_VCT[dim];
	X1 = new TMP_VCT[dim];

	columnL = new vector<TMP_MT> [dim];
	columnU = new vector<TMP_MT> [dim];

	idColumnDL = new vector<int>[dim];
	idColumnDU = new vector<int>[dim];
}
template <typename TMP_MT, typename TMP_VCT>
SparseSolver<TMP_MT, TMP_VCT>::~SparseSolver(){

	delete[] D, RHS, RHS1, X, X1, columnL, columnU;
	delete[] idColumnDL, idColumnDU, RHSb, RHSLDUx, RHSUx;
}
template <typename TMP_MT, typename TMP_VCT>
int SparseSolver<TMP_MT, TMP_VCT>::getDim(){
	return dim;
}
template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::setD(TMP_MT _mt, int i){
	D[i] = _mt;
}
template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::setL(TMP_MT _mt, int _row, int _column) {
	columnL[_column].push_back(_mt);
	idColumnDL[_column].push_back(_row);
}
template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::setU(TMP_MT _mt, int _row, int _column) {
	columnU[_column].push_back(_mt);
	idColumnDU[_column].push_back(_row);
}
template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::setRHSb(TMP_VCT _vt, int _column) {
	RHSb[_column] =(_vt);
}
template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::solveILU(){
	for (int i = 0; i < dim; i++) {
		X[i] = pOldX[i];
	}
	for (int interval = 0; interval < LU_INTERVAL; interval++) {
		cout << "interval: " << interval << endl;
		for (int i = 0; i < dim; i++) {
			RHSUx[i] = TMP_VCT::Zero();
			RHSLDUx[i] = TMP_VCT::Zero();
		}
		for (int i = 0; i < dim; i++) {
			if (idColumnDU[i].size() == 0) { continue; }
			for (int iCU = 0; iCU < idColumnDU[i].size(); iCU++) {
				RHSUx[idColumnDU[i][iCU]] += columnU[i][iCU] * X[i];
			}
		}
		for (int i = 0; i < dim; i++) {
			RHSUx[i] = D[i].inverse()*RHSUx[i];
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
			for (int iCL = 0; iCL < idColumnDL[i].size(); iCL++) {
				RHS[(idColumnDL[i])[iCL]] -= (columnL[i][iCL])*(D[i].inverse()) *RHS[i];
			}
		}
		for (int i = 0; i < dim; i++) {
			X1[i] = (D[i].inverse())*RHS[i];//!!!
			RHS1[i] = D[i] * X1[i];
		}
		for (int i = dim - 1; i != -1; i--) {
			if (idColumnDU[i].size() == 0) { continue; }
			for (int iCU = 0; iCU < idColumnDU[i].size(); iCU++) {
				RHS1[idColumnDU[i][iCU]] -= (columnU[i][iCU])*(D[i].inverse()) * RHS1[i];
			}

		}
		for (int i = 0; i < dim; i++) {
			X[i] = ( D[i].inverse()) * RHS1[i];
		}
	}
}

template <typename TMP_MT, typename TMP_VCT>
void SparseSolver<TMP_MT, TMP_VCT>::solveGS(){


}
template <typename TMP_MT, typename TMP_VCT>
TMP_VCT* SparseSolver<TMP_MT, TMP_VCT>::getItBeginX(){
	return X;
}





