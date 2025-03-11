#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "../include/CONST.h"
#include "../work/FUNCTION.h"


template <typename TMP_MT, typename TMP_VCT>class SparseSolver
{
public:

	SparseSolver(int, TMP_VCT*);
	~SparseSolver();
	int getDim();
	void setD(TMP_MT, int );
	void setL(TMP_MT, int, int);
	void setU(TMP_MT, int, int);
	void setRHSb(TMP_VCT, int);
	void solveILU();
	void solveGS();
	TMP_VCT* getItBeginX();
private:
	int dim;
	TMP_VCT* pOldX;
	TMP_MT* D;
	TMP_VCT* RHSb;
	TMP_VCT* X;	
	vector<TMP_MT> * columnL;
	vector<TMP_MT> * columnU;
	vector<int>* idColumnDL;
	vector<int>* idColumnDU;

	TMP_VCT* RHS;
	TMP_VCT* X1;
	TMP_VCT* RHSUx;
	TMP_VCT* RHSLDUx;
	TMP_VCT* RHS1;



	//vector<TMP_MT> D;
	//vector<TMP_VCT> RHS;
	//vector<TMP_VCT> RHSb;
	//vector<TMP_VCT> RHSLDUx;
	//vector<TMP_VCT> RHS1;
	//vector<TMP_VCT> X;
	//vector<TMP_VCT> X1;
	//vector<vector<TMP_MT>>  columnL;
	//vector<vector<TMP_MT>> columnU;
	//vector<vector<int>> idColumnDL;
	//vector<vector<int>> idColumnDU;
};

