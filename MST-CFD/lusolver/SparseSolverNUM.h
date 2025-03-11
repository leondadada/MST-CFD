#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "../include/CONST.h"
#include "../work/FUNCTION.h"


class SparseSolverNUM
{
public:

	SparseSolverNUM();

	SparseSolverNUM(int);
	SparseSolverNUM(int,int ,VCTDIMU*);
	SparseSolverNUM(int, NUM*);
	~SparseSolverNUM();
	int getDim();
	void addD(NUM, int);
	void setD(NUM, int);
	NUM getD(int _pst);
	void setL(NUM, int, int);
	void setU(NUM, int, int);
	void setELE(NUM, int, int);
	void addELE(NUM , int , int );
	void setRHSb(NUM, int);
	void solveILUSGS();
	void solveGS();
	void reload(int , int , VCTDIMU* );
	void reload(int , NUM* );
	void reload(int _dim);
	void reloadUVP(int _dim, VCTDIMU * _PUVT);
	void reloadRHS(int _dim, int _flagX, VCTDIMU* _pOldX);
	NUM * getPNewX();
private:
	int dim;
	NUM* D;
	NUM* RHSb;
	NUM* X;	
	vector<NUM> * columnL;
	vector<NUM> * columnU;
	vector<int>* idColumnDL;
	vector<int>* idColumnDU;
	void giveNew();
	void deleteAll();
	NUM* RHS;
	NUM* X1;
	NUM* RHSLDUx;
	NUM* RHSUx;
	NUM* RHS1;


};

