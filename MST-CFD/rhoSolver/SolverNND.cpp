#include "SolverNnd.h"
#include "../include/CONST.h"


void SolverNnd::set(Cell* _pLeftCell, Cell* _pRightCell){
	this->leftCell = *_pLeftCell;
	this->rightCell = *_pRightCell;
	

}


VCTDIMU  SolverNnd::solverx() {
	
	return pMatrix;
}
VCTDIMU SolverNnd::solvery() {
	return pMatrix;

}
