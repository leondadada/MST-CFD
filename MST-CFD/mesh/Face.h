#pragma once
#include <Eigen/Eigen>
#include "../include/CONST.h"
#include "./Cell.h"
#include "./Node.h"

using namespace Eigen;
class Cell;
class Node;
class Face
{
public:
	void setCenter();
	void setDirect(VCTDIM);
	void setId(int);
	void setType(int);
	void setFlagLeftRight(int, bool);

	void addNbCell(Cell*);
	void updateEta0();
	void updateCenterMiddleDef();
	void updateBondEta0();
	void tranNbCells();
	void transDirectAndCells();
	void addNbNode(Node*);
	void interpolationUOfCell();
	NUM getEta0();
	VCTDIM getCenter();
	VCTDIM getCenterMiddleDef();
	VCTDIM getDirect();
	std::vector<Cell*>::iterator getBeginItPNbCells();
	std::vector<Cell*>::iterator getEndItPNbCells();
	std::vector<Node*>::iterator getBeginItPNbNodes();
	std::vector<Node*>::iterator getEndItPNbNodes();
	int getNumOfpNbNodes();
	int getNumOfpNbCells();
	int getDirectAndCells();
	bool* getFlagLeftRight();
	int getType();
	int getId();
	NUM getArea();

private:
	Matrix<NUM, DIM, 1> center;//face's center (x,y,z)
	Matrix<NUM, DIM, 1> centerMiddleDef;//face's middle center (x,y,z)
	Matrix<NUM, DIM, 1> directArea;//with area multiply inside
	NUM eta0 = 1;//distance f-cell0/(f-cell0+f-cell1);
	NUM area;//
	int type;
	int id;
	int directAndCells = 1;//true means first cell in pCellVector is 
	std::vector<Cell*> pNbCells;//nabougher cells' pointers 
	std::vector<Node*> pNbNodes;//nabougher faces' pointers
	bool flagLeftRight[DIM];//if 0Cell is left?
};
