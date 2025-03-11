#pragma once
#include <Eigen/Eigen>
#include <iostream>
#include "math.h"
#include "../include/CONST.h"
#include "./Face.h"
#include "./Node.h"

using namespace Eigen;

class Node;
class Face;
class Cell
{
public:

	void setCenter();
	void setVolume();
	void setFlagBoundary();
	void addNbFaces(Face*);
	void addNbCells(Cell*);
	void addNbNodes(Node*);
	void addDirectOfNbFaces(int,VCTDIM);
	void addNbCellsL(Cell*);//20230217-L
	void addNbCellsU(Cell*);//20230217-U
	void creatDirectOfNbFaces();
	void creatEta0C();
	void addEta0C(int,NUM);
	void setId(int);
	int getId();
	NUM getVolume();
	bool getFalgBoundary();

	Matrix<NUM, DIM, 1> getCenter();
	std::vector<Cell*>::iterator getBeginItPNbCells();
	std::vector<Face*>::iterator getBeginItPNbFaces();
	std::vector<Node*>::iterator getBeginItPNbNodes();
	VCTDIM* getBeginItDirectOfNbFaces();
	NUM* getBeginItEta0C();
	std::vector<Cell*>::iterator getBeginItPNbCellsL();
	std::vector<Cell*>::iterator getBeginItPNbCellsU();
	int getNumOfNbFaces();
	int getNumOfNbCells();
	int getNumOfNbNodes();
private:
	NUM volume;
	int id;
	bool flagBoundary = 0;
	std::vector<Cell*> pNbCellL;//20230217-nbCellPoint-L
	std::vector<Cell*> pNbCellU;//20230217-nbCellPoint-U

	Matrix<NUM, DIM, 1> center;
	std::vector<Cell*> pNbCells;// 
	std::vector<Node*> pNbNodes;// 
	std::vector<Face*> pNbFaces;//nabougher faces' pointers
	VCTDIM* directOfNbFaces;
	NUM* eta0C;
};

