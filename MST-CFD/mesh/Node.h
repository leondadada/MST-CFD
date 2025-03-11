#pragma once
#include <Eigen/Eigen>
#include "../include/CONST.h"
#include "./Cell.h"
#include "./Face.h"

using namespace Eigen;
class Cell;
class Face;
class Node
{
public:
	Node();
	Node(Matrix<NUM,DIM,1>);
	void addNbFace(Face*);
	Matrix<NUM, DIM, 1> getCenter();
	Matrix<NUM, DIMU, 1> getU();
	void setPoint(int,NUM);
	void setId(int);
	int getId();
	void setCenter(VCTDIM);
	std::vector<Face*>::iterator getBeginItPNbFaces();
	std::vector<Face*>::iterator getEndItPNbFaces();
	int getNumOfpNbFaces();
private:
	int id;
	std::vector<Face*> pNbFaces;
	Matrix<NUM, DIM, 1> center;
	Matrix<NUM, DIMU, 1> U = Matrix<NUM, DIMU, 1>::Zero();//primitive  variables  <u,v,w,p,T>
	
};

