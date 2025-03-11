#include "./Face.h"
#include "./Cell.h"
#include "./Node.h"
#include <iostream>
#include <Eigen/Eigen>
using namespace Eigen;

void Face::setCenter(){
	VCTDIM sum = VCTDIM::Zero();
	for (int iNbNode = 0; iNbNode < pNbNodes.size(); iNbNode++) {
		sum += pNbNodes[iNbNode]->getCenter();
	}
	this->center = sum / (NUM)pNbNodes.size();
#if  DIM == 2
	{		
		VCTDIM f = pNbNodes[0]->getCenter() - pNbNodes[1]->getCenter();
		area = f.norm();
		directArea << -f(1) , f(0) ;
	}
#else 
	{//dim ==3
		switch (pNbNodes.size()){
		case 3: {
			VCTDIM f1 = pNbNodes[1]->getCenter() - pNbNodes[0]->getCenter();
			VCTDIM f2 = pNbNodes[2]->getCenter() - pNbNodes[0]->getCenter();
			directArea = 0.5 * f1.cross(f2);
			area = directArea.norm();
			break;
		}
		case 4: {// 4 line on face
			VCTDIM f1 = pNbNodes[1]->getCenter() - pNbNodes[0]->getCenter();
			VCTDIM f2 = pNbNodes[2]->getCenter() - pNbNodes[0]->getCenter();
			directArea = f1.cross(f2);
			area = directArea.norm();
			break;
		}

		}
			
		
	}
#endif

}
void Face::setDirect(VCTDIM _direct){
	this->directArea = _direct;
}

void Face::setId(int _id) {
	this->id = _id;
}
void Face::setType(int _type){
	type = _type;
}
int Face::getId() {
	return this->id;
}
void Face::addNbCell(Cell* _pNbCell){
	this->pNbCells.push_back(_pNbCell);
}

void Face::updateEta0(){
	if (this->pNbCells.size() == 2) {
		this->eta0 = (pNbCells[1]->getCenter() - this->center).norm() / ((pNbCells[0]->getCenter() - this->center).norm() + (pNbCells[1]->getCenter() - this->center).norm());
	}
	else {
		this->eta0 = 1;
	}
}
void Face::updateCenterMiddleDef(){
	centerMiddleDef = center - 0.5*(pNbCells[0]->getCenter() + pNbCells[1]->getCenter());
}
void Face::tranNbCells(){
	Cell* temp = this->pNbCells[0];
	if(pNbCells.size() == 2){
		this->pNbCells[0] = this->pNbCells[1];
		this->pNbCells[1] = temp;
	}
	else {
		addNbCell(temp);
		this->pNbCells[0] = NULL;
	}
}

void Face::transDirectAndCells(){
	this->directAndCells = -1;
}
void Face::setFlagLeftRight(int _i, bool _flag){
	flagLeftRight[_i] = _flag;
}
void Face::addNbNode(Node* _pNbNode){
	this->pNbNodes.push_back(_pNbNode);
}
NUM Face::getEta0(){
	return this->eta0;
}
VCTDIM Face::getCenter(){
	return this->center;
}
VCTDIM Face::getCenterMiddleDef(){
	return this->centerMiddleDef;
}
VCTDIM Face::getDirect() {
	return this->directArea;
}
std::vector<Cell*>::iterator Face::getBeginItPNbCells() {
	return this->pNbCells.begin();
}
std::vector<Cell*>::iterator Face::getEndItPNbCells(){
	return this->pNbCells.end();
}
std::vector<Node*>::iterator Face::getBeginItPNbNodes() {
	return this->pNbNodes.begin();
}
std::vector<Node*>::iterator Face::getEndItPNbNodes() {
	return this->pNbNodes.end();
}
int Face::getNumOfpNbNodes(){
	return this->pNbNodes.size();
}
int Face::getNumOfpNbCells(){
	return this->pNbCells.size();
}
int Face::getDirectAndCells(){
	return this->directAndCells;
}
bool* Face::getFlagLeftRight(){
	return flagLeftRight;
}
int Face::getType(){
	return type;
}
NUM Face::getArea(){
	return this->area;
}

