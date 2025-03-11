#include "./Cell.h"

#include <Eigen/Eigen>
using namespace Eigen;

void Cell::setCenter() {
	Matrix<NUM, DIM, 1> sum = Matrix<NUM, DIM, 1>::Zero();
	auto it = this->pNbFaces.begin();
	while (it != this->pNbFaces.end()) {
		sum += (*it)->getCenter();
		it++;
	}
	this->center = sum / (NUM)pNbFaces.size();//center 
}
void Cell::setVolume() {

	if (DIM == 2) {
		if (this->pNbFaces.size() == 3) {

			NUM length[3] = { 0 };
			for (int i = 0; i < 3; i++) {
				length[i] = pNbFaces[i]->getDirect().norm();
			}
			NUM s = 0.5*(length[0] + length[1] + length[2]);
			this->volume = sqrt(s*(s - length[0])*(s - length[1])*(s - length[2]));
			//std::cout << (s - length[0]) << std::endl;
		}
		if (this->pNbFaces.size() == 4) {
			NUM middleLine = 0;

			VCTDIM vct1 = pNbFaces[0]->getCenter() - pNbFaces[1]->getCenter();
			VCTDIM vct2 = pNbFaces[2]->getCenter() - pNbFaces[3]->getCenter();
			if (abs(vct1[0] * vct2[1] - vct2[0] * vct1[1]) < EOR2) {

				// std::cout << this->id<< "\t!"<< std::endl;
				middleLine = 2 * (pNbFaces[0]->getCenter() - pNbFaces[1]->getCenter()).norm();
				NUM s = 0.5*(pNbFaces[0]->getDirect().norm() + pNbFaces[1]->getDirect().norm() + middleLine);
				this->volume = sqrt(s*(s - pNbFaces[0]->getDirect().norm())*(s - pNbFaces[1]->getDirect().norm())*(s - middleLine));
				NUM s2 = 0.5*(pNbFaces[2]->getDirect().norm() + pNbFaces[3]->getDirect().norm() + middleLine);
				this->volume += sqrt(s2*(s2 - pNbFaces[2]->getDirect().norm())*(s2 - pNbFaces[3]->getDirect().norm())*(s2 - middleLine));
			}
			else {
				middleLine = 2 * (pNbFaces[0]->getCenter() - pNbFaces[2]->getCenter()).norm();
				NUM s = 0.5*(pNbFaces[0]->getDirect().norm() + pNbFaces[2]->getDirect().norm() + middleLine);
				this->volume = sqrt(s*(s - pNbFaces[0]->getDirect().norm())*(s - pNbFaces[2]->getDirect().norm())*(s - middleLine));
				NUM s2 = 0.5*(pNbFaces[1]->getDirect().norm() + pNbFaces[3]->getDirect().norm() + middleLine);
				this->volume += sqrt(s2*(s2 - pNbFaces[1]->getDirect().norm())*(s2 - pNbFaces[3]->getDirect().norm())*(s2 - middleLine));
			}
		}

	}
	else {//dim == 3
		volume = 0;
		for (int iNbFace = 0; iNbFace < pNbFaces.size(); iNbFace++){
			VCTDIM unitDirectFace = pNbFaces[iNbFace]->getDirect()/ pNbFaces[iNbFace]->getDirect().norm();//pNbFaces[iNbFace]->getDirect() .Unit()
			VCTDIM tempVector = center-pNbFaces[iNbFace]->getCenter();
			volume += pNbFaces[iNbFace]->getArea()*tempVector.transpose()*unitDirectFace;
		}
	}
	
}

void Cell::setFlagBoundary(){
	flagBoundary = true;
}
bool Cell::getFalgBoundary(){
	return flagBoundary;
}
NUM Cell::getVolume() {
	return this->volume;
}
void Cell::addNbFaces(Face* _pNbFace) {
	this->pNbFaces.push_back(_pNbFace);
}
void Cell::addNbCells(Cell* _pNbCell ){
	this->pNbCells.push_back(_pNbCell);
}
void Cell::addNbNodes(Node * _pNbNode){
	this->pNbNodes.push_back(_pNbNode);
}
void Cell::addDirectOfNbFaces(int _position,VCTDIM _flag){
	directOfNbFaces[_position] = _flag;
}
void Cell::addNbCellsL(Cell* _pNbCell) {
	this->pNbCellL.push_back(_pNbCell);
}
void Cell::addNbCellsU(Cell* _pNbCell) {
	this->pNbCellU.push_back(_pNbCell);
}
void Cell::creatDirectOfNbFaces(){
	directOfNbFaces = new VCTDIM[pNbFaces.size()];
}
void Cell::creatEta0C(){
	eta0C = new NUM[pNbFaces.size()];
}
void Cell::addEta0C(int _position, NUM _et){
	eta0C[_position] = (_et);
}
void Cell::setId(int _id) {
	this->id = _id;
}
int Cell::getId() {
	return this->id;
}
Matrix<NUM, DIM, 1> Cell::getCenter() {
	return this->center;
}
std::vector<Cell*>::iterator Cell::getBeginItPNbCells() {
	return this->pNbCells.begin();
}
std::vector<Face*>::iterator Cell::getBeginItPNbFaces() {
	return this->pNbFaces.begin();
}
std::vector<Node*>::iterator Cell::getBeginItPNbNodes(){
	return this->pNbNodes.begin();
}
VCTDIM* Cell::getBeginItDirectOfNbFaces(){
	return this->directOfNbFaces;
}
std::vector<Cell*>::iterator Cell::getBeginItPNbCellsL() {
	return pNbCellL.begin();
}
std::vector<Cell*>::iterator Cell::getBeginItPNbCellsU() {
	return pNbCellU.begin();
}
NUM* Cell::getBeginItEta0C(){
	return eta0C;
}
int Cell::getNumOfNbFaces() {
	return this->pNbFaces.size();
}
int Cell::getNumOfNbCells() {
	return this->pNbCells.size();
}
int Cell::getNumOfNbNodes() {
	return this->pNbNodes.size();
}



