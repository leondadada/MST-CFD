#include "./Node.h"
#include "./Face.h"
#include "./Cell.h"

#include <Eigen/Eigen>
using namespace Eigen;

Node::Node(){
}
Node::Node( VCTDIM _center){
	setCenter(_center);
}
void Node::addNbFace(Face* _pFace) {
	this->pNbFaces.push_back(_pFace);
}
VCTDIM Node::getCenter(){
	return this->center;
}
VCTDIMU Node::getU(){
	return this->U;
}
void Node::setId(int _id) {
	this->id = _id;
}
int Node::getId() {
	return this->id;
}
void Node::setCenter(VCTDIM _center){
	this->center = _center;
}
std::vector<Face*>::iterator Node::getBeginItPNbFaces(){
	return this->pNbFaces.begin();
}
std::vector<Face*>::iterator Node::getEndItPNbFaces(){
	return this->pNbFaces.end();
}
int Node::getNumOfpNbFaces(){
	return this->pNbFaces.size() ;
}
void Node::setPoint(int _xyz,NUM _x) {
	this->center(_xyz) = _x;
}





