#include "./FacesInf.h"

void FacesInf::setNumber(int _num) {
	this->number = _num;
}
void FacesInf::setType(int _num) {
	this->type = _num;
}
void FacesInf::setStart(int _num) {
	this->start = _num;
}
void FacesInf::setEnd(int _num) {
	this->end = _num;
}
void FacesInf::setName(string _name){
	this->name = _name;
}
int FacesInf::getNumber() {
	return this->number;
}
string FacesInf::getName(){
	return this->name ;
}
int FacesInf::getType(){
	return this->type;
}
int FacesInf::getStart(){
	return this->start;
}
int FacesInf::getEnd(){
	return this->end;
}






