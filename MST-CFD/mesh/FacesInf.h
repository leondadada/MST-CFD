#pragma once
#include "../include/CONST.h"
#include "./Face.h"
#include "./Node.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>
#include "math.h"
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;
class FacesInf
{
public :
	void setNumber(int);
	void setType(int);
	void setStart(int);
	void setEnd(int);
	void setName(string);
	int getType();
	int getStart();
	int getEnd();
	int getNumber();
	string getName();
private:
	string name;
	int number;
	int type;
	int start;
	int end;
	
};
