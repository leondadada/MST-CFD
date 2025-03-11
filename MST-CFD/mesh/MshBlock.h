#pragma once
#include "../include/CONST.h"
#include "../work/FUNCTION.h"
#include "./Cell.h"
#include "./Node.h"
#include "./Face.h"
#include "./FacesInf.h"
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include "math.h"

using namespace std;

class MshBlock 
{
public:
	MshBlock();
	void readMsh(string);
	Node* getBeginItNodesList();
	Cell* getBeginItCellsList();
	Face* getBeginItFacesList();
	Cell* getBeginItFakeCellsList();
	std::vector<FacesInf>::iterator getBeginItFacesInfList();
	std::vector<FacesInf>::iterator getEndItFacesInfList();

	int* getBeginItInnerIdList();
	int* getBeginItBond1IdList();
	int* getBeginItBond2IdList();



	int getNumOfNodes();
	int getNumOfCells();
	int getNumOfFaces();
	int getNumOfIntFaces();
	int getNumOfFacesInfs();
	int getNumOfInnerCells();
	int getNumOfBond1Cells();
	int getNumOfBond2Cells();
	
private:
	string mshInformation = "";
	string mshData;//all xyz relationship data 
	string stringLine;
	string stringLineData;
	string mshName;

	int dimMesh;
	int splitInformationData2D();
	int dataToMesh2D();
	int splitInformationData3D();
	int dataToMesh3D();
	int numOfNodes;
	int numOfCells;
	int numOfIntFaces;
	int numOfFaces;
	Node* nodesList;//nodes
	Cell* cellsList;//cells 
	Face* facesList;//interior faces
	std::vector<FacesInf> facesInfList;//boundary facesList,first will be int face list
	
	Cell* fakeCellsList;//fake cells 
	int innerIdListSize;
	int bond1IdListSize;
	int bond2IdListSize;
	std::vector<int> bond1IdListVector;
	std::vector<int> bond2IdListVector;
	int* innerIdList;
	int* bond1IdList;
	int* bond2IdList;

};