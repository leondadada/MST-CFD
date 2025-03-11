#include "./MshBlock.h"

MshBlock::MshBlock() {
}
void MshBlock::readMsh(string _mshName) {
	this->mshName = _mshName;
	cout << "\n~~~~~~~~~~~~~~ IO test2 start~~~~~~~~~\n";
	if (DIM == 2){
		splitInformationData2D();
		dataToMesh2D();
	}
	else {
		splitInformationData3D();
		dataToMesh3D();
	}


}



Node* MshBlock::getBeginItNodesList() {
	return nodesList;
}
Cell* MshBlock::getBeginItCellsList() {
	return cellsList;
}
Face* MshBlock::getBeginItFacesList() {
	return facesList;
}
Cell* MshBlock::getBeginItFakeCellsList(){
	return fakeCellsList;
}
std::vector<FacesInf>::iterator MshBlock::getBeginItFacesInfList() {
	return this->facesInfList.begin();
}
std::vector<FacesInf>::iterator MshBlock::getEndItFacesInfList() {
	return this->facesInfList.end();
}

int* MshBlock::getBeginItInnerIdList(){
	return innerIdList;
}
int* MshBlock::getBeginItBond1IdList(){
	return bond1IdList;
}
int* MshBlock::getBeginItBond2IdList(){
	return bond2IdList;
}
int MshBlock::getNumOfNodes() {
	return this->numOfNodes;
}
int MshBlock::getNumOfCells() {
	return this->numOfCells;
}
int MshBlock::getNumOfFaces() {
	return this->numOfFaces;
}
int MshBlock::getNumOfIntFaces() {
	return this->numOfIntFaces;
}
int MshBlock::getNumOfFacesInfs(){
	return this->facesInfList.size();
}
int MshBlock::getNumOfInnerCells(){
	return innerIdListSize;
}
int MshBlock::getNumOfBond1Cells(){
	return bond1IdListSize;
}
int MshBlock::getNumOfBond2Cells(){
	return bond2IdListSize;
}

int MshBlock::splitInformationData2D(){
	cout << "\n~~~~~~~~~~~~~~ MshBlock::splitInformationData() ~~~~~~~~~~~~" << endl;
	ifstream mshFile;
	mshFile.open(mshName, ios::in);
	if (!mshFile.is_open())	{
		cout << "\n!!!!!!!!Didn\'t open mshFilel!!!!!!!!" << endl;
		return 0;
	}

	while (!mshFile.eof()){
		char buffer[512];
		mshFile.getline(buffer, 512);
		stringLine = buffer;
		switch (stringLine[0]  ){
			case '(':{mshInformation += "\n"+stringLine ;break;}
			case ')':{mshInformation += stringLine ;break;}
			default :{mshData += "\n"+stringLine;}
		}
	}

	stringstream meshInformation ;
	meshInformation<<this->mshInformation;
	mshInformation = "";//need to read by line ,this is nacessary

	int headFlag = 0;
	getline(meshInformation,stringLine);
	while (!meshInformation.eof()){
		getline(meshInformation,stringLine);
		switch (stringLine[1]  ){
			case ')':{mshInformation += stringLine ;break;}
			default :{mshInformation += "\n"+stringLine;}
		}
	}
	//cout << mshData << endl;
	return 0;
}
int MshBlock::dataToMesh2D(){
	cout<<"\n~~~~~~~~~~~~~~ MshBlock::dataToMesh() ~~~~~~"<<endl;

	stringstream meshInformation ;
	stringstream meshData;
	stringstream stringstreamLine;
	meshInformation<<this->mshInformation;
	meshData<<this->mshData;

	int flagFace = 0;
	int flagNode = 0;
	int flagCell = 0;

	vector<int> labelSectionClass;
	vector<int> labelFaceClass;
	string txtInformation;
	string tempString;
	vector<string> mshNodesPosition;//all xyz 
	vector<string> mshIntFaces;//all interior  faces
	vector<string> mshBoundaryFaces;//all boundary faces
	
	getline(meshInformation,stringLine);//delete the " " on the first line
	while (!meshInformation.eof()){
		getline(meshInformation,stringLine);

		switch (stringLine[1]){
		case '0':{
			txtInformation = betweenQuotes(stringLine);
			break;
		}
		case '2':{
			//#define DIM (stringLine[3]-'0');//not sure will it be right
			if (stringLine[3] == 2){
				dimMesh = 2;
			}
			else {
				dimMesh = 3;
			}
		}
		case '1':{
			string inBrackets = betweenSecBrackets(stringLine);
			stringstream inBracketsStream;
			inBracketsStream << inBrackets;

			if(inBrackets[0] == '0') {//sum inf of node \cell\face
				if (stringLine[2] == '0') {
					string nodesString;
					getline(inBracketsStream, nodesString, ' ');//need record the number of section
					getline(inBracketsStream, nodesString, ' ');
					getline(inBracketsStream, nodesString, ' ');
					numOfNodes = hexStringToInt(nodesString);
					nodesList = new Node[numOfNodes];
					for (int iNodeId = 0; iNodeId < numOfNodes; iNodeId++){
						nodesList[iNodeId].setId(iNodeId);//id will start from 0;
					}
				}
				if (stringLine[2] == '2') {
					string cellsString;
					getline(inBracketsStream, cellsString, ' ');//need record the number of section
					getline(inBracketsStream, cellsString, ' ');
					getline(inBracketsStream, cellsString, ' ');
					numOfCells = hexStringToInt(cellsString);
					cellsList = new Cell[numOfCells];
					for (int iCellId = 0; iCellId < numOfCells; iCellId++) {
						cellsList[iCellId].setId(iCellId);//id will start from 0;
					}
				}
				if (stringLine[2] == '3') {
					string facesString;
					getline(inBracketsStream, facesString, ' ');//need record the number of section
					getline(inBracketsStream, facesString, ' ');
					getline(inBracketsStream, facesString, ' ');
					numOfFaces = hexStringToInt(facesString);
					facesList = new Face[numOfFaces];
				}
				break;
			}//
			switch (stringLine[2]) {
			case '0': {
				
				string stringNum;
				getline(meshData, stringLineData);
				for (; flagNode < numOfNodes; flagNode++) {
					getline(meshData, stringLineData);
					stringstream streamLineData;
					streamLineData << stringLineData;
					for (int iDim = 0; iDim < DIM; iDim++) {
						getline(streamLineData, stringNum, ' ');
						nodesList[flagNode].setPoint(iDim, (NUM)stod(stringNum));
					}
					
				}
				break;
			}
			case '3': {
				FacesInf x;
				x.setName(txtInformation);
				string inBrackets = betweenSecBrackets(stringLine);
				stringstream inBracketsStream;
				inBracketsStream << inBrackets;
				string facesString;
				getline(inBracketsStream, facesString, ' ');//need record the number of section
				x.setNumber(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//
				x.setStart(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//
				x.setEnd(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//type
				x.setType(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//num of nbNode
				int numNbNode = hexStringToInt(facesString);
				facesInfList.push_back(x);//add inf to list 
				if (x.getType() == 2) {
					numOfIntFaces = x.getEnd();
					//cout << "numOfIntFaces"<<numOfIntFaces;
				}
				string stringNum;
				for (; flagFace < x.getEnd(); flagFace++) {
					getline(meshData, stringLineData);
					stringstream streamLineData;
					streamLineData << stringLineData;
					for (int iNbNode = 0; iNbNode < numNbNode; iNbNode++) {
						getline(streamLineData, stringNum, ' '); //add nbNode on face 
						facesList[flagFace].addNbNode(&(nodesList[hexStringToInt(stringNum) - 1]));
						nodesList[hexStringToInt(stringNum) - 1].addNbFace(&(facesList[flagFace]));
					}
					getline(streamLineData, stringNum, ' '); //
					facesList[flagFace].addNbCell(&(cellsList[hexStringToInt(stringNum) - 1]));
					cellsList[hexStringToInt(stringNum) - 1].addNbFaces(&(facesList[flagFace]));
					facesList[flagFace].setType(x.getType());

					if (x.getType() != 2) {//bond face cell
						cellsList[hexStringToInt(stringNum) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum) - 1]));
						if (!cellsList[hexStringToInt(stringNum) - 1].getFalgBoundary()) {
							bond1IdListVector.push_back(hexStringToInt(stringNum) - 1);
							// cout << hexStringToInt(stringNum) - 1 << " ";
							cellsList[hexStringToInt(stringNum) - 1].setFlagBoundary();
						}

					}
					else {// only int face have two cell sides
						string stringNum1 = stringNum;
						getline(streamLineData, stringNum, ' '); //cout << stringNum << "///" << endl;
						facesList[flagFace].addNbCell(&(cellsList[hexStringToInt(stringNum) - 1]));
						cellsList[hexStringToInt(stringNum) - 1].addNbFaces(&(facesList[flagFace]));
						cellsList[hexStringToInt(stringNum1) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum) - 1]));
						cellsList[hexStringToInt(stringNum) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum1) - 1]));
						// cout << cellsList[hexStringToInt(stringNum1) - 1].getId() << endl;
					}
					facesList[flagFace].setId(flagFace);//id will start from 0;

				}
				break;
			}
			default: {break; }
			}
		}
		case '3': {
		}
		}
	}
	fakeCellsList = new Cell[numOfFaces-numOfIntFaces];//make fake cells on boundary 
	for (int i = 0; i < numOfFaces; i++){
		facesList[i].setCenter();
	}//update center of all faces
	for (int i = 0; i < numOfCells; i++){
		cellsList[i].setCenter();
		cellsList[i].setVolume();
		//cout << itCell->getVolume() << endl;		
	}//update center of all cells 
	for (int i = 0; i < numOfFaces; i++){
		if ((facesList[i].getDirect()).transpose() * (facesList[i].getCenter() - (*facesList[i].getBeginItPNbCells())->getCenter()) < 0) {
			facesList[i].transDirectAndCells();//find if the direct point from cell[0] to cell[1];
			for (int iDim = 0; iDim < DIM ; iDim++) {
				//if   FlagLeftRight = true
				if (facesList[i].getDirectAndCells() * facesList[i].getDirect()[iDim] >= 0) {
					facesList[i].setFlagLeftRight(iDim, true);
				}
				else {
					facesList[i].setFlagLeftRight(iDim,false );
				}
			}			
		}
		else {
			for (int iDim = 0; iDim < DIM; iDim++) {
				if (facesList[i].getDirectAndCells() * facesList[i].getDirect()[iDim] >= 0) {
					facesList[i].setFlagLeftRight(iDim, false);
				}
				else {
					facesList[i].setFlagLeftRight(iDim, true);
				}
			}
		}
		// cout <<" "<< itFace->getDirect().transpose()*((itFace->getBeginItPNbNodes())[0]->getCenter() - (itFace->getBeginItPNbNodes())[1]->getCenter());
	}
	
	for (int i = 0; i != numOfCells; i++) {
		cellsList[i].creatDirectOfNbFaces();
		auto itBeginNbFace = cellsList[i].getBeginItPNbFaces();
		for (int iPNbFace = 0; iPNbFace < cellsList[i].getNumOfNbFaces(); iPNbFace++) {
			if (&(cellsList[i]) == *(itBeginNbFace[iPNbFace]->getBeginItPNbCells())) {
				cellsList[i].addDirectOfNbFaces(iPNbFace, itBeginNbFace[iPNbFace]->getDirectAndCells()*itBeginNbFace[iPNbFace]->getDirect());
			}
			else {
				cellsList[i].addDirectOfNbFaces(iPNbFace, -itBeginNbFace[iPNbFace]->getDirectAndCells()*itBeginNbFace[iPNbFace]->getDirect());
			}
		}
	}//update DirectOfNbFaces of all cells 
	for (int i= 0; i < numOfFaces; i++) {
		facesList[i].updateEta0();
	}
	
	for (int i = 0; i < numOfCells; i++) {
		cellsList[i].creatEta0C();
		for (int iNbFace = 0; iNbFace < cellsList[i].getNumOfNbFaces(); iNbFace++){
			if (cellsList[i].getBeginItPNbFaces()[iNbFace]->getBeginItPNbCells()[0] == &(cellsList[i])) {
				cellsList[i].addEta0C(iNbFace,cellsList[i].getBeginItPNbFaces()[iNbFace]->getEta0() );
			}
			else {
				cellsList[i].addEta0C(iNbFace,1-cellsList[i].getBeginItPNbFaces()[iNbFace]->getEta0());

			}
		}				
	}
	for (int iCell = 0; iCell < numOfCells; iCell++) {
		//cout << "iCell" <<iCell<< endl;
		auto itNbFace = cellsList[iCell].getBeginItPNbFaces();
		for (int iNbFace = 0; iNbFace<cellsList[iCell].getNumOfNbFaces(); iNbFace++) {
			//cout << "iNbFace" <<(*iNbFace )->getId()<< endl;
			auto itNbNodeOfFace = itNbFace[iNbFace]->getBeginItPNbNodes();
			for (int iNbNodeOfFace = 0; iNbNodeOfFace< itNbFace[iNbFace]->getNumOfpNbNodes() ; iNbNodeOfFace++) {
				int flagHaveNode = 0;
				auto itNodeOfCell = cellsList[iCell].getBeginItPNbNodes();
				for (int iNbNode = 0 ; iNbNode < cellsList[iCell].getNumOfNbNodes(); iNbNode++) {
					
					if (itNodeOfCell[iNbNode]->getId() == itNbNodeOfFace[iNbNodeOfFace]->getId()) {
						flagHaveNode = 1;
					}
				}
				if (flagHaveNode == 0) {
					cellsList[iCell].addNbNodes(itNbNodeOfFace[iNbNodeOfFace]);
					//cout << "addNbNodes " << (*iNbNodeOfFace)->getId() << endl;
				}
				//cout << "flagHaveNode = "<< flagHaveNode << endl;
			}
		}
		if (cellsList[iCell].getNumOfNbNodes() == 4) {
			auto pBeginNbNode = cellsList[iCell].getBeginItPNbNodes();
			VCTDIM vct1 = pBeginNbNode[0]->getCenter() - pBeginNbNode[1]->getCenter();
			VCTDIM vct2 = pBeginNbNode[2]->getCenter() - pBeginNbNode[3]->getCenter();
			if (vct1.transpose()*vct2 < 0) {
				auto temp = pBeginNbNode[3];
				pBeginNbNode[3] = pBeginNbNode[2];
				pBeginNbNode[2] = temp;
				cout << cellsList[iCell].getCenter() << endl;
			}
		}
	}
	
	for (int iBond1 = 0; iBond1 < bond1IdListVector.size(); iBond1++) {//set bond2idlist
		int iIdBond1 = bond1IdListVector[iBond1];
		for (int iNbCell = 0; iNbCell < cellsList[iIdBond1].getNumOfNbCells(); iNbCell++){
			if (cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->getFalgBoundary() == 0) {
				bond2IdListVector.push_back(cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->getId());
				cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->setFlagBoundary();
				//cout << cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->getId() << endl;
			}
			//cout << cellsList[iIdBond1].getBeginItPNbFaces()[iNbCell]->getType() <<" ";
		}
		//cout << endl;
	}
	innerIdList = new int[numOfCells];
	int flagInnerCells = 0;
	for (int iCell = 0; iCell < numOfCells; iCell++){//set inneridlist		
		if (cellsList[iCell].getFalgBoundary() == 0) {
			innerIdList[flagInnerCells] = (iCell);
			//cellsList[iCell].setFlagBoundary();
			flagInnerCells++;
		}
	}

	innerIdListSize = flagInnerCells;
	bond1IdListSize = bond1IdListVector.size();
	bond2IdListSize = bond2IdListVector.size();
	bond1IdList = new int[bond1IdListVector.size()];
	bond2IdList = new int[bond2IdListVector.size()];

	for (int iBondCell1 = 0; iBondCell1 < bond1IdListVector.size(); iBondCell1++){
		bond1IdList[iBondCell1] = bond1IdListVector[iBondCell1];
	}
	for (int iBondCell2 = 0; iBondCell2 < bond2IdListVector.size(); iBondCell2++) {
		bond2IdList[iBondCell2] = bond2IdListVector[iBondCell2];
	}
	
	bond1IdListVector.clear();
	bond2IdListVector.clear();

	
	
	if (1) {

		fstream mshInformationFile;
		fstream mshDataFile;
		mshInformationFile.open( mshName+".mshinf", ios::out | ios::trunc);
		mshDataFile.open( mshName+".mshdat", ios::out | ios::trunc);
		mshInformationFile<<mshInformation;
		mshDataFile<<mshData;//write document ,not neccessary

	}
	


	







	return 0; 
}



int MshBlock::splitInformationData3D() {
	cout << "\n~~~~~~~~~~~~~~ MshBlock::splitInformationData() ~~~~~~~~~~~~" << endl;
	ifstream mshFile;
	mshFile.open(mshName, ios::in);
	if (!mshFile.is_open()) {
		cout << "\n!!!!!!!!Didn\'t open mshFilel!!!!!!!!" << endl;
		return 0;
	}

	while (!mshFile.eof()) {
		char buffer[512];
		mshFile.getline(buffer, 512);
		stringLine = buffer;
		switch (stringLine[0]) {
		case '(': {mshInformation += "\n" + stringLine; break; }
		case ')': {mshInformation += stringLine; break; }
		default: {mshData += "\n" + stringLine; }
		}
	}

	stringstream meshInformation;
	meshInformation << this->mshInformation;
	mshInformation = "";//need to read by line ,this is nacessary

	int headFlag = 0;
	getline(meshInformation, stringLine);
	while (!meshInformation.eof()) {
		getline(meshInformation, stringLine);
		switch (stringLine[1]) {
		case ')': {mshInformation += stringLine; break; }
		default: {mshInformation += "\n" + stringLine; }
		}
	}
	//cout << mshData << endl;
	return 0;
}
int MshBlock::dataToMesh3D() {
	cout << "\n~~~~~~~~~~~~~~ MshBlock::dataToMesh() ~~~~~~" << endl;

	stringstream meshInformation;
	stringstream meshData;
	stringstream stringstreamLine;
	meshInformation << this->mshInformation;
	meshData << this->mshData;

	int flagFace = 0;
	int flagNode = 0;
	int flagCell = 0;

	vector<int> labelSectionClass;
	vector<int> labelFaceClass;
	string txtInformation;
	string tempString;
	vector<string> mshNodesPosition;//all xyz 
	vector<string> mshIntFaces;//all interior  faces
	vector<string> mshBoundaryFaces;//all boundary faces

	getline(meshInformation, stringLine);//delete the " " on the first line
	while (!meshInformation.eof()) {
		getline(meshInformation, stringLine);

		switch (stringLine[1]) {
		case '0': {
			txtInformation = betweenQuotes(stringLine);
			break;
		}
		case '2': {
			//#define DIM (stringLine[3]-'0');//not sure will it be right
		}
		case '1': {
			string inBrackets = betweenSecBrackets(stringLine);
			stringstream inBracketsStream;
			inBracketsStream << inBrackets;

			if (inBrackets[0] == '0') {//sum inf of node \cell\face
				if (stringLine[2] == '0') {
					string nodesString;
					getline(inBracketsStream, nodesString, ' ');//need record the number of section
					getline(inBracketsStream, nodesString, ' ');//start
					getline(inBracketsStream, nodesString, ' ');//end
					numOfNodes = hexStringToInt(nodesString);
					nodesList = new Node[numOfNodes];
					for (int iNodeId = 0; iNodeId < numOfNodes; iNodeId++) {
						nodesList[iNodeId].setId(iNodeId);//id will start from 0;
					}
				}
				if (stringLine[2] == '2') {//number of cells 
					string cellsString;
					getline(inBracketsStream, cellsString, ' ');//need record the number of section
					getline(inBracketsStream, cellsString, ' ');//start
					getline(inBracketsStream, cellsString, ' ');//end
					numOfCells = hexStringToInt(cellsString);
					cellsList = new Cell[numOfCells];
					for (int iCellId = 0; iCellId < numOfCells; iCellId++) {
						cellsList[iCellId].setId(iCellId);//id will start from 0;
					}
				}
				if (stringLine[2] == '3') {//all faces
					string facesString;
					getline(inBracketsStream, facesString, ' ');//need record the number of section
					getline(inBracketsStream, facesString, ' ');
					getline(inBracketsStream, facesString, ' ');
					numOfFaces = hexStringToInt(facesString);
					facesList = new Face[numOfFaces];
				}
				break;
			}//
			switch (stringLine[2]) {
			case '0': {//node inf 

				string stringNum;
				getline(meshData, stringLineData);
				for (; flagNode < numOfNodes; flagNode++) {
					getline(meshData, stringLineData);
					stringstream streamLineData;
					streamLineData << stringLineData;
					for (int iDim = 0; iDim < DIM; iDim++){
						getline(streamLineData, stringNum, ' ');
						nodesList[flagNode].setPoint(iDim, (NUM)stod(stringNum));
					}

				}
				break;
			}
			case '3': {
				FacesInf x;
				x.setName(txtInformation);
				string inBrackets = betweenSecBrackets(stringLine);
				stringstream inBracketsStream;
				inBracketsStream << inBrackets;
				string facesString;
				getline(inBracketsStream, facesString, ' ');//need record the number of section
				x.setNumber(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//
				x.setStart(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//
				x.setEnd(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//type
				x.setType(hexStringToInt(facesString));
				getline(inBracketsStream, facesString, ' ');//num of nbNode
				int numNbNode = hexStringToInt(facesString);
				facesInfList.push_back(x);//add inf to list 
				if (x.getType() == 2) {
					numOfIntFaces = x.getEnd();
					//cout << "numOfIntFaces"<<numOfIntFaces;
				}
				string stringNum;
				for (; flagFace < x.getEnd(); flagFace++) {
					getline(meshData, stringLineData);
					stringstream streamLineData;
					streamLineData << stringLineData;
					for (int iNbNode = 0; iNbNode < numNbNode; iNbNode++){
						getline(streamLineData, stringNum, ' '); //add nbNode on face 
						facesList[flagFace].addNbNode(&(nodesList[hexStringToInt(stringNum) - 1]));
						nodesList[hexStringToInt(stringNum) - 1].addNbFace(&(facesList[flagFace]));
					}					
					getline(streamLineData, stringNum, ' '); //
					facesList[flagFace].addNbCell(&(cellsList[hexStringToInt(stringNum) - 1]));
					cellsList[hexStringToInt(stringNum) - 1].addNbFaces(&(facesList[flagFace]));
					facesList[flagFace].setType(x.getType());

					if (x.getType() != 2) {//bond face cell
						cellsList[hexStringToInt(stringNum) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum) - 1]));
						if (!cellsList[hexStringToInt(stringNum) - 1].getFalgBoundary()) {
							bond1IdListVector.push_back(hexStringToInt(stringNum) - 1);
							// cout << hexStringToInt(stringNum) - 1 << " ";
							cellsList[hexStringToInt(stringNum) - 1].setFlagBoundary();
						}

					}
					else {// only int face have two cell sides
						string stringNum1 = stringNum;
						getline(streamLineData, stringNum, ' '); //cout << stringNum << "///" << endl;
						facesList[flagFace].addNbCell(&(cellsList[hexStringToInt(stringNum) - 1]));
						cellsList[hexStringToInt(stringNum) - 1].addNbFaces(&(facesList[flagFace]));
						cellsList[hexStringToInt(stringNum1) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum) - 1]));
						cellsList[hexStringToInt(stringNum) - 1].addNbCells(&(cellsList[hexStringToInt(stringNum1) - 1]));
						// cout << cellsList[hexStringToInt(stringNum1) - 1].getId() << endl;
					}
					facesList[flagFace].setId(flagFace);//id will start from 0;
					
				}
				break;
			}
			default: {break; }
			}
		}
		case '3': {
		}
		}
	}
	fakeCellsList = new Cell[numOfFaces - numOfIntFaces];//make fake cells on boundary 
	for (int i = 0; i < numOfFaces; i++) {
		facesList[i].setCenter();
	}//update center of all faces
	for (int i = 0; i < numOfCells; i++) {
		cellsList[i].setCenter();
		cellsList[i].setVolume();
	}//update center of all cells 
	for (int i = 0; i < numOfFaces; i++) {
		if ((facesList[i].getDirect()).transpose() * (facesList[i].getCenter() - (*facesList[i].getBeginItPNbCells())->getCenter()) < 0) {
			facesList[i].transDirectAndCells();//find if the direct point from cell[0] to cell[1];
			for (int iDim = 0; iDim < DIM; iDim++) {
				if (facesList[i].getDirectAndCells() * facesList[i].getDirect()[iDim] >= 0) {
					facesList[i].setFlagLeftRight(iDim, true);
				}
				else {
					facesList[i].setFlagLeftRight(iDim, false);
				}
			}
		}
	}

	for (int i = 0; i < numOfCells; i++) {
		cellsList[i].creatDirectOfNbFaces();
		auto itBeginNbFace = cellsList[i].getBeginItPNbFaces();
		for (int iPNbFace = 0; iPNbFace  < cellsList[i].getNumOfNbFaces(); iPNbFace++) {
			if (&(cellsList[i]) == *(itBeginNbFace[iPNbFace]->getBeginItPNbCells())) {
				cellsList[i].addDirectOfNbFaces(iPNbFace,itBeginNbFace[iPNbFace]->getDirectAndCells()*itBeginNbFace[iPNbFace]->getDirect());
			}
			else {
				cellsList[i].addDirectOfNbFaces(iPNbFace ,-itBeginNbFace[iPNbFace]->getDirectAndCells()*itBeginNbFace[iPNbFace]->getDirect());
			}
		}		
	}//update DirectOfNbFaces of all cells 
	for (int i = 0; i < numOfFaces; i++) {
		facesList[i].updateEta0();
	}
	for (int i = 0; i < numOfIntFaces; i++) {
		facesList[i].updateCenterMiddleDef();
	}//center middle defferecen calculation 
	
	for (int i = 0; i < numOfCells; i++) {
		cellsList[i].creatEta0C();
		for (int iNbFace = 0; iNbFace < cellsList[i].getNumOfNbFaces(); iNbFace++) {

			if (cellsList[i].getBeginItPNbFaces()[iNbFace]->getBeginItPNbCells()[0] == &(cellsList[i])) {
				cellsList[i].addEta0C(iNbFace, cellsList[i].getBeginItPNbFaces()[iNbFace]->getEta0());
			}
			else {
				cellsList[i].addEta0C(iNbFace, 1 - cellsList[i].getBeginItPNbFaces()[iNbFace]->getEta0());

			}
		}	
	}
	for (int iCell = 0; iCell < numOfCells; iCell++) {
		//cout << "iCell" <<iCell<< endl;
		auto itNbFace = cellsList[iCell].getBeginItPNbFaces();
		for (int iNbFace = 0; iNbFace<cellsList[iCell].getNumOfNbFaces(); iNbFace++) {
			//cout << "iNbFace" <<(*iNbFace )->getId()<< endl;
			auto itNbNodeOfFace = itNbFace[iNbFace]->getBeginItPNbNodes();
			for (int iNbNodeOfFace = 0; iNbNodeOfFace< itNbFace[iNbFace]->getNumOfpNbNodes(); iNbNodeOfFace++) {
				int flagHaveNode = 0;
				auto itNodeOfCell = cellsList[iCell].getBeginItPNbNodes();
				for (int iNbNode = 0; iNbNode < cellsList[iCell].getNumOfNbNodes(); iNbNode++) {
					if (itNodeOfCell[iNbNode]->getId() == itNbNodeOfFace[iNbNodeOfFace]->getId()) {
						flagHaveNode = 1;
					}
				}
				if (flagHaveNode == 0) {
					cellsList[iCell].addNbNodes(itNbNodeOfFace[iNbNodeOfFace]);
				}
			}
		}

		if (cellsList[iCell].getNumOfNbNodes() == 4) {//struct structure grid mesh 
			auto pBeginNbNode = cellsList[iCell].getBeginItPNbNodes();
			VCTDIM vct1 = pBeginNbNode[1]->getCenter() - pBeginNbNode[0]->getCenter();
			VCTDIM vct2 = pBeginNbNode[2]->getCenter() - pBeginNbNode[0]->getCenter();
			VCTDIM vct3 = pBeginNbNode[3]->getCenter() - pBeginNbNode[0]->getCenter();
			double angle12 = acos(vct1.dot(vct2) / (vct1.norm() * vct2.norm()));
			double angle23 = acos(vct3.dot(vct2) / (vct3.norm() * vct2.norm()));
			double angle13 = acos(vct1.dot(vct3) / (vct1.norm() * vct3.norm()));

			if (angle13 < angle12 || angle13 < angle23) {//true:not 0123
				if (angle12 > angle23) {
					auto temp = pBeginNbNode[3];
					pBeginNbNode[3] = pBeginNbNode[2];
					pBeginNbNode[2] = temp;
				}
				else {
					auto temp = pBeginNbNode[1];
					pBeginNbNode[1] = pBeginNbNode[2];
					pBeginNbNode[2] = temp;
				}


			}

		}
	}

	for (int iBond1 = 0; iBond1 < bond1IdListVector.size(); iBond1++) {//set bond2idlist
		int iIdBond1 = bond1IdListVector[iBond1];
		for (int iNbCell = 0; iNbCell < cellsList[iIdBond1].getNumOfNbCells(); iNbCell++) {
			if (cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->getFalgBoundary() == 0) {
				bond2IdListVector.push_back(cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->getId());
				cellsList[iIdBond1].getBeginItPNbCells()[iNbCell]->setFlagBoundary();
			}
		}
	}
	innerIdList = new int[numOfCells];
	int flagInnerCells = 0;
	for (int iCell = 0; iCell < numOfCells; iCell++) {//set inneridlist		
		if (cellsList[iCell].getFalgBoundary() == 0) {
			innerIdList[flagInnerCells] = (iCell);
			cellsList[iCell].setFlagBoundary();
			flagInnerCells++;
		}
	}
	for (int iCell = 0; iCell < numOfCells; iCell++) {//20230217-set L-U nbCells
		auto itBeginNbCell = cellsList[iCell].getBeginItPNbCells();
		for (int iNbCell = 0; iNbCell < cellsList[iCell].getNumOfNbFaces(); iNbCell++) {
			if (cellsList[iCell].getId() > cellsList[iCell].getBeginItPNbCells()[iNbCell]->getId()) {
				cellsList[iCell].addNbCellsL(cellsList[iCell].getBeginItPNbCells()[iNbCell]);
			}
			else {
				cellsList[iCell].addNbCellsU(cellsList[iCell].getBeginItPNbCells()[iNbCell]);
			}
		}
	}

	innerIdListSize = flagInnerCells;
	bond1IdListSize = bond1IdListVector.size();
	bond2IdListSize = bond2IdListVector.size();
	bond1IdList = new int[bond1IdListVector.size()];
	bond2IdList = new int[bond2IdListVector.size()];

	for (int iBondCell1 = 0; iBondCell1 < bond1IdListVector.size(); iBondCell1++) {
		bond1IdList[iBondCell1] = (int)bond1IdListVector[iBondCell1];
	}
	for (int iBondCell2 = 0; iBondCell2 < bond2IdListVector.size(); iBondCell2++) {
		bond2IdList[iBondCell2] = (int)bond2IdListVector[iBondCell2];
	}

	bond1IdListVector.clear();
	bond2IdListVector.clear();




	if (1) {

		fstream mshInformationFile;
		fstream mshDataFile;
		mshInformationFile.open(mshName + ".mshinf", ios::out | ios::trunc);
		mshDataFile.open(mshName + ".mshdat", ios::out | ios::trunc);
		mshInformationFile << mshInformation;
		mshDataFile << mshData;//write document ,not neccessary

	}

	return 0;
}

