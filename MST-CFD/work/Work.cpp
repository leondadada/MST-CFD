#include "Work.h"

Work::~Work() {

}

void Work::work(string _mshName) {

	mshAddress = "./";
	mshName = _mshName;
	
	string title;
	string num2;
	string num3;
	string num4;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvEP";

	ss2 << STEP_TIME; ss2 >> num2;
	ss4 << inletu; ss4 >> num4;
	title = mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "-log.lhblog";

	fLog.open(title, ios::out | ios::trunc);
	mesh.readMsh(mshAddress + "msh/" + mshName);
	

	allData.createAllData(mesh.getNumOfCells(), mesh.getNumOfFaces());

	Time time1(&mesh, &fLog,&allData);
	VCTDIMU iniQ;
#if DIM == 2
	iniQ << inirho, inirho*iniu, inirho*iniv, iniE;
#else 
	iniQ << inirho, inirho*iniu, inirho*iniv, inirho*iniw, iniE;
#endif
	time1.initialization(iniQ);
	

#if RHO_P == 0 
		writedataRhoBasedMshNodePlt(&mesh, allData.getP1OldCellQs(), 0);
	
#else 
		writedataPBasedMshNodePlt(&mesh, allData.getPOldCellPUVT(), 0);
	
#endif

#if READ_PLT != 0
		readPbasedCellPlt(_mshName);
#endif
		
	for (; true; t++) {
		cout << "FVM_2D_MSH_" << STEP_TIME << " : " << t << "-->" << t + 1 << endl;
	
		time1.goNextTimeStep();
		cout << "time now :"<<t<<"-->"<<t + 1 << endl;
		if (READ_PLT != 0) {
			if (t < READ_PLT) {
				continue;
			}
		}

		if ((t + 1) % (int)round(SAVE_TIME * STEP_TIME) == 0) {
			/*save per SAVE_TIMEs */
#if RHO_P == 0	
			writedataRhoBasedMshNodePlt(&mesh, allData.getP1OldCellQs(), (t + 1));

#else 
			writedataPBasedMshNodePlt(&mesh, allData.getPOldCellPUVT(), (t + 1));

#endif
			if ((t + 1) % (int)round(SAVE_TIME * STEP_TIME*10) == 0){
				writedataCellContinue(_mshName, t+1);

			}
		}//save data .csv


		if ((t + 1) % (int)(TIME * STEP_TIME) == 0) {

			cout << "\n\n\ncontinue to calculate???  \n0:no break out!!!     \nother number:yes ,go on\n";
			int judge;
			cin >> judge;
			if (judge == 0) {
				break;
			}
		}
	}
	

}



void Work::writedataPBasedMshNodePlt(MshBlock * _pMesh, VCTDIMU * _p1PUVT, int _time_now) {
	//_numm choose which property to save
	string title;
	string num2;
	string num3;
	string num4;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvEP";
	ss2 << STEP_TIME; ss2 >> num2;
	ss3 << _time_now; ss3 >> num3;
	ss4 << inletu; ss4 >> num4;
	int numOfCells = mesh.getNumOfCells();
	int numOfFaces = mesh.getNumOfFaces();
	int numOfNodes = mesh.getNumOfNodes();

	title = mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "_t" + num3 + ".plt";
	
	cout << "writing_node: " << title << endl;
	fstream f;
	f.open(title, ios::out | ios::trunc);
	f << "\"TITLE = \"Example: Variable and Connectivity List Sharing\""<<endl;
#if DIM == 2 //dim ==2
		f << "VARIABLES = \"X\", \"Y\", \"p\" , \"u\" , \"v\" , \"T\" , \"rho\" , \"Ma\" " << endl;
		f << "ZONE T=\""<< t <<"\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETRIANGLE" << endl;
	
#else 
		f << "VARIABLES = \"X\", \"Y\", \"Z\", \"p\" , \"u\" , \"v\" , \"T\" , \"rho\" , \"Ma\" " << endl;
		f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETETRAHEDRON" << endl;
	
#endif
	VCTDIMU* pFacePUVT = new VCTDIMU[numOfFaces];
	VCTDIMU* pNodePUVT = new VCTDIMU[numOfNodes];
	
	for (int inf = 0; inf < _pMesh->getNumOfFacesInfs(); inf++) {
		auto itBeginFacesInf = _pMesh->getBeginItFacesInfList();
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update	
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVT[i] = _pMesh->getBeginItFacesList()[i].getEta0()*_p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				pFacePUVT[i] += (1 - _pMesh->getBeginItFacesList()[i].getEta0())*_p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[1]->getId()];
			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVT[i] = _p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {

				if (FLAG_WALLSMOOTH == 1) {
					pFacePUVT[i] = _p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				}
				else {
					pFacePUVT[i] = _p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
					pFacePUVT[i][1] = 0;
					pFacePUVT[i][2] = 0;
				}// if wall smooth is avaliable
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFacePUVT[i] = _p1PUVT[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		}
	}//all face (interior & boundary) update

	auto itBeginNodes = _pMesh->getBeginItNodesList();
	for (int i = 0; i < numOfNodes; i++) {
		pNodePUVT[i] = VCTDIMU::Zero();
		NUM d = 0;
		for (auto itFaceNb = itBeginNodes[i].getBeginItPNbFaces(); itFaceNb < itBeginNodes[i].getEndItPNbFaces(); itFaceNb++) {
			pNodePUVT[i] += 1 / mesh.getBeginItFacesList()[i].getArea() *pFacePUVT[(*itFaceNb)->getId()];
			d += 1 / mesh.getBeginItFacesList()[i].getArea();
		}
		pNodePUVT[i] /= d;
		for (int iDim = 0; iDim < DIM; iDim++){
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << itBeginNodes[i].getCenter()[iDim] << " ";
		}
		f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodePUVT[i][0] << " ";
		for(int iDim = 0; iDim < DIM; iDim++) {
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodePUVT[i][iDim+1] << " ";
		}
		f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodePUVT[i][DIMU - 1] << " " << getPBasedRho( pNodePUVT[i])<<" "<< getPBasedMa(pNodePUVT[i]);
		f << endl;
	}
	for (int i = 0; i < numOfCells; i++) {
		auto pBeginNbNode = _pMesh->getBeginItCellsList()[i].getBeginItPNbNodes();
		for (auto iNbNode = 0; iNbNode < _pMesh->getBeginItCellsList()[i].getNumOfNbNodes(); iNbNode++) {
			f << pBeginNbNode[iNbNode]->getId() + 1 << ' ';
		}
		f << endl;
	}



	f.close();

	delete[] pNodePUVT;
	delete[] pFacePUVT;
}
void Work::writedataRhoBasedMshNodePlt(MshBlock * _pMesh, VCTDIMU * _p1CellQ, int _time_now) {
	//_numm choose which property to save
	string title;
	string num2;
	string num3;
	string num4;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvTP";
	ss2 << STEP_TIME; ss2 >> num2;
	ss3 << _time_now; ss3 >> num3;
	ss4 << inletu; ss4 >> num4;
	int numOfCells = mesh.getNumOfCells();
	int numOfFaces = mesh.getNumOfFaces();
	int numOfNodes = mesh.getNumOfNodes();

	title = mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "_t" + num3 + ".plt";

	cout << "writing_node: " << title << endl;
	//fLog << "writing_node: " << title << endl;
	fstream f;
	f.open(title, ios::out | ios::trunc);
	f << "\"TITLE = \"Example: Variable and Connectivity List Sharing\"" << endl;
#if DIM == 2 //dim ==2
		f << "VARIABLES = \"X\", \"Y\", \"rho\" , \"u\" , \"v\" , \"T\" , \"p\" , \"Ma\" " << endl;

#if FELNUM == 3 
			f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETRIANGLE" << endl;
#else 
			f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FEQUADRILATERAL" << endl;
#endif
	
#else 
		f << "VARIABLES = \"X\", \"Y\", \"Z\", \"rho\" , \"u\" , \"v\" , \"T\" , \"p\" , \"Ma\" " << endl;
		f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETETRAHEDRON" << endl;
	
#endif

	VCTDIMU* pFaceQ = new VCTDIMU[numOfFaces];
	VCTDIMU* pNodeQ = new VCTDIMU[numOfNodes];



	for (int inf = 0; inf < _pMesh->getNumOfFacesInfs(); inf++) {
		auto itBeginFacesInf = _pMesh->getBeginItFacesInfList();
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _pMesh->getBeginItFacesList()[i].getEta0()*_p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				pFaceQ[i] += (1- _pMesh->getBeginItFacesList()[i].getEta0())*_p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[1]->getId()];
			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {

				if (FLAG_WALLSMOOTH == 1) {
					pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				}
				else {
					pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
					pFaceQ[i][1] = 0;
					pFaceQ[i][2] = 0;
				}// if wall smooth is avaliable
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		}
	}//all face (interior & boundary) update

	auto itBeginNodes = _pMesh->getBeginItNodesList();
	for (int i = 0; i < numOfNodes; i++) {
		pNodeQ[i] = VCTDIMU::Zero();
		NUM d = 0;
		for (auto itFaceNb = itBeginNodes[i].getBeginItPNbFaces(); itFaceNb < itBeginNodes[i].getEndItPNbFaces(); itFaceNb++) {
			pNodeQ[i] += 1 / mesh.getBeginItFacesList()[i].getArea() *pFaceQ[(*itFaceNb)->getId()];
			d += 1 / mesh.getBeginItFacesList()[i].getArea();
		}
		pNodeQ[i] /= d;
		for (int iDim = 0; iDim < DIM; iDim++) {
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << itBeginNodes[i].getCenter()(iDim) << " ";
		}
		f  << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodeQ[i][0] << " ";
		for (int iDim = 0; iDim < DIM; iDim++) {
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodeQ[i][iDim + 1]/ pNodeQ[i][0] << " ";
		}
		f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << getT( pNodeQ[i]) << " " << getP(pNodeQ[i]) << " " << getMa(pNodeQ[i]) << " ";
		f << "\n";
	}
	for (int i = 0; i < numOfCells; i++) {
		auto pBeginNbNode = _pMesh->getBeginItCellsList()[i].getBeginItPNbNodes();
		for (auto iNbNode = 0; iNbNode < _pMesh->getBeginItCellsList()[i].getNumOfNbNodes(); iNbNode++) {
			f << pBeginNbNode[iNbNode]->getId() + 1 << ' ';
		}
		f << "\n";
	}


	f.close();

	delete[] pNodeQ;
	delete[] pFaceQ;
}
void Work::writedataRhoBasedMshNodeQUADPlt(MshBlock * _pMesh, VCTDIMU * _p1CellQ, int _time_now) {
	//_numm choose which property to save
	string title;
	string num2;
	string num3;
	string num4;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvEP";
	ss2 << STEP_TIME; ss2 >> num2;
	ss3 << _time_now; ss3 >> num3;
	ss4 << inletu; ss4 >> num4;
	int numOfCells = mesh.getNumOfCells();
	int numOfFaces = mesh.getNumOfFaces();
	int numOfNodes = mesh.getNumOfNodes();

	title = mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "_t" + num3 + ".plt";

	cout << "writing_node: " << title << endl;
	//fLog << "writing_node: " << title << endl;
	fstream f;
	f.open(title, ios::out | ios::trunc);
	
	f << "\"TITLE = \"Example: Variable and Connectivity List Sharing\"" << endl;
	if (DIM == 2) {//dim ==2
		f << "VARIABLES = \"X\", \"Y\", \"rho\" , \"u\" , \"v\" , \"T\" , \"p\" , \"Ma\" " << endl;
		f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETETRAHEDRON" << endl;
	}
	else {
		f << "VARIABLES = \"X\", \"Y\", \"Z\", \"rho\" , \"u\" , \"v\" , \"T\" , \"p\" , \"Ma\" " << endl;
		f << "ZONE T=\"" << t << "\", DATAPACKING=POINT, NODES=" << numOfNodes << ", ELEMENTS=" << numOfCells << ", ZONETYPE=FETETRAHEDRON" << endl;
	}

	VCTDIMU* pFaceQ = new VCTDIMU[numOfFaces];
	VCTDIMU* pNodeQ = new VCTDIMU[numOfNodes];


	for (int inf = 0; inf < _pMesh->getNumOfFacesInfs(); inf++) {
		auto itBeginFacesInf = _pMesh->getBeginItFacesInfList();
		switch (itBeginFacesInf[inf].getType()) {
		case 2: {//interior face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _pMesh->getBeginItFacesList()[i].getEta0()*_p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				pFaceQ[i] += (1 - _pMesh->getBeginItFacesList()[i].getEta0())*_p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[1]->getId()];
			}
			continue;
		}
		case 10: {//inlet face update 
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		case 3: {//wall face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {

				if (FLAG_WALLSMOOTH == 1) {
					pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
				}
				else {
					pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
					pFaceQ[i][1] = 0;
					pFaceQ[i][2] = 0;
				}// if wall smooth is avaliable
			}
			continue;
		}
		case 5: {//outlet face update
			for (int i = itBeginFacesInf[inf].getStart() - 1; i < itBeginFacesInf[inf].getEnd(); i++) {
				pFaceQ[i] = _p1CellQ[_pMesh->getBeginItFacesList()[i].getBeginItPNbCells()[0]->getId()];
			}
			continue;
		}
		}
	}//all face (interior & boundary) update

	auto itBeginNodes = _pMesh->getBeginItNodesList();
	for (int i = 0; i < numOfNodes; i++) {
		pNodeQ[i] = VCTDIMU::Zero();
		NUM d = 0;
		for (auto itFaceNb = itBeginNodes[i].getBeginItPNbFaces(); itFaceNb < itBeginNodes[i].getEndItPNbFaces(); itFaceNb++) {
			pNodeQ[i] += 1 / mesh.getBeginItFacesList()[i].getArea() *pFaceQ[(*itFaceNb)->getId()];
			d += 1 / mesh.getBeginItFacesList()[i].getArea();
		}
		pNodeQ[i] /= d;
		for (int iDim = 0; iDim < DIM; iDim++) {
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << itBeginNodes[i].getCenter()[iDim] << " ";
		}
		f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodeQ[i][0];
		for (int iDim = 0; iDim < DIM; iDim++) {
			f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << pNodeQ[i][iDim + 1] / pNodeQ[i][0] << " ";
		}
		f << fixed << setprecision(accuracyFile) << setw(accuracyFile) << getT(pNodeQ[i]) << " " << getP(pNodeQ[i]) << " " << getMa(pNodeQ[i]);

	}	
	for (int i = 0; i < numOfCells; i++) {
		auto iPNode = _pMesh->getBeginItCellsList()[i].getBeginItPNbNodes();
		f << (iPNode[0])->getId() + 1 << ' ';
		f << (iPNode[1])->getId() + 1 << ' ';
		f << (iPNode[3])->getId() + 1 << ' ';
		f << (iPNode[2])->getId() + 1 << ' ';
		f << endl;
	}


	f.close();

	delete[] pNodeQ;
	delete[] pFaceQ;
}

void Work::readPbasedCellPlt(string _mshName) {
	this->mshAddress = "./";
	this->mshName = _mshName;

	string title;
	string num1;
	string num2;
	string num3;
	string num4;
	stringstream ss1;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvEP";

	ss1 << READ_PLT; ss1 >> num1;
	ss2 << STEP_TIME; ss2 >> num2;
	ss4 << inletu; ss4 >> num4;
	title = this->mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "_t" + num1 + ".cellPUVT";

	fstream readFile;
	readFile.open(title, ios::in);
	if (!readFile.is_open()) {
		cout << "\n!!!!!!!!Didn\'t open readFilel!!!!!!!!" << endl;
		cout << title << endl;
		//return ;
	}

	string stringLine;		
	for (int iCell = 0; iCell < mesh.getNumOfCells(); iCell++) {

		string stringLineData;
		getline(readFile, stringLineData);
		stringstream streamLineData;
		streamLineData << stringLineData;
		string stringNum;
		for (size_t iDimu = 0; iDimu < DIMU; iDimu++){
						
			getline(streamLineData, stringNum, ' ');
			allData.getPNewCellPUVT()[iCell][iDimu] = (NUM)stod(stringNum);

		}
		



	}
	
	cout << stringLine << endl;

}
void Work::writedataCellContinue(string _mshName,int _numm) {
	this->mshAddress = "./";
	this->mshName = _mshName;

	string title;
	string num1;
	string num2;
	string num3;
	string num4;
	stringstream ss1;
	stringstream ss2;
	stringstream ss3;
	stringstream ss4;
	string dataclass = "ruvEP";

	ss1 << _numm; ss1 >> num1;
	ss2 << STEP_TIME; ss2 >> num2;
	ss4 << inletu; ss4 >> num4;
	title = this->mshAddress + "result/" + mshName + "_TIME" + num2 + "_u" + num4 + "_t" + num1 + ".cellPUVT";

	fstream writeFile;
	writeFile.open(title, ios::out | ios::trunc);
	for (int iCell = 0; iCell < mesh.getNumOfCells(); iCell++){
		for (int iDimu = 0; iDimu < DIMU; iDimu++){
			writeFile << fixed << setprecision(accuracyFile) << setw(accuracyFile) << allData.getPNewCellPUVT()[iCell][iDimu] << " ";
		}
		writeFile << endl;
	}



}



