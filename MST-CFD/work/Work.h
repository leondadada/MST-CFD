#pragma once 
#include "../include/CONST.h"
#include "../mesh/MshBlock.h"
#include "./FUNCTION.h"
#include "../time/Time.h"
#include "../data/AllData.h"
#include <iomanip>

using namespace std;

class Work 
{
public:
	~Work();
	void work(string _mshName);
private:
	int t = READ_PLT;
	AllData allData;
	fstream fLog;
	MshBlock mesh;
	void writedataNode(MshBlock*, NUM, int);
	void writedataCell(MshBlock*, VCTDIMU*, int, int);
	void writedataCellValue(MshBlock*, VCTDIMU*, int, int);
	void writedataCellValueAll(MshBlock *, VCTDIMU *, int);
	void writedataPBasedCellValueAll(MshBlock * _pMesh, VCTDIMU * _p1PUVT, int _time_now);
	void writedataPBasedMshNodePlt(MshBlock * _pMesh, VCTDIMU * _p1PUVT, int _time_now);
	void writedataRhoBasedMshNodePlt(MshBlock * _pMesh, VCTDIMU * _p1PUVT, int _time_now);
	void writedataRhoBasedMshNodeQUADPlt(MshBlock * _pMesh, VCTDIMU * _p1CellQ, int _time_now);
	void readPbasedNodePlt(string _mshName);
	void readPbasedCellPlt(string);
	void writedataCellContinue(string _mshName, int _numm);
	string mshAddress;
	string mshName;
	


};



