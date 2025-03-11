#include "./FUNCTION.h"

NUM getht(VCTDIMU _mx) {
	NUM p = (_mx(DIMU-1) - 0.5 * (_mx(1) * _mx(1) + _mx(2) * _mx(2)) / _mx(0)) * (GAMMA - 1);
	NUM ht = (_mx(DIMU-1) + p) / _mx(0) ;
	return ht;
}
NUM getT(VCTDIMU _mx) {
	NUM T = (_mx(DIMU - 1) - 0.5 * (_mx(1) * _mx(1) + _mx(2) * _mx(2)) / _mx(0)) / _mx(0) / CV;
	return T;
}
NUM getP(VCTDIMU _mx) {
	NUM p = (_mx(DIMU - 1) - 0.5 *  (_mx(1) * _mx(1) + _mx(2) * _mx(2)) / _mx(0))* (GAMMA - 1);
	return p;
}
NUM getMa(VCTDIMU _mx) {
	NUM p = (_mx(DIMU - 1) - 0.5 *  (_mx(1) * _mx(1) + _mx(2) * _mx(2)) / _mx(0))* (GAMMA - 1);
	NUM ma = sqrt((_mx(1)*_mx(1)+ _mx(2)*_mx(2)) /(GAMMA*p*_mx(0)));
	return ma;
}
NUM getPBasedMa(VCTDIMU _puvt) {
	NUM ma = sqrt((_puvt(1)*_puvt(1) + _puvt(2)*_puvt(2)) *getPBasedRho(_puvt) / (GAMMA*_puvt[0]));
	return ma;
}
NUM getPBasedRho(VCTDIMU _PUVT) {
	return _PUVT(0) / (CP-CV) / _PUVT(DIMU-1);
}

VCTDIMU getFluxFFromQ(VCTDIMU _Q) {
	VCTDIMU fluxF;
	fluxF << _Q(1), _Q(1)*_Q(1) / _Q(0) + getP(_Q), _Q(1)*_Q(2) / _Q(0), _Q(1) / _Q(0)*(_Q(DIMU - 1) + getP(_Q));
	return fluxF;
}
VCTDIMU getFluxGFromQ(VCTDIMU _Q) {
	VCTDIMU fluxG;
	fluxG << _Q(2), _Q(2)*_Q(1) / _Q(0), _Q(2)*_Q(2) / _Q(0) + getP(_Q), _Q(2) / _Q(0)*(_Q(DIMU - 1) + getP(_Q));
	return fluxG;
}


int hexStringToInt(string _hexString) {
	int numInt = 0;
	for (int i = _hexString.size(); i>0; i--) {
		int bit = _hexString.size() - i;
		if (_hexString[bit] >= 'a' && _hexString[bit] <= 'f') {
			numInt += round(pow(16, i - 1)) *(_hexString[bit] - 'a' + 10);
		}
		else {
			if (_hexString[bit] >= '0' && _hexString[bit] <= '9') {
				numInt += round(pow(16, i - 1)) * (_hexString[bit] - '0');
			}
		}
	}
	return numInt;
}
string betweenQuotes(string _stringLine) {
	string stringBetween;
	stringstream stringstream;
	stringstream << _stringLine;
	getline(stringstream, stringBetween, '\"');
	getline(stringstream, stringBetween, '\"');
	return stringBetween;
}//only one quotes,return inside quotes
string betweenSecBrackets(string _stringLine) {
	string stringBetween;
	stringstream stringstream;
	stringstream << _stringLine;
	getline(stringstream, stringBetween, '(');
	getline(stringstream, stringBetween, '(');
	getline(stringstream, stringBetween, ')');
	return stringBetween;
}//only one Brackets,return inside Brackets
NUM stodd(string _NUMnum) {
	string c = _NUMnum;
	NUM cc = 0;
	cc = (NUM)atof(c.c_str());
	return cc;
}
