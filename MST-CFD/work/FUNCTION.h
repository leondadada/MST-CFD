#pragma once 
#include "../include/CONST.h"
#include <Eigen/Eigen>
#include <math.h>
#include <string>
#include <iostream>

using namespace Eigen;
using namespace std;

NUM getht(VCTDIMU _mx);
NUM getP(VCTDIMU _mx);
NUM getT(VCTDIMU _mx);
NUM getMa(VCTDIMU);
NUM getPBasedMa(VCTDIMU);
NUM getPBasedRho(VCTDIMU);
VCTDIMU getFluxFFromQ(VCTDIMU);
VCTDIMU getFluxGFromQ(VCTDIMU);
int hexStringToInt(string);
string betweenQuotes(string);
string betweenSecBrackets(string);
NUM stodd(string );


