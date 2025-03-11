#include "solverRoe.h"

void solverRoe::set(VCTDIMU _LeftQ, VCTDIMU _RightQ){
	leftQ = _LeftQ;
	rightQ = _RightQ;

	D = sqrt(abs(rightQ(0) / leftQ(0)));//abs may cause problem
	roe_rho = sqrt(rightQ(0) * leftQ(0)); 

	for (int iDim = 0; iDim < DIM; iDim++){
		roeU(iDim) = (leftQ(iDim+1) / leftQ(0) + D * rightQ(iDim + 1) / rightQ(0)) / (1. + D);
	}
	roe_ht = (getht(leftQ) + D * getht(rightQ)) / (1. + D);
	roe_a = sqrt(abs((GAMMA - 1.) * (roe_ht - 0.5 * (roeU.norm()*roeU.norm()))));//there may be problem when abs() is used
	
}

VCTDIMU  solverRoe::solverx(){
	//calculate roe flux on face on X	
	Eigen::Matrix<NUM, DIMU, DIMU> roe_A_eigenleft = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	Eigen::Matrix<NUM, DIMU, DIMU> roe_absA_eigen = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	Eigen::Matrix<NUM, DIMU, DIMU> roe_absA = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	roe_absA_eigen(0, 0) = entropyRefine(abs(this->roeU(0) - this->roe_a));
	roe_absA_eigen(1, 1) = entropyRefine(abs(this->roeU(0)));
	roe_absA_eigen(2, 2) = entropyRefine(abs(this->roeU(0)));
	roe_absA_eigen(3, 3) = entropyRefine(abs(this->roeU(0) + this->roe_a));
	roe_A_eigenleft << 1, 1, 0, 1, roeU(0) - roe_a, roeU(0), 0, roeU(0) + roe_a, roeU(1), roeU(1), 1, roeU(1), roe_ht - roeU(0) * roe_a, 0.5 * (roeU(0) * roeU(0) + roeU(1) * roeU(1)), roeU(1), roe_ht + roeU(0) * roe_a;
	Eigen::Matrix<NUM, DIMU, DIMU> roe_A_eigenright = roe_A_eigenleft.inverse();
	roe_absA = roe_A_eigenleft * roe_absA_eigen * roe_A_eigenright;
	
	VCTDIMU F_left;
	VCTDIMU F_right;
	VCTDIMU delta_U;

	delta_U = rightQ - leftQ;
	F_left << leftQ(1), leftQ(1)* leftQ(1) / (leftQ(0) + EOR) + getP(leftQ), leftQ(1)* leftQ(2) / (leftQ(0) + EOR), getht(leftQ)* leftQ(1);
	F_right << rightQ(1), rightQ(1)* rightQ(1) / (rightQ(0) + EOR) + getP(rightQ), rightQ(1)* rightQ(2) / (rightQ(0) + EOR), getht(rightQ)* rightQ(1);
	//all the solution of roe X direct
	VCTDIMU roe_fcx;
	roe_fcx = 0.5 * (F_left + F_right) - 0.5 * roe_absA * delta_U;
	return roe_fcx;
}
VCTDIMU  solverRoe::solvery(){
	//calculate roe flux on face on Y
	Eigen::Matrix<NUM, DIMU, DIMU> roe_A_eigenleft = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	Eigen::Matrix<NUM, DIMU, DIMU> roe_absA_eigen = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	Eigen::Matrix<NUM, DIMU, DIMU> roe_absA = Eigen::Matrix<NUM, DIMU, DIMU>::Zero();
	roe_absA_eigen(0, 0) = entropyRefine(abs(this->roeU(1) - this->roe_a));
	roe_absA_eigen(1, 1) = entropyRefine(abs(this->roeU(1)));
	roe_absA_eigen(2, 2) = entropyRefine(abs(this->roeU(1)));
	roe_absA_eigen(3, 3) = entropyRefine(abs(this->roeU(1) + this->roe_a));
	roe_A_eigenleft << 1, 1, 0, 1, roeU(0), roeU(0), 1, roeU(0), roeU(1) - roe_a, roeU(1), 0, roeU(1) + roe_a, roe_ht - roeU(1) * roe_a, 0.5 * (roeU(0) * roeU(0) + roeU(1) * roeU(1)), roeU(0), roe_ht + roeU(1) * roe_a;
	Eigen::Matrix<NUM, DIMU, DIMU> roe_A_eigenright = roe_A_eigenleft.inverse();
	roe_absA = roe_A_eigenleft * roe_absA_eigen * roe_A_eigenright;

	VCTDIMU F_left;
	VCTDIMU F_right;
	VCTDIMU delta_U;

	delta_U = rightQ - leftQ; 
	F_left << leftQ(2), leftQ(1)* leftQ(2) / (leftQ(0)), leftQ(2)* leftQ(2) / (leftQ(0)) + getP(leftQ), getht(leftQ)* leftQ(2);
	F_right << rightQ(2), rightQ(1)* rightQ(2) / rightQ(0), rightQ(2)* rightQ(2) / (rightQ(0)) + getP(rightQ), getht(rightQ)* rightQ(2);
	//all the solution of roe Y direct
	VCTDIMU roe_fcy;
	VCTDIMU roeUy;
	roe_fcy = 0.5 * (F_left + F_right) - 0.5 * roe_absA * delta_U;

	return roe_fcy;
}
VCTDIMU  solverRoe::solverAll(int _xyz) {
	//calculate roe flux on face on All 
	MTDIMU_DIMU roe_A_eigenleft = MTDIMU_DIMU::Zero();
	MTDIMU_DIMU roe_absA_eigen = MTDIMU_DIMU::Zero();
	MTDIMU_DIMU roe_absA = MTDIMU_DIMU::Zero();
	roe_absA_eigen(0, 0) = entropyRefine(abs(roeU(_xyz) - roe_a));
	for (int iDim = 0; iDim < DIM; iDim++) {
		roe_absA_eigen(iDim+1, iDim+1) = entropyRefine(abs(roeU(_xyz)));
	}
	roe_absA_eigen(DIM + 1, DIM + 1) = entropyRefine(abs(roeU(_xyz) + roe_a));
	VCTDIMU F_left;
	VCTDIMU F_right;
	VCTDIMU delta_U;

#if DIM == 2
	if (_xyz == 0){//x
		roe_A_eigenleft << 1, 1, 0, 1, roeU(0) - roe_a, roeU(0), 0, roeU(0) + roe_a, roeU(1), roeU(1), 1, roeU(1), roe_ht - roeU(0) * roe_a, 0.5 * (roeU(0) * roeU(0) + roeU(1) * roeU(1)), roeU(1), roe_ht + roeU(0) * roe_a;
		F_left << leftQ(1), leftQ(1)* leftQ(1) / (leftQ(0) + EOR) + getP(leftQ), leftQ(1)* leftQ(2) / (leftQ(0) + EOR), getht(leftQ)* leftQ(1);
		F_right << rightQ(1), rightQ(1)* rightQ(1) / (rightQ(0) + EOR) + getP(rightQ), rightQ(1)* rightQ(2) / (rightQ(0) + EOR), getht(rightQ)* rightQ(1);

	}
	else {//y
		roe_A_eigenleft << 1, 1, 0, 1, roeU(0), roeU(0), 1, roeU(0) , roeU(1)- roe_a, roeU(1), 0, roeU(1)+ roe_a, roe_ht - roeU(1) * roe_a, 0.5 * (roeU(0) * roeU(0) + roeU(1) * roeU(1)), roeU(0), roe_ht + roeU(1) * roe_a;
		F_left << leftQ(2), leftQ(1)* leftQ(2) / (leftQ(0) + EOR), leftQ(2)* leftQ(2) / (leftQ(0) + EOR) + getP(leftQ), getht(leftQ)* leftQ(2);
		F_right << rightQ(2), rightQ(1)* rightQ(2) / (rightQ(0) + EOR), rightQ(2)* rightQ(2) / (rightQ(0) + EOR) + getP(rightQ), getht(rightQ)* rightQ(2);

	}

#else//DIM == 3


#endif //all roe calculation 

	MTDIMU_DIMU roe_A_eigenright = roe_A_eigenleft.inverse();
	roe_absA = roe_A_eigenleft * roe_absA_eigen * roe_A_eigenright;

	delta_U = rightQ - leftQ;
	//all the solution of roe X direct
	VCTDIMU roe_fcx;
	roe_fcx = 0.5 * (F_left + F_right) - 0.5 * roe_absA * delta_U;
	return roe_fcx;
}


NUM solverRoe::entropyRefine(NUM _x) {
	NUM entropyError = 0.125;
	if (_x > entropyError) {
		return _x;
	}
	else {
		return (_x * _x + entropyError * entropyError) / 2 / entropyError;
	}

}
