#pragma once

#define DIM 2 //dimension 
#define DIMU DIM + 2
#define FELNUM 4 // triangle grid or quad?
#define ACCURACY 1 // 1 means 1st-order upwind//2 means 2order upwind 
#define NUM_CPU_THREADS 8
//#define TIMESOLVER PSolver//PSolver// RhoSolver
#define RHO_P 1 //0 means rho 
#define RHOSOLVER SolverAusm
#define PSOLVER SolverSIMPLE //SolverCoupled
#define SPARSESOLVER SimplicialLDLT

#define FLAGVISCID 0 // 0 means no viscid
#define FLAGPSUEDO 0 // 0 means no psuedotime
#define FLAG_WALLSMOOTH 1 //0 means u==0
#define FLAG_PBASED_EXPLICIT 1

//BiCGSTAB 
//ConjugateGradient 
//SimplicialLDLT  
//SimplicialCholesky 


#define READ_PLT 0 //which time to start?


#define NUM double // float 
#define VCTDIMU Matrix<NUM, DIMU, 1>
#define VCTDIM Matrix<NUM, DIM, 1>
#define VCTDIM1 Matrix<NUM, DIM+1, 1>
#define MTDIM1 Matrix<NUM, DIM+1, DIM+1>
#define MTDIM_DIMU Matrix<NUM, DIM, DIMU>
#define MTDIMU_DIM Matrix<NUM, DIMU, DIM>
#define MTDIM_DIM Matrix<NUM, DIM, DIM>
#define MTDIMU_DIMU Matrix<NUM, DIMU, DIMU>

#define CP 1002.12
#define CV 715.8
#define RGAS CP-CV 
#define GAMMA 1.4
#define EOR 1e-10
#define EOR2 1e-7
#define accuracyFile 15
#define INF 1e8
#define VISCIDMU 1.7894e-05
#define VISCIDLAMBDA -0.666667*VISCIDMU
#define TEMPK 0.0242


#define TIME  8e-1 //total time 
#define STEP_TIME  4e+3 //time steps per s
#define SAVE_TIME  2.5e-3 /*save per SAVE_TIMEs */
#define GRAD_INTERVAL 1
#define PT_STEP 5 // peusdo time steps 
#define PT_STEP_TIME   PT_STEP
#define LU_INTERVAL 5 // must larger than 1,at least 1.
#define SIMPLE_INTERVAL 5
#define COUPLED_INTERVAL 10
#define SIMPLE_ALPHA_P  0.3  // cof of p =p* + alpha*p'
#define SIMPLE_ALPHA_U SIMPLE_ALPHA_P // cof of p =p* + alpha*p'





//----------------------------------------
 //below is unscaled condition 

#define iniT  1/286.32
#define inirho  1 // p = 1
#define iniu 0   // 340 : ma = 1
#define iniv  0
#define iniw 0
#define iniE  inirho*(iniT*CV+0.5*(iniu*iniu+iniv*iniv)) //consider


#define inletT  iniT
#define inletrho  inirho // 340 ma = 1
#define inletu  iniu
#define inletv  iniv
#define inletw  iniw
#define inletE  inletrho*(inletT*CV+0.5*(inletu*inletu+inletv*inletv)) //consider inlet P =0


