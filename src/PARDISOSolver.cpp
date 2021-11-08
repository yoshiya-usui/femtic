//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//-------------------------------------------------------------------------------------------------------
#include <iostream>
#include <stdio.h>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"

#include "PARDISOSolver.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"

// Default constructer
PARDISOSolver::PARDISOSolver():
	m_maxfct(1),
	m_mnum(1),
	m_mtype(PARDISOSolver::COMPLEX_AND_SYMMETRIC_MATRIX),
	m_numEquations(NULL),
	m_msglvl(0),
	m_peakMemorySymbolicFactorization(0.0),
	m_permanetMemorySymbolicFactorization(0.0),
	m_memoryForNumericalFactorizationIncore(0.0),
	m_memoryForNumericalFactorizationOutcore(0.0),
	m_solutionStage(PARDISOSolver::MEMORY_RELEASED)
{
	for( int i = 0; i < 64; ++i ){
      m_pt[i] = NULL; // Initialize
	  m_iparm[i] = NULL; // Initialize
    }
}

// Default constructer
PARDISOSolver::PARDISOSolver( const long long int matrixType ):
	m_maxfct(1),
	m_mnum(1),
	m_mtype(matrixType),
	m_numEquations(NULL),
	m_msglvl(0),
	m_peakMemorySymbolicFactorization(0.0),
	m_permanetMemorySymbolicFactorization(0.0),
	m_memoryForNumericalFactorizationIncore(0.0),
	m_memoryForNumericalFactorizationOutcore(0.0),
	m_solutionStage(PARDISOSolver::MEMORY_RELEASED)
{
	for( int i = 0; i < 64; ++i ){
      m_pt[i] = NULL; // Initialize
	  m_iparm[i] = NULL; // Initialize
    }
}

// Destructer
PARDISOSolver::~PARDISOSolver(){

	if( m_solutionStage != PARDISOSolver::MEMORY_RELEASED ){ // Release memory of PARDISO solver
		releaseMemory();
	}

}

// Initialize PARDISO solver
void PARDISOSolver::initialize( const std::string& oocHeaderName, const long long int imode, const long long int type ){

	//if( numThreads < 0 ){
	//	OutputFiles::m_logFile << "Error : Total number of threads must be greater than or equals to 1 !! numThreads = " << numThreads << std::endl;
	//	exit(1);
	//}

	m_mtype = type;

	if( m_solutionStage != PARDISOSolver::MEMORY_RELEASED ){ // Release memory of PARDISO solver
		releaseMemory();
	}

	for( int i = 0; i < 64; ++i ){
		m_iparm[i] = NULL; // Initialize
    }

	pardisoinit( m_pt, &m_mtype, m_iparm );

	m_iparm[0] = 1; // Do not use default parameters
	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	const int numThreads = pAnalysisControl->getNumThreads();
	if( numThreads == 1 ){
		m_iparm[1] = 2; // METIS
	}else{
		m_iparm[1] = 3; // Parallel version of the nested dissection algorithm
	}
	m_iparm[3] = 0; // Do not perform preconditioned CGS/CG iterations
	m_iparm[4] = 0; // Do not use user permutation vector
	m_iparm[5] = 0; // Solution vector is returned to array x
	m_iparm[7] = 0; // Two steps of iterative refinements if pivots are perturbed at the numerical factorization stage
	m_iparm[9] = 8; // Small pivots are perturbed with eps = 10^(-8) ( default value of symmetric indefinite matrices )
	m_iparm[10] = 0; // Do not perform scaling
	m_iparm[11] = 0; // Solve normally Ax=b
	m_iparm[12] = 0; // Do not perform scaling symmetric weighted matching
	m_iparm[17] = 0; // Do not report the number of non-zero elements in the factors
	m_iparm[18] = 0; // Do not report Mflops that are necessary to factor the matrix A
	m_iparm[20] = 1; // 1x1 and 2x2 Bunch and Kaufman pivoting during the factorization stage
	m_iparm[23] = 0; // Use 1x1 and 2x2 Bunch and Kaufman pivoting during the factorization stage
	m_iparm[24] = 0; // Use the parallel algorithm for solve step
	m_iparm[26] = 0; // Do not check the sparse matrix representation
	m_iparm[27] = 0; // Double precision 
	m_iparm[30] = 0; // Do not assume sparse right-hand sides and sparse solution as sparce
	m_iparm[34] = 1; // Zero-based indexing
	//m_iparm[59] = 0; // In-core PARDISO
	m_iparm[59] = imode; // PARDISO mode

	//// Specifies the number of threads to use
	//mkl_set_num_threads( numThreads );

#ifdef _INTEL_LT_21
	// Set the header name of out-of-core files
	const PARDISO_ENV_PARAM param = PARDISO_OOC_FILE_NAME;
	pardiso_setenv( m_pt, &param, oocHeaderName.c_str() );
#endif	
	m_solutionStage = PARDISOSolver::INITIALIZED;

}

// Analysis phase of PARDISO solver
void PARDISOSolver::analysis( const long long int nEq, long long int* rowIndex, long long int* columns ){

	if( m_solutionStage < PARDISOSolver::INITIALIZED ){
		OutputFiles::m_logFile << "Error : Forward solver has not been initialized yet." << std::endl;
		exit(1);
	}else if( m_solutionStage == PARDISOSolver::ANALYZED ){
		OutputFiles::m_logFile << "Warning : Analysis phase has already been performed." << std::endl;
	}

	m_numEquations = nEq;

	long long int phase = 11;
	double ddum;
	long long int idum;
	long long int error;
	long long int nrhs = 1;

	pardiso_64( m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_numEquations, &ddum, rowIndex, columns,
		&idum, &nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &error);
	if (error != 0)
    {
 	//	OutputFiles::m_logFile << "Error : Error during analysis phase of forward solver. : error = " << error << std::endl;
		//exit(1);
		outputErrorMessages( error );
    }

	m_peakMemorySymbolicFactorization = m_iparm[14];
	m_permanetMemorySymbolicFactorization = m_iparm[15];
	m_memoryForNumericalFactorizationIncore = m_iparm[16];
	m_memoryForNumericalFactorizationOutcore = m_iparm[62];

	m_solutionStage = PARDISOSolver::ANALYZED;

}

// Release memory of PARDISO solver
void PARDISOSolver::releaseMemory(){

	//if( m_solutionStage == PARDISO::MEMORY_RELEASED ){
	//	OutputFiles::m_logFile << "Warning : Memory has already been released." << std::endl;
	//	return;
	//}
	if( m_solutionStage < PARDISOSolver::ANALYZED ){
		return;
	}

	long long int phase = -1;
	long long int idum;
	double ddum;
	long long int error;
	long long int nrhs = 1;
	pardiso_64( m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_numEquations, &ddum, &idum, &idum,
		&idum, &nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &error);
	if (error != 0)
    {
 	//	OutputFiles::m_logFile << "Error : Error during memory release phase of forward solver. : error = " << error << std::endl;
		//exit(1);
 		outputErrorMessages( error );
   }

	m_solutionStage = PARDISOSolver::MEMORY_RELEASED;

}

// Get memory required by PARDISO solver
void PARDISOSolver::writeMemoryRequired() const{

	if( m_solutionStage < PARDISOSolver::ANALYZED ){
		OutputFiles::m_logFile << "Error : Memory required by forward solver is obtained at the analysis phase : stage = " << m_solutionStage << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	OutputFiles::m_logFile << "# Total peak memory required during the analysis and symbolic factorization phase    : " << m_peakMemorySymbolicFactorization << " [KByte]" << std::endl;
	OutputFiles::m_logFile << "# Permanent memory required from the analysis phase to the solve phases              : " << m_permanetMemorySymbolicFactorization << " [KByte]" << std::endl;
	OutputFiles::m_logFile << "# Total memory consumed by in-core forward solver for internal float point arrays    : " << m_memoryForNumericalFactorizationIncore << " [KByte]" << std::endl;
	OutputFiles::m_logFile << "# Minimum memory consumed by out-core forward solver for internal float point arrays : " << m_memoryForNumericalFactorizationOutcore << " [KByte]" << std::endl;
	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

}

// Get stage of PARDISO solver
int PARDISOSolver::getSolutionStage() const{
	return m_solutionStage;
}

// Set stage of PARDISOSolver solver
void PARDISOSolver::setSolutionStage( const int stage ){

	if( stage < MEMORY_RELEASED || stage > SOLVED ){
		OutputFiles::m_logFile << "Error : Stage number is wrong. : stage = " << stage << std::endl;
		exit(1);		
	}

	m_solutionStage = stage;

}

// Output error messages
void PARDISOSolver::outputErrorMessages( const int ier ) const{

	switch (ier){
		case -1:
			OutputFiles::m_logFile << "Error : Some parameters passed to forward solver may be wrong." << std::endl;
			break;
		case -2:
			OutputFiles::m_logFile << "Error : Insufficient memory for forward solver." << std::endl;
			break;
		case -3:
			OutputFiles::m_logFile << "Error : Some problems occur in reordering." << std::endl;
			break;
		case -4:
			OutputFiles::m_logFile << "Error : Zero pivot is found." << std::endl;
			break;
		case -5:
			OutputFiles::m_logFile << "Error : Internal error of forward solver." << std::endl;
			break;
		case -6:
			OutputFiles::m_logFile << "Error : Reordering failed." << std::endl;
			break;
		case -7:
			OutputFiles::m_logFile << "Error : Matrix is singular." << std::endl;
			break;
		case -8:
			OutputFiles::m_logFile << "Error : 32-bit integer overflow problem." << std::endl;
			break;
		case -9:
			OutputFiles::m_logFile << "Error : Insufficient memory for out-of-core mode of forward solver." << std::endl;
			break;
		case -10:
			OutputFiles::m_logFile << "Error : Fail to open out-of-core file of forward solver." << std::endl;
			break;
		case -11:
			OutputFiles::m_logFile << "Error : Fail to read/write out-of-core file of forward solver." << std::endl;
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown error of forward solver. ier = " << ier << std::endl;
			break;
	}

	exit(1);

}
