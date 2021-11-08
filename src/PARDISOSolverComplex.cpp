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

#include "PARDISOSolverComplex.h"
#include "OutputFiles.h"

// Default constructer
PARDISOSolverComplex::PARDISOSolverComplex():
	PARDISOSolver(PARDISOSolver::COMPLEX_AND_SYMMETRIC_MATRIX)
{
}

// Constructer
PARDISOSolverComplex::PARDISOSolverComplex(const long long int matrixType):
	PARDISOSolver(matrixType)
{
}

// Destructer
PARDISOSolverComplex::~PARDISOSolverComplex(){
}

// Numerical factorization phase of PARDISO solver
void PARDISOSolverComplex::numericalFactorization( long long int* rowIndex, long long int* columns, std::complex<double>* values ){

	if( m_solutionStage < PARDISOSolver::ANALYZED ){
		OutputFiles::m_logFile << "Error : Matrix has not been analized by forward solver yet." << std::endl;
		exit(1);
	}

	long long int phase = 22;
	double ddum;
	long long int idum;
	long long int error;
	long long int nrhs = 1;
	pardiso_64( m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_numEquations, values, rowIndex, columns,
		&idum, &nrhs, m_iparm, &m_msglvl, &ddum, &ddum, &error);
	if (error != 0)
    {
 		outputErrorMessages( error );
    }

	m_solutionStage = PARDISOSolver::FACTORIZED;

}

// Solve phase of PARDISO solver
void PARDISOSolverComplex::solve( long long int* rowIndex, long long int* columns, std::complex<double>* values, long long int nrhs, std::complex<double>* rhsValues, std::complex<double>* solution ){

	if( m_solutionStage < PARDISOSolver::FACTORIZED ){
		OutputFiles::m_logFile << "Error : Matrix has not been numerically factorized by forward solver yet." << std::endl;
		exit(1);
	}

	long long int phase = 33;
	long long int idum;
	long long int error;
	pardiso_64( m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_numEquations, values, rowIndex, columns,
		&idum, &nrhs, m_iparm, &m_msglvl, rhsValues, solution, &error);
	if (error != 0)
    {
   		outputErrorMessages( error );
	}

	m_solutionStage = PARDISOSolver::SOLVED;

}

