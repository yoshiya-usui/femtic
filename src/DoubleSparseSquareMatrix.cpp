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
#include <stddef.h> // For null pointer
#include <stdlib.h> // For exit
#include <iostream>
#include "DoubleSparseSquareMatrix.h"
#include "OutputFiles.h"
#include <assert.h>
#include <math.h>

#ifdef _DEBUG_WRITE_FOR_BOTTOM_RESISTIVITY
#ifdef _LINUX
#include <sys/time.h>
#include <sys/resource.h>
#endif
#endif

//Default Constructer
DoubleSparseSquareMatrix::DoubleSparseSquareMatrix():
	DoubleSparseMatrix()
{}

// Constructer
DoubleSparseSquareMatrix::DoubleSparseSquareMatrix( const int nEq, const int nRhs ):
	DoubleSparseMatrix( nEq, nEq, nRhs )
{
	//if( nEq <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of equation specified is less than or equals to zero. : nEq = " << nEq << std::endl;
	//	exit(1);		
	//}
	//
	//if( nRhs <= 0 ){
	//	OutputFiles::m_logFile << "Error : Total number of right-hand side vectors is specified to be less than or equals to zero. : nRhs = " << nRhs << std::endl;
	//	exit(1);		
	//}
	assert( nEq > 0 );
	assert( nRhs > 0 );
}

// Destructer
DoubleSparseSquareMatrix::~DoubleSparseSquareMatrix(){
	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}
}

// Set number of rows and columns
void DoubleSparseSquareMatrix::setNumRowsAndColumns( const int nrows, const int ncols ){

	//if( nrows != ncols ){
	//	OutputFiles::m_logFile << "Error : Number of rows and the one of columns are different for square matrix. : nrows = " << nrows << ", ncols = " << ncols << std::endl;
	//	exit(1);		
	//}
	assert( nrows == ncols );

	DoubleSparseMatrix::setNumRowsAndColumns( nrows, ncols );

}

// Set Degree of equation
// Note : This function must be called BEFORE the matrix is converted into CRS format
void DoubleSparseSquareMatrix::setDegreeOfEquation( const int nEq ){

	//if( nEq <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of equation specified is less than or equals to zero. : nEq = " << nEq << std::endl;
	//	exit(1);		
	//}
	assert( nEq > 0 );

	setNumRowsAndColumns( nEq, nEq );
}

//Initialize matrix and right-hand side vectors
void DoubleSparseSquareMatrix::initializeMatrixAndRhsVectors( const int nEq, const int nRhs ){

	//if( nEq <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of equation is specified to be less than or equals to zero. : nEq = " << nEq << std::endl;
	//	exit(1);		
	//}
	//
	//if( nRhs <= 0 ){
	//	OutputFiles::m_logFile << "Error : Total number of right-hand side vectors is specified to be less than or equals to zero. : nRhs = " << nRhs << std::endl;
	//	exit(1);		
	//}
	assert( nEq > 0 );
	assert( nRhs > 0 );

	releaseMemoryMatrixSolver();

	DoubleSparseMatrix::initializeMatrixAndRhsVectors( nEq, nEq, nRhs );

}

//Anaysis phase of matrix solver
void DoubleSparseSquareMatrix::analysisPhaseMatrixSolver(){
	assert( m_hasConvertedToCRSFormat );
	m_pardisoSolver.analysis( m_numRows, m_rowIndex, m_columns );
}

//Numerical factorization phase of matrix solver
void DoubleSparseSquareMatrix::factorizationPhaseMatrixSolver(){
	assert( m_hasConvertedToCRSFormat );
	m_pardisoSolver.numericalFactorization( m_rowIndex, m_columns, m_values );
}

//Solve phase of matrix solver with a specified number of right-hand side
void DoubleSparseSquareMatrix::solvePhaseMatrixSolver( double* solution, const long long iRhsStart ,const int nRhs ){
	assert( m_hasConvertedToCRSFormat );
	const long long index = static_cast<long long>(m_numRows) * iRhsStart;
	m_pardisoSolver.solve( m_rowIndex, m_columns, m_values, nRhs, &m_rightHandSideVector[index], solution );
}

//Solve phase of matrix solver
void DoubleSparseSquareMatrix::solvePhaseMatrixSolver( double* solution ){
	assert( m_hasConvertedToCRSFormat );
	m_pardisoSolver.solve( m_rowIndex, m_columns, m_values, m_numRightHandSideVectors, m_rightHandSideVector, solution );
}

//Solve phase of matrix solver
void DoubleSparseSquareMatrix::solvePhaseMatrixSolver( const int nrhs, double* rhs, double* solution ){
	assert( m_hasConvertedToCRSFormat );
	m_pardisoSolver.solve( m_rowIndex, m_columns, m_values, nrhs, rhs, solution );
}

//Solve phase of matrix solver by the conjugate gradient method with the point Jacobi preconditioner
//@note Matrix should be symmetric
void DoubleSparseSquareMatrix::solvePhaseMatrixSolverByPCGPointJacobi(const int nrhs, double* rhs, double* solution) const{
	assert(m_hasConvertedToCRSFormat);

	const int maxIterationNumber = m_numRows;
	const double eps = 1.0e-20;
	double* invDiagonals = new double[m_numRows];
	double* workP = new double[m_numRows];
	double* workR = new double[m_numRows];// Residuals
	double* workQ = new double[m_numRows];
	double* workX = new double[m_numRows];// Solution vector
	double* workZ = new double[m_numRows];

	for (int irow = 0; irow < m_numRows; ++irow)
	{
		for (int j = m_rowIndex[irow]; j < m_rowIndex[irow + 1]; ++j)
		{
			if (irow == m_columns[j])
			{
				invDiagonals[irow] = 1.0 / m_values[j];
			}
		}
	}
	for (int irhs = 0; irhs < nrhs; ++irhs)
	{
		// Initial solution is a zero vector
		for (int irow = 0; irow < m_numRows; ++irow)
		{
			workX[irow] = 0.0;
		}
		// [r0] = [b] - [A][x0]
		double normOfRhsVector(0.0);
		for (int irow = 0; irow < m_numRows; ++irow)
		{
			const long long int index = static_cast<long long int>(irow) + static_cast<long long int>(irhs) * static_cast<long long int>(m_numRows);
			normOfRhsVector += rhs[index] * rhs[index];
			workR[irow] = rhs[index];
		}
		int iter = 0;
		double rhoPre(0.0);
		for (; iter < maxIterationNumber; ++iter)
		{
			// [z] = [M]^-1[r]
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				workZ[irow] = invDiagonals[irow] * workR[irow];
			}
			// rho = [r]T[z]
			double rho(0.0);
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				rho += workR[irow] * workZ[irow];
			}
			if (iter == 0)
			{
				// [p0] - [z0]
				for (int irow = 0; irow < m_numRows; ++irow)
				{
					workP[irow] = workZ[irow];
				}
			}
			else
			{
				// [p] = [z] + beta*[p]
				const double beta = rho / rhoPre;
				for (int irow = 0; irow < m_numRows; ++irow)
				{
					workP[irow] = workZ[irow] + beta * workP[irow];
				}
			}
			// [q] = [A][p]
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				workQ[irow] = 0.0;
				for (int j = m_rowIndex[irow]; j < m_rowIndex[irow + 1]; ++j)
				{
					workQ[irow] += m_values[j] * workP[m_columns[j]];
				}
			}
			// alpha = rho / [p]T[q]
			double pq(0.0);
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				pq += workP[irow] * workQ[irow];
			}
			const double alpha = rho / pq;
			// [x] = [x] + alpha * [p]
			// [r] = [r] - alpha * [q]
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				workX[irow] += alpha * workP[irow];
				workR[irow] -= alpha * workQ[irow];
			}
			// Check convergence
			double normOfResidualVector(0.0);
			for (int irow = 0; irow < m_numRows; ++irow)
			{
				normOfResidualVector += workR[irow] * workR[irow];
			}
			if( sqrt(normOfResidualVector/ normOfRhsVector) < eps )
			{
				break;
			}
			rhoPre = rho;
		}
		if (iter >= maxIterationNumber) {
			OutputFiles::m_logFile << "Error : PCG solver is not converged !!" << std::endl;
			exit(1);
		}
		else {
			OutputFiles::m_logFile << "# PCG solver is converged after " << iter << " iterations." << std::endl;
		}
		for (int irow = 0; irow < m_numRows; ++irow)
		{
			const long long int index = static_cast<long long int>(irow) + static_cast<long long int>(irhs) * static_cast<long long int>(m_numRows);
			solution[index] = workX[irow];
		}
	}

#ifdef _DEBUG_WRITE_FOR_BOTTOM_RESISTIVITY
#ifdef _LINUX
	{
		struct rusage r;
		if (getrusage(RUSAGE_SELF, &r) != 0) {
			/*Failure*/
		}
		OutputFiles::m_logFile << "maxrss= " << r.ru_maxrss << std::endl;
	}
#endif
#endif

	delete[] invDiagonals;
	delete[] workP;
	delete[] workR;
	delete[] workQ;
	delete[] workX;
	delete[] workZ;

}

//Release memory of matrix solver
void DoubleSparseSquareMatrix::releaseMemoryMatrixSolver(){
	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}
}

//Get memory required by matrix solver
void DoubleSparseSquareMatrix::writeMemoryRequiredByMatrixSolver(){
	m_pardisoSolver.writeMemoryRequired();
}

//Release memory
void DoubleSparseSquareMatrix::releaseMemory(){

	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}
	DoubleSparseMatrix::releaseMemory();

}

// Get Degree of equation
int DoubleSparseSquareMatrix::getDegreeOfEquation() const{
	return m_numRows;
}

//Copy constructer
DoubleSparseSquareMatrix::DoubleSparseSquareMatrix(const DoubleSparseSquareMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class DoubleSparseSquareMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
DoubleSparseSquareMatrix& DoubleSparseSquareMatrix::operator=(const DoubleSparseSquareMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class DoubleSparseSquareMatrix is not implemented." << std::endl;
	exit(1);
}
