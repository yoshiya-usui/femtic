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
#ifndef DBLDEF_DOUBLE_SPARSE_SQUARE_MATRIX
#define DBLDEF_DOUBLE_SPARSE_SQUARE_MATRIX

#include <string>
#include "PARDISOSolverDouble.h"
#include "DoubleSparseMatrix.h"

class DoubleSparseSquareMatrix : public DoubleSparseMatrix{

public:

	//Default Constructer
	explicit DoubleSparseSquareMatrix();

	//Constructer
	explicit DoubleSparseSquareMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~DoubleSparseSquareMatrix();

	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Set Degree of equation
	void setDegreeOfEquation( const int nEq );

	//Initialize matrix and right-hand side vectors
	void initializeMatrixAndRhsVectors( const int nEq, const int nRhs );

	//Initialize matrix solver
	virtual void initializeMatrixSolver( const std::string& oocHeaderName, const int imode ) = 0;

	//Anaysis phase of matrix solver
	void analysisPhaseMatrixSolver();

	//Numerical factorization phase of matrix solver
	void factorizationPhaseMatrixSolver();

	//Solve phase of matrix solver with a specified number of right-hand side
	void solvePhaseMatrixSolver( double* solution, const long long iRhsStart ,const int nRhs );

	//Solve phase of matrix solver
	void solvePhaseMatrixSolver( double* solution );

	//Solve phase of matrix solver
	void solvePhaseMatrixSolver( const int nrhs, double* rhs, double* solution );

	//Solve phase of matrix solver by the conjugate gradient method with the point Jacobi preconditioner
	//@note Matrix should be symmetric
	void solvePhaseMatrixSolverByPCGPointJacobi(const int nrhs, double* rhs, double* solution) const;

	//Release memory of matrix solver
	void releaseMemoryMatrixSolver();

	//Get memory required by matrix solver
	void writeMemoryRequiredByMatrixSolver();

	//Release memory
	virtual void releaseMemory();

	//Get Degree of equation
	int getDegreeOfEquation() const;

protected:
	//PARDISO solver
	PARDISOSolverDouble m_pardisoSolver;

private:
	//Copy constructer
	DoubleSparseSquareMatrix(const DoubleSparseSquareMatrix &matrix );

	// Assignment operator
	DoubleSparseSquareMatrix& operator=(const DoubleSparseSquareMatrix& rhs);

};

#endif
