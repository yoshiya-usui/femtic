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
#ifndef DBLDEF_COMPLEX_SPARSE_SQUARE_MATRIX
#define DBLDEF_COMPLEX_SPARSE_SQUARE_MATRIX

#include <complex>

#include "PARDISOSolverComplex.h"
#include "ComplexSparseMatrix.h"

class ComplexSparseSquareMatrix : public ComplexSparseMatrix {

public:

	//Default Constructer
	explicit ComplexSparseSquareMatrix();

	//Constructer
	explicit ComplexSparseSquareMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~ComplexSparseSquareMatrix();

	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Set Degree of equation
	void setDegreeOfEquation( const int nEq );

	////Input component of matrix.
	//void setStructureByTripletFormat( const int row, const int col );
	
	////Convert matrix from triplet format to CRS format
	//void convertToCRSFormat();

	////Add non-zero value to matrix
	//void addNonZeroValues( const int row, const int col, const std::complex<double> val );

	////Zero clear non-zero values of matrix stored by CSR format
	//void zeroClearNonZeroValues();

	////Add non-zero value to the right hand side vector
	//void addRightHandSideVector( const int row, const std::complex<double> val, const int irhs = 0 );

	////Zero clear non-zero values of the right hand side vector
	//void zeroClearRightHandSideVector();

	//Initialize matrix and right-hand side vectors
	void initializeMatrixAndRhsVectors( const int nEq, const int nRhs );

	//Initialize matrix solver
	virtual void initializeMatrixSolver( const std::string& oocHeaderName, const int imode ) = 0;

	//Anaysis phase of matrix solver
	void analysisPhaseMatrixSolver();

	//Numerical factorization phase of matrix solver
	void factorizationPhaseMatrixSolver();

	//Solve phase of matrix solver with a specified number of right-hand side
	void solvePhaseMatrixSolver( std::complex<double>* solution, const long long iRhsStart ,const int nRhs );

	//Solve phase of matrix solver
	void solvePhaseMatrixSolver( std::complex<double>* solution );

	//Release memory of matrix solver
	void releaseMemoryMatrixSolver();

	//Get memory required by matrix solver
	void writeMemoryRequiredByMatrixSolver();

	//Release memory
	virtual void releaseMemory();

	//Get Degree of equation
	int getDegreeOfEquation() const;

	////Get total number of right hand side vector 
	//int getNumRightHandSideVectors() const;

	////Return whether matrix has already been converted to CRS format
	//bool hasConvertedToCRSFormat() const;

	////Reallocate memory for right hand side vectors
	//void reallocateMemoryForRightHandSideVectors( const int nrhs );

	////Set number of right hand side vectors
	//void setNumRightHandSideVectors( const int nrhs );

	////Multiply matrix by inputed vector and add calculated vector to another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndAddResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const;

	////Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndSubtractResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const;

	////Substitute right-hand side vector to another vector
	////void substituteRhsVector( std::complex<double>* vecOut ) const;

	////Debug write the matrix componets
	//void debugWriteMatrix() const;

	////Debug write the componets of right hand side vector 
	//void debugWriteRightHandSide() const;

protected:
	//PARDISO solver
	PARDISOSolverComplex m_pardisoSolver;

private:
	//Copy constructer
	ComplexSparseSquareMatrix(const ComplexSparseSquareMatrix &matrix );

	// Assignment operator
	ComplexSparseSquareMatrix& operator=(const ComplexSparseSquareMatrix& rhs);

	////Delete the matrix of triplet ( Coordinate ) format
	//void deleteTripletMatrix();

	////Total number of equations
	//int m_numEquations;

	////Total number of non-zero elements in the matrix
	//int m_numNonZeros;

	////Total number of right hand side vectors
	//int m_numRightHandSideVectors;

	//////Rows in which non-zero compnents exist of triplet ( Coordinate ) format
	////std::vector<int> m_rowsTriplet;
	//////Columns in which non-zero compnents exist of triplet ( Coordinate ) format
	////std::vector<int> m_columnsTriplet;
	//////Values of non-zero compnents of triplet ( Coordinate ) format
	////std::vector< std::complex<double> > m_valuesTriplet;

	////Matrix stored by triplet ( Coordinate ) format
	//std::set<int>* m_matrixTripletFormat;

	////Flag indicating whether matrix has already been converted to CRS format
	//bool m_hasConvertedToCRSFormat;

	////Row indices of the compressed row storage format
	//int* m_rowIndex;

	////Columns in which non-zero compnents exist
	//int* m_columns;

	////Values of non-zero compnents
	//std::complex<double>* m_values;

	////Right hand side vector
	//std::complex<double>* m_rightHandSideVector;

};

#endif
