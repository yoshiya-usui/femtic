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
#include "ComplexSparseSquareMatrix.h"
#include "OutputFiles.h"
#include <assert.h>

//Default Constructer
ComplexSparseSquareMatrix::ComplexSparseSquareMatrix():
	ComplexSparseMatrix()
{}

// Constructer
ComplexSparseSquareMatrix::ComplexSparseSquareMatrix( const int nEq, const int nRhs ):
	ComplexSparseMatrix( nEq, nEq, nRhs )
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
ComplexSparseSquareMatrix::~ComplexSparseSquareMatrix(){

	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}

}

// Set number of rows and columns
void ComplexSparseSquareMatrix::setNumRowsAndColumns( const int nrows, const int ncols ){

	//if( nrows != ncols ){
	//	OutputFiles::m_logFile << "Error : Number of rows and the one of columns are different for square matrix. : nrows = " << nrows << ", ncols = " << ncols << std::endl;
	//	exit(1);		
	//}
	assert( nrows == ncols );

	ComplexSparseMatrix::setNumRowsAndColumns( nrows, ncols );

}

// Set Degree of equation
// Note : This function must be called BEFORE the matrix is converted into CRS format
void ComplexSparseSquareMatrix::setDegreeOfEquation( const int nEq ){

	//if( nEq <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of equation specified is less than or equals to zero. : nEq = " << nEq << std::endl;
	//	exit(1);		
	//}
	assert( nEq > 0 );

	setNumRowsAndColumns( nEq, nEq );
}
//
//// Set matrix structure ( locations of non-zero components ) by triplet format
//// Note : This function must be called BEFORE the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::setStructureByTripletFormat( const int row, const int col ){
//
//	if( m_hasConvertedToCRSFormat == true ){
//		//Matrix has already been converted to CRS format
//		OutputFiles::m_logFile <<  "Error : Matrix has already been converted to CRS format." << std::endl;
//		exit(1);
//	}
//	
//	if( row > m_numEquations - 1 || row < 0 || col > m_numEquations - 1 || col < 0 ){
//		OutputFiles::m_logFile <<  "Error : Row or column is out of the range of the matrix. row : " << row << " col : " << col << std::endl;
//		exit(1);
//	}
//
//	if( row <= col ){// Only upper triangle components are stored
//		m_matrixTripletFormat[row].insert( col );
//	}
//}
//
//// Convert matrix from triplet format to CRS format
//// Note : This function must be called BEFORE the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::convertToCRSFormat(){
//
//	if( m_hasConvertedToCRSFormat == true ){
//		//Matrix has already been converted to CRS format
//		OutputFiles::m_logFile <<  "Warning : Matrix has already been converted to CRS format." << std::endl;
//	}
//	else{
//
//		if( m_rowIndex != NULL ){
//			delete[] m_rowIndex;
//			m_rowIndex = NULL;
//		}
//		m_rowIndex = new int[ m_numEquations + 1 ];//Row indices of the compressed row storage format
//		for( int i = 0; i < m_numEquations + 1; ++i ){
//			m_rowIndex[i] = NULL; // Initialize
//		}
//
//		// Calculate total number of non-zero components and row indices of the compressed row storage format
//		m_rowIndex[0] = 0;
//		int nNonZeros(0);
//		for( int irow = 0; irow < m_numEquations; ++irow ){
//			const unsigned int nColNonZeros = static_cast<unsigned int>( m_matrixTripletFormat[irow].size() );
//			nNonZeros += nColNonZeros;
//			m_rowIndex[ irow + 1 ] = nNonZeros;
//		}
//		m_numNonZeros = nNonZeros;
//				
//		if( m_columns != NULL ){
//			delete[] m_columns;
//			m_columns = NULL;
//		}
//		if( m_values != NULL ){
//			delete[] m_values;
//			m_values = NULL;
//		}
//		m_columns = new int[ m_numNonZeros ];//Columns in which non-zero compnents exist
//		m_values = new std::complex<double>[ m_numNonZeros ];//Values of non-zero compnents
//		for( int i = 0; i < m_numNonZeros; ++i ){
//			m_columns[i] = NULL; // Initialize
//			m_values[i] = std::complex<double>(0.0,0.0); // Initialize
//		}
//
//		//Calculate columns in which non-zero compnents exist
//		int iNonZero(0);
//		for( int irow = 0; irow < m_numEquations; ++irow ){
//			for( std::set<int>::iterator it = m_matrixTripletFormat[irow].begin(); it != m_matrixTripletFormat[irow].end(); ++it ){
//				m_columns[iNonZero] = *it;
//				++iNonZero;
//			}
//		}
//		if( iNonZero != m_numNonZeros ){
//			OutputFiles::m_logFile << "Error : Total number of non-zero componets is wrong." << std::endl;
//			exit(1);
//		}
//
//		// For debug >>>>>
//#ifdef _DEBUG_WRITE
//		for( int i = 0; i < m_numEquations + 1 ; ++i ){
//			std::cout << "i : " << i << " m_rowIndex[i]  : " << m_rowIndex[i] << std::endl;
//		}
//		for( int i = 0; i < m_numNonZeros; ++i ){
//			std::cout << "i : " << i << " m_columns[i]  : " << m_columns[i] << std::endl;
//		}
//#endif
//		// For debug <<<<<
//
//		m_hasConvertedToCRSFormat = true;
//		deleteTripletMatrix();
//	}	
//
//}
//
//// Add non-zero value to matrix
//// Note : This function must be called AFTER the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::addNonZeroValues( const int row, const int col, const std::complex<double> val ){
//
//	if( m_hasConvertedToCRSFormat == false ){
//		//Matrix has not yet been converted to CRS format
//		OutputFiles::m_logFile << "Error :Matrix has not yet been converted to CRS format." << std::endl;
//		exit(1);
//	}
//
//	if( row > m_numEquations - 1 || row < 0 || col > m_numEquations - 1 || col < 0 ){
//		OutputFiles::m_logFile <<  "Error : Row or column is out of the range of the matrix. row : " << row << " col : " << col << std::endl;
//		exit(1);
//	}
//
//	if( pardisoSolver.getSolutionStage() >= PARDISOSolver::FACTORIZED ){
//		zeroClearNonZeroValues();
//		pardisoSolver.setSolutionStage( PARDISOSolver::ANALYZED );
//	}
//
//	//----- Search the column to which value is added -----
//	//for( int inum = m_rowIndex[row]; inum < m_rowIndex[row+1]; ++inum ){
//	//	if( col == m_columns[inum] ){
//	//		m_values[inum] += val;
//	//		return;
//	//	}
//	//}
//	int low = m_rowIndex[row];
//	int high = m_rowIndex[row+1] - 1;
//	while( low <= high ){// binary search
//		const int mid = ( low + high ) / 2;
//		if( m_columns[mid] == col ){// Find column location
//			m_values[mid] += val;
//			return;
//		}else if( m_columns[mid] < col ){
//			low = mid + 1;
//		}else{
//			high = mid - 1;
//		}		
//	}
//
//	//Do not find corresponding location in the matrix strucuture of CSR format.
//	OutputFiles::m_logFile << "Error : Location of non-zero value is improper. row = " << row << " , column = " <<  col << std::endl;
//	exit(1);
//		
//}
//
////Zero clear non-zero values of matrix stored by CSR format
//// Note : This function must be called AFTER the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::zeroClearNonZeroValues(){
//
//	if( m_hasConvertedToCRSFormat == false ){
//		//Matrix has not yet been converted to CRS format
//		OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
//		exit(1);
//	}
//
//	for( int i = 0; i < m_numNonZeros; ++i ){
//		m_values[i] = std::complex<double>(0.0,0.0); // Zero clear
//	}
//
//}
//
//
////Add non-zero value to the right hand side vector
//void ComplexSparseSquareMatrix::addRightHandSideVector( const int row, const std::complex<double> val, const int irhs ){
//
//	if( row < 0 || row >= m_numEquations ){
//		OutputFiles::m_logFile << "Error : No-zero value is tried to set to an improper location of right hand side vector. : row = " << row << std::endl;
//		exit(1);		
//	}
//	
//	if( irhs < 0 || irhs >= m_numRightHandSideVectors ){
//		OutputFiles::m_logFile << "Error : Number of right hand side vector is out of range. : irhs = " << irhs << std::endl;
//		exit(1);		
//	}
//
//	if( pardisoSolver.getSolutionStage() >= PARDISOSolver::SOLVED ){
//		zeroClearRightHandSideVector();
//		pardisoSolver.setSolutionStage( PARDISOSolver::FACTORIZED );
//	}
//
//	m_rightHandSideVector[ row + m_numEquations * irhs ] += val;
//}
//
////Zero clear non-zero values of the right hand side vector
//void ComplexSparseSquareMatrix::zeroClearRightHandSideVector(){
//
//	for( int i = 0; i < m_numEquations*m_numRightHandSideVectors; ++i ){
//		m_rightHandSideVector[i] = std::complex<double>(0.0,0.0); // Zero clear
//	}
//
//}
//

//Initialize matrix and right-hand side vectors
void ComplexSparseSquareMatrix::initializeMatrixAndRhsVectors( const int nEq, const int nRhs ){

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

	ComplexSparseMatrix::initializeMatrixAndRhsVectors( nEq, nEq, nRhs );

}
//
//// Initialize matrix solver
//void ComplexSparseSquareMatrix::initializeMatrixSolver( const std::string oocHeaderName, const int imode ){
//
//	m_pardisoSolver.initialize( oocHeaderName, imode, PARDISOSolver::COMPLEX_AND_UNSYMMETRIC_MATRIX );
//
//}

// Anaysis phase of matrix solver
// [Note] : This function must be called AFTER the matrix is converted into CRS format
void ComplexSparseSquareMatrix::analysisPhaseMatrixSolver(){

	//if( m_hasConvertedToCRSFormat == false ){
	//	//Matrix has not yet been converted to CRS format
	//	OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
	//	exit(1);
	//}

	assert( m_hasConvertedToCRSFormat );

	m_pardisoSolver.analysis( m_numRows, m_rowIndex, m_columns );
}

//Numerical factorization phase of matrix solver
void ComplexSparseSquareMatrix::factorizationPhaseMatrixSolver(){

	//if( m_hasConvertedToCRSFormat == false ){
	//	//Matrix has not yet been converted to CRS format
	//	OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
	//	exit(1);
	//}

	assert( m_hasConvertedToCRSFormat );

	m_pardisoSolver.numericalFactorization( m_rowIndex, m_columns, m_values );
}

//Solve phase of matrix solver with a specified number of right-hand side
void ComplexSparseSquareMatrix::solvePhaseMatrixSolver( std::complex<double>* solution, const long long iRhsStart ,const int nRhs ){

	assert( m_hasConvertedToCRSFormat );

	const long long index = static_cast<long long>(m_numRows) * iRhsStart;
	m_pardisoSolver.solve( m_rowIndex, m_columns, m_values, nRhs, &m_rightHandSideVector[index], solution );
}

//Solve phase of matrix solver
void ComplexSparseSquareMatrix::solvePhaseMatrixSolver( std::complex<double>* solution ){

	//if( m_hasConvertedToCRSFormat == false ){
	//	//Matrix has not yet been converted to CRS format
	//	OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
	//	exit(1);
	//}

	assert( m_hasConvertedToCRSFormat );

	m_pardisoSolver.solve( m_rowIndex, m_columns, m_values, m_numRightHandSideVectors, m_rightHandSideVector, solution );
}

//Release memory of matrix solver
void ComplexSparseSquareMatrix::releaseMemoryMatrixSolver(){

	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}

}

// Get memory required by matrix solver
void ComplexSparseSquareMatrix::writeMemoryRequiredByMatrixSolver(){
	m_pardisoSolver.writeMemoryRequired();
}

//Release memory
void ComplexSparseSquareMatrix::releaseMemory(){

	if( m_pardisoSolver.getSolutionStage() > PARDISOSolver::MEMORY_RELEASED ){
		m_pardisoSolver.releaseMemory();
	}
	ComplexSparseMatrix::releaseMemory();

}

// Get Degree of equation
int ComplexSparseSquareMatrix::getDegreeOfEquation() const{
	return m_numRows;
}

////Get total number of right hand side vector 
//int ComplexSparseSquareMatrix::getNumRightHandSideVectors() const{
//	return m_numRightHandSideVectors;
//}
//
//////Return whether matrix structure has already been set or not
////// [note] : If matrix solver other than PARDISO is implimented, the procedure must be changed
////bool ComplexSparseSquareMatrix::hasMatrixStructureSet() const{
////	if( pardisoSolver.getSolutionStage() < PARDISO::ANALYZED ){
////		return false;
////	}
////	return true;
////}
//
////Return whether matrix has already been converted to CRS format
//bool ComplexSparseSquareMatrix::hasConvertedToCRSFormat() const{
//	return m_hasConvertedToCRSFormat;
//}
//
////Reallocate memory for right hand side vector
//void ComplexSparseSquareMatrix::reallocateMemoryForRightHandSideVectors( const int nrhs ){
//
//	if( nrhs <= 0){
//		OutputFiles::m_logFile << "Error : Number of right-hand sides vectors is less than or equals to zero. rhs = " << nrhs << std::endl;
//		exit(1);
//	}
//
//	if( m_numEquations <= 0){
//		OutputFiles::m_logFile << "Error : Total number of equations is less than or equals to zero. m_numEquations = " << m_numEquations << std::endl;
//		exit(1);
//	}
//
//	if( m_rightHandSideVector != NULL ){
//		delete[] m_rightHandSideVector;
//		m_rightHandSideVector = NULL;
//	}
//
//	m_numRightHandSideVectors = nrhs;
//	m_rightHandSideVector = new std::complex<double>[m_numEquations*m_numRightHandSideVectors];
//
//	for( int i = 0; i < m_numEquations*m_numRightHandSideVectors; ++i ){
//		m_rightHandSideVector[i] = std::complex<double>(0.0,0.0); // Initialize
//	}
//
//}
//
////Set number of right hand side vectors
//void ComplexSparseSquareMatrix::setNumRightHandSideVectors( const int nrhs ){
//	m_numRightHandSideVectors = nrhs;
//}
//
////Multiply matrix by inputed vector and add result to another vector for the case the indexes of inputed vector does not correspond to colunm number
//void ComplexSparseSquareMatrix::multiplyMatrixByVectorAndAddResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const{
//
//	for( int i = 0; i < m_numEquations; ++i ){
//		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
//			vecOut[i] += m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
//		}
//	}
//
//}
//
////Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
//void ComplexSparseSquareMatrix::multiplyMatrixByVectorAndSubtractResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const{
//
//	for( int i = 0; i < m_numEquations; ++i ){
//		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
//			vecOut[i] -= m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
//		}
//	}
//
//}
//
////Substitute right-hand side vector to another vector
//void ComplexSparseSquareMatrix::substituteRhsVector( std::complex<double>* vecOut ) const{
//
//	for( int j = 0; j < m_numRightHandSideVectors; ++j ){
//		for( int i = 0; i < m_numEquations; ++i ){
//			vecOut[i+m_numEquations*j] = m_rightHandSideVector[i+m_numEquations*j];
//		}
//	}
//}
//
//// Debug write the matrix componets
//// [Note] : This function must be called AFTER the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::debugWriteMatrix() const{
//
//	if( m_hasConvertedToCRSFormat == false ){
//		//Matrix has not yet been converted to CRS format
//		OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
//		exit(1);
//	}
//
//	for( int i = 0; i < m_numEquations; ++i ){
//		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
//			std::cout << "row col val " << i << " " << m_columns[j] << " " << m_values[j] << std::endl;
//		}
//	}
//
//}
//
//// Debug write the componets of right hand side vector 
//void ComplexSparseSquareMatrix::debugWriteRightHandSide() const{
//
//	for( int j = 0; j < m_numRightHandSideVectors; ++j ){
//		for( int i = 0; i < m_numEquations; ++i ){
//			std::cout << "row irhs " << i << " " << j << " " << m_rightHandSideVector[i+m_numEquations*j] << std::endl;
//		}
//	}
//
//}

//Copy constructer
ComplexSparseSquareMatrix::ComplexSparseSquareMatrix(const ComplexSparseSquareMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class ComplexSparseSquareMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
ComplexSparseSquareMatrix& ComplexSparseSquareMatrix::operator=(const ComplexSparseSquareMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class ComplexSparseSquareMatrix is not implemented." << std::endl;
	exit(1);
}

//// Delete the matrix of triplet ( Coordinate ) format
//// Note : This function must be called AFTER the matrix is converted into CRS format
//void ComplexSparseSquareMatrix::deleteTripletMatrix(){
//
//	if( m_hasConvertedToCRSFormat == false ){
//		//Matrix has not yet been converted to CRS format
//		OutputFiles::m_logFile << "Warning : Matrix has not yet been converted to CRS format." << std::endl;
//	}
//	else{
//		//m_rowsTriplet.clear();
//		//m_columnsTriplet.clear();
//		//m_valuesTriplet.clear();
//
//		for( int i = 0; i < m_numEquations; ++i ){
//			m_matrixTripletFormat[i].clear();
//		}
//		delete[] m_matrixTripletFormat;
//		m_matrixTripletFormat = NULL;
//	}
//}