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
#include "DoubleSparseMatrix.h"
#include "OutputFiles.h"
#include <assert.h>

//Default Constructer
DoubleSparseMatrix::DoubleSparseMatrix():
	m_numRows(0),
	m_numColumns(0),
	m_numNonZeros(0),
	m_numRightHandSideVectors(1),
	m_hasConvertedToCRSFormat(false),
	m_rowIndex(NULL),
	m_columns(NULL),
	m_values(NULL),
	m_rightHandSideVector(NULL),
	m_matrixTripletFormat(NULL)
{}

// Constructer
DoubleSparseMatrix::DoubleSparseMatrix( const int nrows, const int ncols, const int nrhs ):
	m_numRows(nrows),
	m_numColumns(ncols),
	m_numNonZeros(0),
	m_numRightHandSideVectors(nrhs),
	m_hasConvertedToCRSFormat(false),
	m_rowIndex(NULL),
	m_columns(NULL),
	m_values(NULL),
	m_rightHandSideVector(NULL),
	m_matrixTripletFormat( new std::map< int, double >[nrows] )
{
	assert( nrows > 0 );
	assert( ncols > 0 );
	assert( nrhs > 0 );

	const int num = m_numRows*m_numRightHandSideVectors;
	m_rightHandSideVector = new double[num];
	for( int i = 0; i < num; ++i ){
		m_rightHandSideVector[i] = 0.0; // Initialize
	}	
}

// Destructer
DoubleSparseMatrix::~DoubleSparseMatrix(){

	releaseMemory();

	if( m_rowIndex != NULL ){
		delete[] m_rowIndex;
		m_rowIndex = NULL;
	}

	if( m_columns != NULL ){
		delete[] m_columns;
		m_columns = NULL;
	}

	if( m_values != NULL ){
		delete[] m_values;
		m_values = NULL;
	}

	if( m_rightHandSideVector != NULL ){
		delete[] m_rightHandSideVector;
		m_rightHandSideVector = NULL;
	}

	if( m_matrixTripletFormat != NULL ){
		deleteTripletMatrix();
		m_matrixTripletFormat = NULL;
	}

}

// Set Degree of equation
// Note : This function must be called BEFORE the matrix is converted into CRS format
void DoubleSparseMatrix::setNumRowsAndColumns( const int nrows, const int ncols ){

	assert( !m_hasConvertedToCRSFormat );

	//Total number of rows
	m_numRows = nrows;

	//Total number of columns
	m_numColumns = ncols;

	if( m_matrixTripletFormat != NULL ){
		delete[] m_matrixTripletFormat;
		m_matrixTripletFormat = NULL;
	}

	m_matrixTripletFormat = new std::map< int, double >[m_numRows];

	if( m_rightHandSideVector != NULL ){
		delete[] m_rightHandSideVector;
		m_rightHandSideVector = NULL;
	}

	const int num = m_numRows*m_numRightHandSideVectors;
	m_rightHandSideVector = new double[num];
	for( int i = 0; i < num; ++i ){
		m_rightHandSideVector[i] = 0.0; // Initialize
	}

}

// Set matrix structure ( locations of non-zero components ) by triplet format
// Note : This function must be called BEFORE the matrix is converted into CRS format
void DoubleSparseMatrix::setStructureByTripletFormat( const int row, const int col ){

	assert( !m_hasConvertedToCRSFormat );
	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( col <= m_numColumns - 1 );
	assert( col >= 0  );

	if( m_matrixTripletFormat[row].find(col) == m_matrixTripletFormat[row].end() ){
		// specified column has not been inserted yet
		m_matrixTripletFormat[row].insert( std::map<int, double >::value_type( col, 0.0 ) );
	}

}

// Set matrix structure ( locations of non-zero components ) and add values by triplet format
void DoubleSparseMatrix::setStructureAndAddValueByTripletFormat( const int row, const int col, const double val ){

	assert( !m_hasConvertedToCRSFormat );
	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( col <= m_numColumns - 1 );
	assert( col >= 0  );

	//m_matrixTripletFormat[row].insert( col );
	if( m_matrixTripletFormat[row].find(col) == m_matrixTripletFormat[row].end() ){
		// specified column has not been inserted yet
		m_matrixTripletFormat[row].insert( std::map<int, double >::value_type( col, val ) );
	}else{
		// specified column has already been inserted
		m_matrixTripletFormat[row][col] += val;
	}

}

// Convert matrix from triplet format to CRS format
// Note : This function must be called BEFORE the matrix is converted into CRS format
void DoubleSparseMatrix::convertToCRSFormat(){

	assert( !m_hasConvertedToCRSFormat );

	if( m_rowIndex != NULL ){
		delete[] m_rowIndex;
		m_rowIndex = NULL;
	}
	m_rowIndex = new long long int[ m_numRows + 1 ];//Row indices of the compressed row storage format
	for( int i = 0; i < m_numRows + 1; ++i ){
		m_rowIndex[i] = 0; // Initialize
	}

	// Calculate total number of non-zero components and row indices of the compressed row storage format
	m_rowIndex[0] = 0;
	int nNonZeros(0);
	for( int irow = 0; irow < m_numRows; ++irow ){
		const int nColNonZeros = static_cast<int>( m_matrixTripletFormat[irow].size() );
		nNonZeros += nColNonZeros;
		m_rowIndex[ irow + 1 ] = nNonZeros;
	}
	m_numNonZeros = nNonZeros;
				
	if( m_columns != NULL ){
		delete[] m_columns;
		m_columns = NULL;
	}
	if( m_values != NULL ){
		delete[] m_values;
		m_values = NULL;
	}
	m_columns = new long long int[ m_numNonZeros ];//Columns in which non-zero compnents exist
	m_values = new double[ m_numNonZeros ];//Values of non-zero compnents

	//Calculate columns in which non-zero compnents exist
	int iNonZero(0);
	for( int irow = 0; irow < m_numRows; ++irow ){
		const std::map< int, double >::iterator itEnd = m_matrixTripletFormat[irow].end();
		for( std::map< int, double >::iterator it = m_matrixTripletFormat[irow].begin(); it != itEnd; ++it ){
			m_columns[iNonZero] = it->first;
			m_values[iNonZero] = it->second;
			++iNonZero;
		}
	}

	if( iNonZero != m_numNonZeros ){
		OutputFiles::m_logFile << "Error : Total number of non-zero componets is wrong." << std::endl;
		exit(1);
	}

	// For debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numRows + 1 ; ++i ){
		std::cout << "i : " << i << " m_rowIndex[i]  : " << m_rowIndex[i] << std::endl;
	}
	for( int i = 0; i < m_numNonZeros; ++i ){
		std::cout << "i : " << i << " m_columns[i]  : " << m_columns[i] << std::endl;
	}
	for( int i = 0; i < m_numNonZeros; ++i ){
		std::cout << "i : " << i << " m_values[i]  : " << m_values[i] << std::endl;
	}
	double* temp = new double[m_numRows*m_numRows];
	for( int irow = 0; irow < m_numRows; ++irow ){
		for( int icol = 0; icol < m_numRows; ++icol ){
			temp[ icol + irow * m_numRows ] = 0.0;
		}
	}
	for( int i = 0; i < m_numRows; ++i ){
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			temp[ m_columns[j] + i * m_numRows ] = m_values[j];
		}
	}
	for( int irow = 0; irow < m_numRows; ++irow ){
		for( int icol = 0; icol < m_numRows; ++icol ){
			std::cout << " " << temp[ icol + irow * m_numRows ];
		}
		std::cout << std::endl;
	}
	delete [] temp;
#endif
	// For debug <<<<<

	m_hasConvertedToCRSFormat = true;
	deleteTripletMatrix();

}

// Add non-zero value to matrix
// Note : This function must be called AFTER the matrix is converted into CRS format
void DoubleSparseMatrix::addNonZeroValues( const int row, const int col, const double val ){

	//if( m_hasConvertedToCRSFormat == false ){
	//	//Matrix has not yet been converted to CRS format
	//	OutputFiles::m_logFile << "Error :Matrix has not yet been converted to CRS format." << std::endl;
	//	exit(1);
	//}

	//if( row > m_numRows - 1 || row < 0 || col > m_numColumns - 1 || col < 0 ){
	//	OutputFiles::m_logFile <<  "Error : Row or column is out of the range of the matrix. row : " << row << " col : " << col << std::endl;
	//	exit(1);
	//}
	assert( m_hasConvertedToCRSFormat );
	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( col <= m_numColumns - 1 );
	assert( col >= 0  );

	//----- Search the column to which value is added -----
	//for( int inum = m_rowIndex[row]; inum < m_rowIndex[row+1]; ++inum ){
	//	if( col == m_columns[inum] ){
	//		m_values[inum] += val;
	//		return;
	//	}
	//}
	int low = m_rowIndex[row];
	int high = m_rowIndex[row+1] - 1;
	while( low <= high ){// binary search
		const int mid = ( low + high ) / 2;
		if( m_columns[mid] == col ){// Find column location
			m_values[mid] += val;
			return;
		}else if( m_columns[mid] < col ){
			low = mid + 1;
		}else{
			high = mid - 1;
		}		
	}

	//Do not find corresponding location in the matrix strucuture of CSR format.
	OutputFiles::m_logFile << "Error : Location of non-zero value is improper. row = " << row << " , column = " <<  col << std::endl;
	exit(1);
		
}

//Zero clear non-zero values of matrix stored by CSR format
// Note : This function must be called AFTER the matrix is converted into CRS format
void DoubleSparseMatrix::zeroClearNonZeroValues(){

	//if( m_hasConvertedToCRSFormat == false ){
	//	//Matrix has not yet been converted to CRS format
	//	OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
	//	exit(1);
	//}
	assert( m_hasConvertedToCRSFormat );

	for( int i = 0; i < m_numNonZeros; ++i ){
		m_values[i] =0.0; // Zero clear
	}

}


//Add non-zero value to the right hand side vector
void DoubleSparseMatrix::addRightHandSideVector( const int row, const double val, const int irhs ){

	//if( row < 0 || row >= m_numRows ){
	//	OutputFiles::m_logFile << "Error : No-zero value is tried to set to an improper location of right hand side vector. : row = " << row << std::endl;
	//	exit(1);		
	//}
	//
	//if( irhs < 0 || irhs >= m_numRightHandSideVectors ){
	//	OutputFiles::m_logFile << "Error : Number of right hand side vector is out of range. : irhs = " << irhs << std::endl;
	//	exit(1);		
	//}

	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( irhs <= m_numRightHandSideVectors - 1 );
	assert( irhs >= 0  );

	m_rightHandSideVector[ row + m_numRows * irhs ] += val;
}

//Zero clear non-zero values of the right hand side vector
void DoubleSparseMatrix::zeroClearRightHandSideVector(){

	const int num = m_numRows*m_numRightHandSideVectors;
	for( int i = 0; i < num; ++i ){
		m_rightHandSideVector[i] = 0.0; // Zero clear
	}

}

//Initialize matrix and right-hand side vectors
void DoubleSparseMatrix::initializeMatrixAndRhsVectors( const int nrows, const int ncols, const int nrhs ){

	//if( nrows <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of rows specified is less than or equals to zero. : nrows = " << nrows << std::endl;
	//	exit(1);		
	//}
	//
	//if( ncols <= 0  ){
	//	OutputFiles::m_logFile << "Error : Total number of columns specified is less than or equals to zero. : ncols = " << ncols << std::endl;
	//	exit(1);		
	//}

	//if( nrhs <= 0 ){
	//	OutputFiles::m_logFile << "Error : Total number of right-hand side vectors is specified to be less than or equals to zero. : nrhs = " << nrhs << std::endl;
	//	exit(1);		
	//}

	assert( nrows > 0 );
	assert( ncols > 0 );
	assert( nrhs > 0 );

	if( m_rowIndex != NULL ){
		delete[] m_rowIndex;
		m_rowIndex = NULL;
	}

	if( m_columns != NULL ){
		delete[] m_columns;
		m_columns = NULL;
	}

	if( m_values != NULL ){
		delete[] m_values;
		m_values = NULL;
	}

	if( m_rightHandSideVector != NULL ){
		delete[] m_rightHandSideVector;
		m_rightHandSideVector = NULL;
	}

	if( m_matrixTripletFormat != NULL ){
		deleteTripletMatrix();
		m_matrixTripletFormat = NULL;
	}

	m_numRows = nrows;
	m_numColumns = ncols;
	m_numNonZeros = 0;
	m_numRightHandSideVectors = nrhs;
	m_hasConvertedToCRSFormat = false;
	m_rowIndex = NULL;
	m_columns = NULL;
	m_values = NULL;
	m_matrixTripletFormat = new std::map< int, double >[nrows];

	const int num = m_numRows*m_numRightHandSideVectors;
	m_rightHandSideVector = new double[num];
	for( int i = 0; i < num; ++i ){
		m_rightHandSideVector[i] = 0.0; // Initialize
	}	

}

// Get total number of rows
int DoubleSparseMatrix::getNumRows() const{
	return m_numRows;
}

// Get total number of columns
int DoubleSparseMatrix::getNumColumns() const{
	return m_numColumns;
}

//Get total number of right hand side vector 
int DoubleSparseMatrix::getNumRightHandSideVectors() const{
	return m_numRightHandSideVectors;
}

//Return whether matrix has already been converted to CRS format
bool DoubleSparseMatrix::hasConvertedToCRSFormat() const{
	return m_hasConvertedToCRSFormat;
}

//Release memory
void DoubleSparseMatrix::releaseMemory(){

	if( m_rowIndex != NULL ){
		delete[] m_rowIndex;
		m_rowIndex = NULL;
	}

	if( m_columns != NULL ){
		delete[] m_columns;
		m_columns = NULL;
	}

	if( m_values != NULL ){
		delete[] m_values;
		m_values = NULL;
	}

	if( m_rightHandSideVector != NULL ){
		delete[] m_rightHandSideVector;
		m_rightHandSideVector = NULL;
	}

	if( m_matrixTripletFormat != NULL ){
		deleteTripletMatrix();
		m_matrixTripletFormat = NULL;
	}

	m_numRows = 0;
	m_numColumns = 0;
	m_numNonZeros = 0;
	m_numRightHandSideVectors = 1;
	m_hasConvertedToCRSFormat = false;

}

//Reallocate memory for right hand side vector
void DoubleSparseMatrix::reallocateMemoryForRightHandSideVectors( const int nrhs ){

	//if( nrhs <= 0){
	//	OutputFiles::m_logFile << "Error : Number of right-hand sides vectors is less than or equals to zero. rhs = " << nrhs << std::endl;
	//	exit(1);
	//}

	//if( m_numRows <= 0){
	//	OutputFiles::m_logFile << "Error : Total number of rows is less than or equals to zero. m_numEquations = " << m_numRows << std::endl;
	//	exit(1);
	//}
	assert( nrhs > 0 );
	assert( m_numRows > 0);

	if( m_rightHandSideVector != NULL ){
		delete[] m_rightHandSideVector;
		m_rightHandSideVector = NULL;
	}

	m_numRightHandSideVectors = nrhs;

	const int num = m_numRows*m_numRightHandSideVectors;
	m_rightHandSideVector = new double[num];
	for( int i = 0; i < num; ++i ){
		m_rightHandSideVector[i] = 0.0; // Initialize
	}

}

////Multiply matrix by inputed vector and add result to another vector for the case the indexes of inputed vector does not correspond to colunm number
//void DoubleSparseMatrix::multiplyMatrixByVectorAndAddResult( const double* const vecIn, const int* convertArray, double* const vecOut ) const{
//
//	for( int i = 0; i < m_numRows; ++i ){
//		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
//			vecOut[i] += m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
//		}
//	}
//
//}
//
////Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
//void DoubleSparseMatrix::multiplyMatrixByVectorAndSubtractResult( const double* const vecIn, const int* convertArray, double* const vecOut ) const{
//
//	for( int i = 0; i < m_numRows; ++i ){
//		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
//			vecOut[i] -= m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
//		}
//	}
//
//}

//Copy right-hand side vector to another vector
void DoubleSparseMatrix::copyRhsVector( double* vecOut ) const{

	for( int j = 0; j < m_numRightHandSideVectors; ++j ){
		for( int i = 0; i < m_numRows; ++i ){
			vecOut[i+m_numRows*j] = m_rightHandSideVector[i+m_numRows*j];
		}
	}
}

// Debug write the matrix componets
// [Note] : This function must be called AFTER the matrix is converted into CRS format
void DoubleSparseMatrix::debugWriteMatrix() const{

	if( m_hasConvertedToCRSFormat == false ){
		//Matrix has not yet been converted to CRS format
		OutputFiles::m_logFile << "Error : Matrix has not yet been converted to CRS format." << std::endl;
		exit(1);
	}

	for( int i = 0; i < m_numRows; ++i ){
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			std::cout << "row col val " << i << " " << m_columns[j] << " " << m_values[j] << std::endl;
		}
	}

}

// Debug write the componets of right hand side vector 
void DoubleSparseMatrix::debugWriteRightHandSide() const{

	for( int j = 0; j < m_numRightHandSideVectors; ++j ){
		for( int i = 0; i < m_numRows; ++i ){
			std::cout << "row irhs " << i << " " << j << " " << m_rightHandSideVector[i+m_numRows*j] << std::endl;
		}
	}

}

//Get row indices of the compressed row storage format
int DoubleSparseMatrix::getRowIndexCRS( const int iRow ) const{
	assert( m_hasConvertedToCRSFormat );
	assert( iRow >= 0 && iRow <= m_numRows );
	return m_rowIndex[iRow];
}

//Get column in which non-zero compnents exist
int DoubleSparseMatrix::getColumnsCRS( const int iNonZero ) const{
	assert( m_hasConvertedToCRSFormat );
	assert( iNonZero >= 0 && iNonZero < m_numNonZeros );
	return m_columns[iNonZero];
}

//Get value of non-zero compnents
double DoubleSparseMatrix::getValueCRS( const int iNonZero ) const{
	assert( m_hasConvertedToCRSFormat );
	assert( iNonZero >= 0 && iNonZero < m_numNonZeros );
	return m_values[iNonZero];
}

//Get right hand side vector
double DoubleSparseMatrix::getRightHandSideVector( const int row, const int irhs ) const{
	assert( row <= m_numRows - 1 );
	assert( row >= 0 );
	assert( irhs <= m_numRightHandSideVectors - 1 );
	assert( irhs >= 0  );
	return m_rightHandSideVector[ row + m_numRows * irhs ];
}

// Calculate matrix-vector product of coefficient matrix and inputted vector
void DoubleSparseMatrix::calcMatrixVectorProduct( const double* invVec, double* outVec ) const{

	assert(!m_matrixTripletFormat);

	for( int i = 0; i < m_numRows; ++i ){
		double work(0.0);
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			 work += m_values[j] * invVec[ m_columns[j] ];
		}
		outVec[i] = work;
	}

}

// Calculate matrix-vector product of transposed coefficient matrix and inputted vector
void DoubleSparseMatrix::calcMatrixVectorProductUsingTransposedMatrix( const double* invVec, double* outVec ) const{

	assert(!m_matrixTripletFormat);

	const int nCol = getNumColumns();
	for( int i = 0; i < nCol; ++i ){
		outVec[i] = 0.0;//Zero clear
	}

	for( int i = 0; i < m_numRows; ++i ){
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){
			const int col = m_columns[j];
			outVec[col] += m_values[j] * invVec[i];
		}
	}

}

// Delete the matrix of triplet ( Coordinate ) format
// Note : This function must be called AFTER the matrix is converted into CRS format
void DoubleSparseMatrix::deleteTripletMatrix(){

	if( m_matrixTripletFormat != NULL ){
		for( int i = 0; i < m_numRows; ++i ){
			m_matrixTripletFormat[i].clear();
		}
		delete[] m_matrixTripletFormat;
		m_matrixTripletFormat = NULL;
	}

}

//Copy constructer
DoubleSparseMatrix::DoubleSparseMatrix(const DoubleSparseMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class DoubleSparseMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
DoubleSparseMatrix& DoubleSparseMatrix::operator=(const DoubleSparseMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class DoubleSparseMatrix is not implemented." << std::endl;
	exit(1);
}
