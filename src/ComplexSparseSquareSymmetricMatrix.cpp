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
#include "ComplexSparseSquareSymmetricMatrix.h"
#include "OutputFiles.h"
#include <assert.h>

//Default Constructer
ComplexSparseSquareSymmetricMatrix::ComplexSparseSquareSymmetricMatrix():
	ComplexSparseSquareMatrix()
{}

// Constructer
ComplexSparseSquareSymmetricMatrix::ComplexSparseSquareSymmetricMatrix( const int nEq, const int nRhs ):
	ComplexSparseSquareMatrix( nEq, nRhs )
{}

// Destructer
ComplexSparseSquareSymmetricMatrix::~ComplexSparseSquareSymmetricMatrix(){
}

// Set matrix structure ( locations of non-zero components ) by triplet format
// Note : This function must be called BEFORE the matrix is converted into CRS format
void ComplexSparseSquareSymmetricMatrix::setStructureByTripletFormat( const int row, const int col ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col ); 

	ComplexSparseMatrix::setStructureByTripletFormat( row, col );

}

// Set matrix structure ( locations of non-zero components ) and add values by triplet format
void ComplexSparseSquareSymmetricMatrix::setStructureAndAddValueByTripletFormat( const int row, const int col, const std::complex<double>& val ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col ); 

	ComplexSparseMatrix::setStructureAndAddValueByTripletFormat( row, col, val );

}

// Add non-zero value to matrix
// Note : This function must be called AFTER the matrix is converted into CRS format
void ComplexSparseSquareSymmetricMatrix::addNonZeroValues( const int row, const int col, const std::complex<double>& val ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col ); 
		
	ComplexSparseMatrix::addNonZeroValues( row, col, val );

}

// Initialize matrix solver
void ComplexSparseSquareSymmetricMatrix::initializeMatrixSolver( const std::string& oocHeaderName, const int imode ){

	m_pardisoSolver.initialize( oocHeaderName, imode, PARDISOSolver::COMPLEX_AND_SYMMETRIC_MATRIX );

}

// Multiply matrix by inputed vector and subtract calculated vector from another vector
// for the case the indexes of inputed vector does not correspond to colunm number
//void ComplexSparseSquareSymmetricMatrix::postmultiplyByVectorAndSubtractResult( const std::complex<double>* const vecIn, const int* convertArray, std::complex<double>* const vecOut ) const{
void ComplexSparseSquareSymmetricMatrix::postmultiplyByVectorAndSubtractResult( const std::complex<double>* const vecIn, const int numCompsCopied, 
	const int* const compsCopied2Full, const int* const full2CompsCopied, std::complex<double>* const vecOut ) const{

	//for( int i = 0; i < m_numRows; ++i ){// Upper triangle part with diagonal components
	//	for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
	//		vecOut[i] -= m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
	//	}
	//}

	//for( int i = 0; i < m_numRows; ++i ){// Lower triangle part without diagonal components
	//	const std::complex<double> val = vecIn[ convertArray[ i ] ];
	//	for( int j = m_rowIndex[i] + 1; j < m_rowIndex[i+1]; ++j ){	// Without diagonal components
	//		vecOut[ m_columns[j] ] -= m_values[j] * val;
	//	}
	//}

	int i(0);
	int irow(0);
	int j(0);
	//for( i = 0; i < m_numRows; ++i ){// Upper triangle part with diagonal components
	for( i = 0; i < numCompsCopied; ++i ){// Upper triangle part with diagonal components
		irow = compsCopied2Full[i];
		for( j = m_rowIndex[irow]; j < m_rowIndex[irow+1]; ++j ){		
			//vecOut[i] -= m_values[j] * vecIn[ convertArray[ m_columns[j] ] ];
			vecOut[i] -= m_values[j] * vecIn[ m_columns[j] ];
		}
	}

	//for( int i = 0; i < m_numRows; ++i ){// Lower triangle part without diagonal components
	for( int i = 0; i < numCompsCopied; ++i ){// Lower triangle part without diagonal components
		irow = compsCopied2Full[i];
		//const std::complex<double> val = vecIn[ convertArray[ i ] ];
		const std::complex<double> val = vecIn[ irow ];
		for( int j = m_rowIndex[irow] + 1; j < m_rowIndex[irow+1]; ++j ){	// Without diagonal components
			vecOut[ full2CompsCopied[ m_columns[j] ] ] -= m_values[j] * val;
		}
	}

}

// Copy constructer
ComplexSparseSquareSymmetricMatrix::ComplexSparseSquareSymmetricMatrix(const ComplexSparseSquareSymmetricMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class ComplexSparseSquareSymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
ComplexSparseSquareSymmetricMatrix& ComplexSparseSquareSymmetricMatrix::operator=(const ComplexSparseSquareSymmetricMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class ComplexSparseSquareSymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

