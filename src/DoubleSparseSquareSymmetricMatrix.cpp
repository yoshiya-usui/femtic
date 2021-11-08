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
#include "DoubleSparseSquareSymmetricMatrix.h"
#include "OutputFiles.h"
#include <assert.h>

//Default Constructer
DoubleSparseSquareSymmetricMatrix::DoubleSparseSquareSymmetricMatrix():
	DoubleSparseSquareMatrix()
{}

// Constructer
DoubleSparseSquareSymmetricMatrix::DoubleSparseSquareSymmetricMatrix( const int nEq, const int nRhs ):
	DoubleSparseSquareMatrix( nEq, nRhs )
{}

// Destructer
DoubleSparseSquareSymmetricMatrix::~DoubleSparseSquareSymmetricMatrix(){
}

// Set matrix structure ( locations of non-zero components ) by triplet format
// Note : This function must be called BEFORE the matrix is converted into CRS format
void DoubleSparseSquareSymmetricMatrix::setStructureByTripletFormat( const int row, const int col ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col );

	DoubleSparseMatrix::setStructureByTripletFormat( row, col );

}

// Set matrix structure ( locations of non-zero components ) and add values by triplet format
void DoubleSparseSquareSymmetricMatrix::setStructureAndAddValueByTripletFormat( const int row, const int col, const double val ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col );

	DoubleSparseMatrix::setStructureAndAddValueByTripletFormat( row, col, val );

}

// Add non-zero value to matrix
// Note : This function must be called AFTER the matrix is converted into CRS format
void DoubleSparseSquareSymmetricMatrix::addNonZeroValues( const int row, const int col, const double val ){

	//if( row > col ){// Only upper triangle components are stored
	//	OutputFiles::m_logFile <<  "Error : Row number must be less than or equal to column number because only upper triangle components can be stored for symmetric matrix. : row = " << row << ", col =  " << col << std::endl;
	//	exit(1);
	//}
	assert( row <= col );

	DoubleSparseMatrix::addNonZeroValues( row, col, val );

}

//Initialize matrix solver
void DoubleSparseSquareSymmetricMatrix::initializeMatrixSolver( const std::string& oocHeaderName, const int imode ){
	m_pardisoSolver.initialize( oocHeaderName, imode, PARDISOSolver::REAL_AND_SYMMETRIC_INDEFINITE );
}

//Copy constructer
DoubleSparseSquareSymmetricMatrix::DoubleSparseSquareSymmetricMatrix(const DoubleSparseSquareSymmetricMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class DoubleSparseSquareSymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
DoubleSparseSquareSymmetricMatrix& DoubleSparseSquareSymmetricMatrix::operator=(const DoubleSparseSquareSymmetricMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class DoubleSparseSquareSymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

