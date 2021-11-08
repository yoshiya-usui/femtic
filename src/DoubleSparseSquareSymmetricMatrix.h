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
#ifndef DBLDEF_DOUBLE_SPARSE_SQUARE_SYMMETRIC_MATRIX
#define DBLDEF_DOUBLE_SPARSE_SQUARE_SYMMETRIC_MATRIX

#include <set>

#include "DoubleSparseSquareMatrix.h"

class DoubleSparseSquareSymmetricMatrix : public DoubleSparseSquareMatrix {

public:

	//Default Constructer
	explicit DoubleSparseSquareSymmetricMatrix();

	//Constructer
	explicit DoubleSparseSquareSymmetricMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~DoubleSparseSquareSymmetricMatrix();
	
	//Input component of matrix.
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const double val );
	
	//Add non-zero value to matrix
	virtual void addNonZeroValues( const int row, const int col, const double val );

	//Initialize matrix solver
	virtual void initializeMatrixSolver( const std::string& oocHeaderName, const int imode );

private:
	//Copy constructer
	DoubleSparseSquareSymmetricMatrix(const DoubleSparseSquareSymmetricMatrix &matrix );

	// Assignment operator
	DoubleSparseSquareSymmetricMatrix& operator=(const DoubleSparseSquareSymmetricMatrix& rhs);

};

#endif
