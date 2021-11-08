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
#ifndef DBLDEF_COMPLEX_SPARSE_SQUARE_UNSYMMETRIC_MATRIX
#define DBLDEF_COMPLEX_SPARSE_SQUARE_UNSYMMETRIC_MATRIX

#include <complex>
#include <set>

#include "ComplexSparseSquareMatrix.h"

class ComplexSparseSquareUnsymmetricMatrix : public ComplexSparseSquareMatrix {

public:

	//Default Constructer
	explicit ComplexSparseSquareUnsymmetricMatrix();

	//Constructer
	explicit ComplexSparseSquareUnsymmetricMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~ComplexSparseSquareUnsymmetricMatrix();
	
	//Initialize matrix solver
	virtual void initializeMatrixSolver( const std::string& oocHeaderName, const int imode );

private:
	//Copy constructer
	ComplexSparseSquareUnsymmetricMatrix(const ComplexSparseSquareUnsymmetricMatrix &matrix );

	// Assignment operator
	ComplexSparseSquareUnsymmetricMatrix& operator=(const ComplexSparseSquareUnsymmetricMatrix& rhs);

};

#endif
