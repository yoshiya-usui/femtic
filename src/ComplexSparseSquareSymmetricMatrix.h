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
#ifndef DBLDEF_COMPLEX_SPARSE_SQUARE_SYMMETRIC_MATRIX
#define DBLDEF_COMPLEX_SPARSE_SQUARE_SYMMETRIC_MATRIX

#include <complex>
#include <set>

#include "ComplexSparseSquareMatrix.h"

class ComplexSparseSquareSymmetricMatrix : public ComplexSparseSquareMatrix {

public:

	//Default Constructer
	explicit ComplexSparseSquareSymmetricMatrix();

	//Constructer
	explicit ComplexSparseSquareSymmetricMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~ComplexSparseSquareSymmetricMatrix();
	
	//Input component of matrix.
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const std::complex<double>& val );
	
	//Add non-zero value to matrix
	virtual void addNonZeroValues( const int row, const int col, const std::complex<double>& val );

	//Initialize matrix solver
	virtual void initializeMatrixSolver( const std::string& oocHeaderName, const int imode );

	//Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void postmultiplyByVectorAndSubtractResult( const std::complex<double>* const vecIn, const int* convertArray, std::complex<double>* const vecOut ) const;
	void postmultiplyByVectorAndSubtractResult( const std::complex<double>* const vecIn, const int numCompsCopied, const int* const compsCopied2Full, const int* const full2CompsCopied, std::complex<double>* const vecOut ) const;

private:
	//Copy constructer
	ComplexSparseSquareSymmetricMatrix(const ComplexSparseSquareSymmetricMatrix &matrix );

	// Assignment operator
	ComplexSparseSquareSymmetricMatrix& operator=(const ComplexSparseSquareSymmetricMatrix& rhs);

};

#endif
