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
#ifndef DBLDEF_ROUGHENING_SQUARE_MATRIX
#define DBLDEF_ROUGHENING_SQUARE_MATRIX

#include "DoubleSparseSquareUnsymmetricMatrix.h"
#include "DoubleSparseSquareSymmetricMatrix.h"

class RougheningSquareMatrix : public DoubleSparseSquareUnsymmetricMatrix {

public:

	//Default Constructer
	explicit RougheningSquareMatrix();

	//Constructer
	explicit RougheningSquareMatrix( const int nEq, const int nRhs = 1 );

	//Destructer
	virtual ~RougheningSquareMatrix();
	
	//Input component of matrix.
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const double val );
	
	//Add non-zero value to matrix
	virtual void addNonZeroValues( const int row, const int col, const double val );
	
	//Make [R]T[R] matrix, where [R] is a constraining matrix
	void makeRTRMatrix( DoubleSparseSquareSymmetricMatrix& RTRMatrix, const double smallValueOnDiagonals = 0.0 ) const;

	// Calculate model roughness
	double calcModelRoughness( const double* modelVec ) const;

	// Calculate model roughness for difference filter
	double calcModelRoughnessForDifferenceFilter( const double* modelVec ) const;

	// Calculate vector of model roughness
	void calcVectorOfModelRoughness( const double* modelVec, double* roughnessVec ) const;

	// Postmultiply diagonal matrix
	void postmultiplyDiagonalMatrix( const double* diagMatrix );

	//Calculaste singular values
	void calcSingularValues() const;

	//Calculaste eigen values
	void calcEigenValues() const;

	// Output roughneing matrix
	void outputRougheningMatrix() const;

private:
	//Copy constructer
	RougheningSquareMatrix(const RougheningSquareMatrix &matrix );

	// Assignment operator
	RougheningSquareMatrix& operator=(const RougheningSquareMatrix& rhs);

};

#endif
