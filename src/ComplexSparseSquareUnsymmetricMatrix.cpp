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
#include "ComplexSparseSquareUnsymmetricMatrix.h"
#include "OutputFiles.h"

//Default Constructer
ComplexSparseSquareUnsymmetricMatrix::ComplexSparseSquareUnsymmetricMatrix():
	ComplexSparseSquareMatrix()
{}

// Constructer
ComplexSparseSquareUnsymmetricMatrix::ComplexSparseSquareUnsymmetricMatrix( const int nEq, const int nRhs ):
	ComplexSparseSquareMatrix( nEq, nRhs )
{}

// Destructer
ComplexSparseSquareUnsymmetricMatrix::~ComplexSparseSquareUnsymmetricMatrix(){
}

// Initialize matrix solver
void ComplexSparseSquareUnsymmetricMatrix::initializeMatrixSolver( const std::string& oocHeaderName, const int imode ){

	m_pardisoSolver.initialize( oocHeaderName, imode, PARDISOSolver::COMPLEX_AND_UNSYMMETRIC_MATRIX );

}

//Copy constructer
ComplexSparseSquareUnsymmetricMatrix::ComplexSparseSquareUnsymmetricMatrix(const ComplexSparseSquareUnsymmetricMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class ComplexSparseSquareUnsymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
ComplexSparseSquareUnsymmetricMatrix& ComplexSparseSquareUnsymmetricMatrix::operator=(const ComplexSparseSquareUnsymmetricMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class ComplexSparseSquareUnsymmetricMatrix is not implemented." << std::endl;
	exit(1);
}

