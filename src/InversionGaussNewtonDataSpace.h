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
#ifndef DBLDEF_INVERSION_GN_DATA_SPACE
#define DBLDEF_INVERSION_GN_DATA_SPACE

#include "Inversion.h"
#include "RougheningMatrix.h"

// Class of inversion using Gauss-Newton method (data space)
class InversionGaussNewtonDataSpace : public Inversion {

public:
	// Constructer
	explicit InversionGaussNewtonDataSpace();

	// Constructer
	explicit InversionGaussNewtonDataSpace( const int nModel, const int nData );

	// Destructer
	virtual ~InversionGaussNewtonDataSpace();

	// Perform inversion
	virtual void inversionCalculation();

	// Perform inversion by the new method
	void inversionCalculationByNewMethod() const;

	// Perform inversion by the new method using inverse of [R]T[R] matrix
	void inversionCalculationByNewMethodUsingInvRTRMatrix() const;

	// Read sensitivity matrix
	void readSensitivityMatrix( const std::string& fileName, int& numData, int& numModel, double*& sensitivityMatrix ) const;

private:
	// Copy constructer
	InversionGaussNewtonDataSpace( const InversionGaussNewtonDataSpace& rhs ){
		std::cerr << "Error : Copy constructer of the class InversionGaussNewtonDataSpace is not implemented." << std::endl;
		exit(1);
	}

	// Copy assignment operator
	InversionGaussNewtonDataSpace& operator=( const InversionGaussNewtonDataSpace& rhs ){
		std::cerr << "Error : Assignment operator of the class InversionGaussNewtonDataSpace is not implemented." << std::endl;
		exit(1);
	}

	// Calculate constraining matrix for difference filter
	void calcConstrainingMatrixForDifferenceFilter( DoubleSparseMatrix& constrainingMatrix ) const;

};

#endif
