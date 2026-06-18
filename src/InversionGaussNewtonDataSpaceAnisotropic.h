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
#ifndef DBLDEF_INVERSION_GN_DATA_SPACE_ANISOTROPIC
#define DBLDEF_INVERSION_GN_DATA_SPACE_ANISOTROPIC

#include "Inversion.h"
#include "RougheningMatrix.h"

// Class of inversion using Gauss-Newton method (data space) for anisotropic conductivity
class InversionGaussNewtonDataSpaceAnisotropic : public Inversion {

public:
	// Constructer
	explicit InversionGaussNewtonDataSpaceAnisotropic();

	// Constructer
	explicit InversionGaussNewtonDataSpaceAnisotropic( const int nModel, const int nData );

	// Destructer
	virtual ~InversionGaussNewtonDataSpaceAnisotropic();

	// Perform inversion
	virtual void inversionCalculation();

private:
	// Copy constructer
	InversionGaussNewtonDataSpaceAnisotropic( const InversionGaussNewtonDataSpaceAnisotropic& rhs ){
		std::cerr << "Error : Copy constructer of the class InversionGaussNewtonDataSpace is not implemented." << std::endl;
		exit(1);
	}

	// Copy assignment operator
	InversionGaussNewtonDataSpaceAnisotropic& operator=( const InversionGaussNewtonDataSpaceAnisotropic& rhs ){
		std::cerr << "Error : Assignment operator of the class InversionGaussNewtonDataSpace is not implemented." << std::endl;
		exit(1);
	}

	// Calculate constraining matrices
	void calcConstrainingMatrices( DoubleSparseMatrix& constrainingMatrixFull, DoubleSparseMatrix& constrainingMatrixWithoutCrossProductTerm, 
		DoubleSparseMatrix& derivativeMatrixOfCrossProductTerm, std::vector<double>& crossProductVector ) const;

	// Add the galvanic distortion term to constraining matrix
	void addGalvanicDistortionTermToConstrainingMatrix(const int numberOfUnfixedResistivityParameters, std::vector< std::vector<int> >& nonZeroCols,
		std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const;

	// Add values to constraing matrix
	void addValuesToConstrainingMatrix(const int numModel, const std::vector< std::vector<int> >& nonZeroCols,
		const std::vector< std::vector<double> >& matValues, const std::vector<double>& rhsValues, DoubleSparseMatrix& constrainingMatrix) const;

	// Calculate predicted vector of weighted observed data for each PE
	void calculatePredictedVectorOfWeightedObservedDataLocal( const int numDataThisPE, const int nFreqThisPE, const int numModel, const bool useBLAS,
		const double* const resultVector, double* predictedVectorLocal) const;

};

#endif
