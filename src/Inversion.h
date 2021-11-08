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
#ifndef DBLDEF_INVERSION
#define DBLDEF_INVERSION

#include <iostream>
#include <complex>
#include "Forward3D.h"
#include "RougheningSquareMatrix.h"
#include "DoubleSparseSquareSymmetricMatrix.h"

// Class of inversion
class Inversion{

public:
	enum InversionMethod{
		GAUSS_NEWTON_MODEL_SPECE = 0,
		GAUSS_NEWTON_DATA_SPECE 
	};

	// Constructer
	explicit Inversion();

	// Constructer
	explicit Inversion( const int nModel, const int nData );

	// Destructer
	virtual ~Inversion();

	// Calculate derivatives of EM field
	void calculateDerivativesOfEMField( Forward3D* const ptrForward3D, const double freq, const int iPol );
	
	// Calculate sensitivity matrix
	void calculateSensitivityMatrix( const int freqIDAmongThisPE, const double freq );
	
	// Allocate memory for sensitivity values
	void allocateMemoryForSensitivityScalarValues();
	
	// Release memory of sensitivity values
	void releaseMemoryOfSensitivityScalarValues();
	
	// Output scaler sensitivity values to vtk file
	void outputSensitivityScalarValuesToVtk() const;

	// Output scaler sensitivity values to binary file
	void outputSensitivityScalarValuesToBinary( const int interNum ) const;

	// Perform inversion
	virtual void inversionCalculation() = 0;

	// Delete out-of-core file all
	void deleteOutOfCoreFileAll();

	// Get number of model
	int getNumberOfModel() const;

	// Output number of model to log file
	void outputNumberOfModel() const;

protected:
	// Calculate constraining matrix
	void calcConstrainingMatrix( DoubleSparseMatrix& constrainingMatrix ) const;

	// Copy model transforming jacobian matrix
	void copyModelTransformingJacobian( const int numBlockNotFixed, const int numModel, double* jacobian ) const;

	// Multiply model transforming jacobian matrix
	void multiplyModelTransformingJacobian( const int numData, const int numModel, const double* jacobian, double* matrix ) const;

private:
	// Copy constructer
	Inversion( const Inversion& rhs ){
		std::cerr << "Error : Copy constructer of the class Inversion is not implemented." << std::endl;
		exit(1);
	}

	// Copy assignment operator
	Inversion& operator=( const Inversion& rhs ){
		std::cerr << "Error : Assignment operator of the class Inversion is not implemented." << std::endl;
		exit(1);
	}

	// Number of model
	int m_numModel;

	// Number of data
	int m_numData;

	// Derivatives of EM field
	std::complex<double>* m_derivativesOfEMField[2];

	// Sensitivity values
	double* m_sensitivityScalarValues;

};

#endif
