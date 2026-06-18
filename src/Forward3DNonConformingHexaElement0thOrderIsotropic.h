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
#ifndef DBLDEF_FORWARD_3D_NONCONFORMING_HEXA_ELEMENT_0TH_ORDER_ISOTROPIC
#define DBLDEF_FORWARD_3D_NONCONFORMING_HEXA_ELEMENT_0TH_ORDER_ISOTROPIC

#include "Forward3DNonConformingHexaElement0thOrder.h"

// Class of 3D isotropic forward calculation using the 0th order nonconforming hexahedral element 
class Forward3DNonConformingHexaElement0thOrderIsotropic : public Forward3DNonConformingHexaElement0thOrder {

public:

	// Constructer
	Forward3DNonConformingHexaElement0thOrderIsotropic();

	// Destructer
	virtual ~Forward3DNonConformingHexaElement0thOrderIsotropic();

	// Calculate electric current density vector
	virtual CommonParameters::ComplexValuedVector calculateElectricCurrentDensityVector(const int iElem) const;

private:

	// Copy constructer
	Forward3DNonConformingHexaElement0thOrderIsotropic(const Forward3DNonConformingHexaElement0thOrderIsotropic& rhs);

	// Copy assignment operator
	Forward3DNonConformingHexaElement0thOrderIsotropic& operator=(const Forward3DNonConformingHexaElement0thOrderIsotropic& rhs);

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	virtual void setNonZeroValues(ComplexSparseSquareSymmetricMatrix& matrix);

	// Calculate vector x of the reciprocity algorithm of Rodi (1976) for isotropic conductivity
	virtual void calVectorXOfReciprocityAlgorithmForIsotropicConductivity(const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows);

	// Calculate vector x of the reciprocity algorithm of Rodi (1976) for anisotropic conductivity
	virtual void calVectorXOfReciprocityAlgorithmForAnisotropicConductivity(const std::complex<double>* const vecIn, const int blkID, const int paramID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows);

};

#endif
