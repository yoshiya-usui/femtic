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
#ifndef DBLDEF_FORWARD_2D
#define DBLDEF_FORWARD_2D

#include <complex>
#include "ComplexSparseSquareSymmetricMatrix.h"

// Class of 2D forward calculation
class Forward2D{

public:

	// Constructer
	explicit Forward2D( const int planeID, const int iPol );

	// Destructer
	~Forward2D();

	// Get result of forward analysis directly from solution array;
	std::complex<double> getSolutionDirectly( const int freedum ) const;

	// Get result of forward analysis from element ID and local node(edge) ID
	std::complex<double> getSolutionFromLocalID( const int iElem, const int localID ) const;

protected:

	// Flag specifing the types of 2D analysis
	enum typeOf2DAnalysis{
		NOT_ASSIGNED = -1,
		NODE_BASED_FIRST_ORDER = 0,
		NODE_BASED_SECOND_ORDER,
		EDGE_BASED_ZEROTH_ORDER,
		EDGE_BASED_FIRST_ORDER,
	};

	struct Matrix2x2{
		double mat11;
		double mat12;
		double mat21;
		double mat22;
	};

	// Whether electric field specified at the top of the model as source
	static const bool m_sourceFieldElectric = true;

	// Whether horizontal electric field specified at the left and right side for analysis with edge-based element
	static const bool m_specifyTEResultToSidesOfEdgeElement = false;

	// ID of boundary plane
	int m_planeID;

	// ID of polarization
	int m_polarization;

	//// Type of 2D analysis
	//int m_typeOf2DAnalysis;

	// Array of the matrix of 2D anaysis ( TM & TE mode )
	ComplexSparseSquareSymmetricMatrix m_matrix2DAnalysis;

	// Whether matrix structure has already been set or not
	bool m_hasMatrixStructureSetAndAnalyzed;

	// Solution vector of 2D analysis
	std::complex<double>* m_solution;

	// Array converting local IDs to global ones
	int** m_IDsLocal2Global;

	// Array converting local IDs to global ones after degeneration
	int** m_IDsLocal2GlobalDegenerated;

	// Whether array converting local edge IDs to global ones and  global ones after degeneration has already been set or not
	bool m_hasAlreadySetIDsLocal2Global;

	// Return TM or TE mode from ploralization and plane ID
	int calcMode() const;

private:

	// Defailt constructer
	Forward2D();

	// Copy constructer
	Forward2D(const Forward2D& rhs);

	// Copy assignment operator
	Forward2D& operator=(const Forward2D& rhs);

	// Set ID of boundary plane
	void setPlaneID( const int planeID );
	
	// Initialize matrix solver
	void initializeMatrixSolver();

};

#endif