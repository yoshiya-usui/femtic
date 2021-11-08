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
#include "Forward2DTriangleElement1stOrderNodeBased.h"
#include "OutputFiles.h"
#include "AnalysisControl.h"
#include <iostream>

// Constructer
Forward2DTriangleElement1stOrderNodeBased::Forward2DTriangleElement1stOrderNodeBased( const int planeID, const int iPol ):
	Forward2DTriangleElementNodeBased( planeID, iPol )
{}

// Destructer
Forward2DTriangleElement1stOrderNodeBased::~Forward2DTriangleElement1stOrderNodeBased(){
}

// Calculate EM fields of boundary planes by 2D forward calculcation with 1st order nodal element
void Forward2DTriangleElement1stOrderNodeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ){

//	const int imode = calcMode();// TM or TE mode
//	if ( imode != CommonParameters::TE_MODE ){
//		OutputFiles::m_logFile << "Error : Only TE mode can be treated in calcEMFieldsOfBoundaryPlanes1stOrderNodeBased ! imode = " << imode << "." << std::endl;
//		exit(1);
//	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "imode " << imode << std::endl;// For debug
//#endif
//
//	if( m_sourceFieldElectric == false ){
//		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in calcEMFieldsOfBoundaryPlanes1stOrderNodeBased !" << std::endl;
//		exit(1);
//	}
//
//	OutputFiles::m_logFile << "# Calculating electric field on a boundary plane with 1st order node-based element." << std::endl;
//
//	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();
//
//	//-------------------------------------------------
//	//--- Total element number on the boudary plane ---
//	//-------------------------------------------------
//	const int nElem = pMeshDataTetraElement->getNumElemOnBoundaryPlanes( m_planeID ); // Total number of elements on the plane
//#ifdef _DEBUG_WRITE
//	std::cout << "nElem " << nElem << std::endl;// For debug
//#endif
//
//	//----------------------------------
//	//--- Set degree of the equation ---
//	//----------------------------------
//	const int nEq = countTotalNodeNumberOnPlane( pMeshDataTetraElement );// Number of equations
//
//	exit(1);

	OutputFiles::m_logFile << "Error : Forward2DTriangleElement1stOrderNodeBased::calcEMFieldsOfBoundaryPlanes has not been implemented yet." << std::endl;
	//exit(1);

}
