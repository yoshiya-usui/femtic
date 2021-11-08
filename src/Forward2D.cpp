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
#include "Forward2D.h"
#include "MeshData.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include <assert.h>

//// Defailt constructer
//Forward2D::Forward2D():
//	m_planeID(0),
//	m_typeOf2DAnalysis(Forward2D::NOT_ASSIGNED)
//{
//	m_hasMatrixStructureSetAndAnalyzed = false;
//	m_solution = NULL;
//	m_IDsLocal2Global = NULL;
//	m_IDsLocal2GlobalDegenerated = NULL;
//	m_hasAlreadySetIDsLocal2Global = false;
//}

// Constructer
Forward2D::Forward2D( const int planeID, const int iPol ):
	m_planeID(planeID),
	m_polarization(iPol)
{
	m_hasMatrixStructureSetAndAnalyzed = false;
	m_solution = NULL;
	m_IDsLocal2Global = NULL;
	m_IDsLocal2GlobalDegenerated = NULL;
	m_hasAlreadySetIDsLocal2Global = false;

	setPlaneID( planeID );

	initializeMatrixSolver();
}

//Destructer
Forward2D::~Forward2D(){

	m_matrix2DAnalysis.releaseMemoryMatrixSolver();

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}

	if( m_IDsLocal2Global != NULL ){
		const int nElem = sizeof( m_IDsLocal2Global ) / sizeof( m_IDsLocal2Global[0] );
		for( int iElem = 0; iElem < nElem; ++iElem ){
			delete [] m_IDsLocal2Global[iElem];
		}
		delete [] m_IDsLocal2Global;
		m_IDsLocal2Global = NULL;
	}

	if( m_IDsLocal2GlobalDegenerated != NULL ){
		const int nElem = sizeof( m_IDsLocal2GlobalDegenerated ) / sizeof( m_IDsLocal2GlobalDegenerated[0] );
		for( int iElem = 0; iElem < nElem; ++iElem ){
			delete [] m_IDsLocal2GlobalDegenerated[iElem];
		}
		delete [] m_IDsLocal2GlobalDegenerated;
		m_IDsLocal2GlobalDegenerated = NULL;
	}

}

// Set ID of boundary plane
void Forward2D::setPlaneID( const int planeID ){
	
	m_planeID = planeID;

	if( m_planeID != MeshData::YZMinus &&
		m_planeID != MeshData::YZPlus  &&
		m_planeID != MeshData::ZXMinus &&
		m_planeID != MeshData::ZXPlus ){
		OutputFiles::m_logFile << "Error : Unsupported plane " <<  m_planeID << "." << std::endl;
		exit(1);
	}

}

// Initialize matrix solver
void Forward2D::initializeMatrixSolver(){

	const int myPE = ( AnalysisControl::getInstance() )->getMyPE();

	m_hasMatrixStructureSetAndAnalyzed = false;
	m_solution = NULL;

	std::ostringstream oocHeaderName;
	oocHeaderName << "ooc_temp_2D_Plane" << m_planeID << "_Pol" << m_polarization << "_PE" << myPE;
	m_matrix2DAnalysis.initializeMatrixSolver( oocHeaderName.str(), PARDISOSolver::INCORE_MODE );

}

// Get result of forward analysis
std::complex<double> Forward2D::getSolutionDirectly( const int freedum ) const {

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : Result of forward analysis is NULL !!" << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );

	return m_solution[freedum];

}

// Get result of forward analysis from element ID and local node(edge) ID
std::complex<double> Forward2D::getSolutionFromLocalID( const int iElem, const int localID ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : Result of forward analysis is NULL !!" << std::endl;
	//	exit(1);
	//}else if( m_IDsLocal2Global[iElem] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	return m_solution[ m_IDsLocal2Global[iElem][localID] ];

}

// Return TM or TE mode from ploralization and plane ID
int Forward2D::calcMode() const{

	if( m_polarization == CommonParameters::EX_POLARIZATION ){ // EX Polarization
		if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){ //YZ Plane
			return CommonParameters::TE_MODE;
		}else{//ZX Plane
			return CommonParameters::TM_MODE;
		}
	}else if( m_polarization == CommonParameters::EY_POLARIZATION ){ // EY Polarization
		if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){ //YZ Plane
			return CommonParameters::TM_MODE;
		}else{//ZX Plane
			return CommonParameters::TE_MODE;
		}
	}else{
		OutputFiles::m_logFile << "Error : Wrong ploralization ID " <<  m_polarization << "." << std::endl;
		exit(1);
	}

}

// Copy constructer
Forward2D::Forward2D(const Forward2D& rhs){
	std::cerr << "Error : Copy constructer of the class Forward2D is not implemented." << std::endl;
	exit(1);
}

// Copy assignment operator
Forward2D& Forward2D::operator=(const Forward2D& rhs){
	std::cerr << "Error : Copy constructer of the class Forward2D is not implemented." << std::endl;
	exit(1);
}
