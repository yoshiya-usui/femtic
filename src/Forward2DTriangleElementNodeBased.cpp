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
#include "Forward2DTriangleElementNodeBased.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "ResistivityBlock.h"

// Constructer
Forward2DTriangleElementNodeBased::Forward2DTriangleElementNodeBased( const int planeID, const int iPol ):
	Forward2DTriangleElement( planeID, iPol )
{}

// Destructer
Forward2DTriangleElementNodeBased::~Forward2DTriangleElementNodeBased(){
}

// Calculate Ex
std::complex<double> Forward2DTriangleElementNodeBased::calcEx( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcEx is not implemented." << std::endl;
	exit(1);
}

// Calculate Ey
std::complex<double> Forward2DTriangleElementNodeBased::calcEy( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcEy is not implemented." << std::endl;
	exit(1);
}

// Calculate Ez
std::complex<double> Forward2DTriangleElementNodeBased::calcEz( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcEz is not implemented." << std::endl;
	exit(1);
}

// Calculate Hx
std::complex<double> Forward2DTriangleElementNodeBased::calcHx( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcHx is not implemented." << std::endl;
	exit(1);
}

// Calculate Hy
std::complex<double> Forward2DTriangleElementNodeBased::calcHy( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcHy is not implemented." << std::endl;
	exit(1);
}

// Calculate Hz
std::complex<double> Forward2DTriangleElementNodeBased::calcHz( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{
	OutputFiles::m_logFile << "Error : Forward2DTriangleElementNodeBased::calcHz is not implemented." << std::endl;
	exit(1);
}