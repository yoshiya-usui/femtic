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
#include "Forward2DSquareElementEdgeBased.h"

//// Defailt constructer
//Forward2DSquareElementEdgeBased::Forward2DSquareElementEdgeBased()
//{}

// Constructer
Forward2DSquareElementEdgeBased::Forward2DSquareElementEdgeBased( const int planeID, const int iPol ):
	Forward2DSquareElement( planeID, iPol )
{}

//Destructer
Forward2DSquareElementEdgeBased::~Forward2DSquareElementEdgeBased()
{}

// Calculate Ex
std::complex<double> Forward2DSquareElementEdgeBased::calcEx( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}else{//XZ Plane
		
		return calcValueElectricFieldHorizontal( iElem, wCoord, hCoord );

	}

}

// Calculate Ey
std::complex<double> Forward2DSquareElementEdgeBased::calcEy( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return calcValueElectricFieldHorizontal( iElem, wCoord, hCoord );

	}else{//XZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}

}

// Calculate Ez
std::complex<double> Forward2DSquareElementEdgeBased::calcEz( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	return calcValueElectricFieldVertical( iElem, wCoord, hCoord );

}

// Calculate Hx
std::complex<double> Forward2DSquareElementEdgeBased::calcHx( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		return calcValueMagneticFieldPerpendicular( freq, iElem, wCoord, hCoord, pMeshDataBrickElement );
	}else{//XZ Plane
		return std::complex<double>( 0.0, 0.0 ); 
	}

}

// Calculate Hy
std::complex<double> Forward2DSquareElementEdgeBased::calcHy( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		return std::complex<double>( 0.0, 0.0 ); 
	}else{//XZ Plane
		return calcValueMagneticFieldPerpendicular( freq, iElem, wCoord, hCoord, pMeshDataBrickElement );
	}

}

// Calculate Hz
std::complex<double> Forward2DSquareElementEdgeBased::calcHz( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	return std::complex<double>( 0.0, 0.0 ); 

}