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
#include "Forward2DSquareElementNodeBased.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "ResistivityBlock.h"

// Local coordinate of vertical direction on top side
const double Forward2DSquareElementNodeBased::m_hLocalTopSide = -1.0;

// Local coordinate of vertical direction on bottom side
const double Forward2DSquareElementNodeBased::m_hLocalBottomSide = 1.0;

//// Defailt constructer
//Forward2DSquareElementNodeBased::Forward2DSquareElementNodeBased()
//{}

// Constructer
Forward2DSquareElementNodeBased::Forward2DSquareElementNodeBased( const int planeID, const int iPol ):
	Forward2DSquareElement( planeID, iPol )
{}

//Destructer
Forward2DSquareElementNodeBased::~Forward2DSquareElementNodeBased()
{}

// Calculate parameter Eta of Rodi(1976)
// [Input] : 1) imode : TM/TE
//           2) freq : Frequency
//           3) iElem : Element ID of the boundary plane, that is a numerical sequence beginning with zero
// [Output] : parameter Eta of Rodi(1976)
std::complex<double> Forward2DSquareElementNodeBased::calcEta( const int imode, const double freq, const int iElem ) const{

	const int elemID = ( ( AnalysisControl::getInstance() )->getPointerOfMeshData() )->getElemBoundaryPlanes( m_planeID, iElem );

	std::complex<double> eta(0,0);
	if( imode == CommonParameters::TM_MODE ){// TM mode
		//--- Calculate sigma
		//ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();
		const double sigma = ( ResistivityBlock::getInstance() )->getConductivityValuesFromElemID(elemID);

		//--- Calculate eta
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		eta = std::complex<double>( sigma, -omega * CommonParameters::epsilon );
	}else if( imode == CommonParameters::TE_MODE ){// TE mode
		//--- Calculate eta
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		eta = std::complex<double>( 0.0, -omega * CommonParameters::mu );
	}else{
		OutputFiles::m_logFile << "Error : Wrong mode" << imode << "." << std::endl;
		exit(1);
	}

	return eta;
}

// Calculate parameter Gamma of Rodi(1976)
// [Input] : 1) imode : TM/TE
//           2) freq : Frequency
//           3) iElem : Element ID of the boundary plane, that is a numerical sequence beginning with zero
// [Output] : parameter Gamma of Rodi(1976)
std::complex<double> Forward2DSquareElementNodeBased::calcGamma( const int imode, const double freq, const int iElem ) const{

	std::complex<double> gamma(0,0);
	if( imode == CommonParameters::TM_MODE ){// TM mode
		gamma = calcEta( CommonParameters::TE_MODE, freq, iElem );
	}else if( imode == CommonParameters::TE_MODE ){// TE mode
		gamma = calcEta( CommonParameters::TM_MODE, freq, iElem );
	}else{
		OutputFiles::m_logFile << "Error : Wrong mode" << imode << "." << std::endl;
		exit(1);
	}

	return gamma;

}

// Calculate Ex
std::complex<double> Forward2DSquareElementNodeBased::calcEx( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return calcValueElectricFieldPerpendicular( iElem, wCoord, hCoord );

	}else{//XZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}

}

// Calculate Ey
std::complex<double> Forward2DSquareElementNodeBased::calcEy( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}else{//XZ Plane
		
		return calcValueElectricFieldPerpendicular( iElem, wCoord, hCoord );

	}

}

// Calculate Ez
std::complex<double> Forward2DSquareElementNodeBased::calcEz( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	return std::complex<double>( 0.0, 0.0 ); 

}

// Calculate Hx
std::complex<double> Forward2DSquareElementNodeBased::calcHx( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}else{//XZ Plane
		
		return calcValueMagneticFieldHorizontal( freq, iElem, wCoord, hCoord, pMeshDataBrickElement );

	}

}

// Calculate Hy
std::complex<double> Forward2DSquareElementNodeBased::calcHy( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		
		return calcValueMagneticFieldHorizontal( freq, iElem, wCoord, hCoord, pMeshDataBrickElement );

	}else{//XZ Plane
		
		return std::complex<double>( 0.0, 0.0 ); 

	}

}

// Calculate Hz
std::complex<double> Forward2DSquareElementNodeBased::calcHz( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	return calcValueMagneticFieldVertical( freq, iElem, wCoord, hCoord, pMeshDataBrickElement );

}

