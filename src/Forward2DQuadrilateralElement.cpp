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
#include "Forward2DQuadrilateralElement.h"
#include "OutputFiles.h"
#include "AnalysisControl.h"
#include <iostream>
#include <vector>
#include <algorithm>

// Constructer
Forward2DQuadrilateralElement::Forward2DQuadrilateralElement( const int planeID, const int iPol ):
	Forward2D( planeID, iPol ),
	m_numEquations(0),
	m_numEquationsDegenerated(0),
	m_numNodeTotal2D(0)
{}

// Destructer
Forward2DQuadrilateralElement::~Forward2DQuadrilateralElement(){

}

// Output EM field and responses of boundary planes obtained by 2D analysis
void Forward2DQuadrilateralElement::output2DResult( const double freq, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if( !ptrAnalysisControl->getIsOutput2DResult() ){
		return;
	}

	const int imode = calcMode();// TM or TE mode

	//---------------------------------
	//----- Input output location -----
	//---------------------------------
	std::ifstream inFile( "output_2D_result.dat", std::ios::in );
	if( inFile.fail() ){
		OutputFiles::m_logFile << "File open error : output_2D_result.dat !!" << std::endl;
	}

	//-----------------------------------------------
	//----- Read data from output_2D_result.dat -----
	//-----------------------------------------------
	int ibuf(0);

	inFile >> ibuf;
	if( ibuf <= 0 ){
		OutputFiles::m_logFile << "Error : Number vertical coordinate for output must be greater than zero!! nCalcPointH = " << ibuf << std::endl;
		exit(1);
	}

	const int nCalcPointH = ibuf;
	double* CalcPointH = new double[nCalcPointH];
	for( int i = 0; i < nCalcPointH; ++i ){
		double buf(0.0);
		inFile >> buf;
		CalcPointH[i] = buf * 1000.0;
	}

	inFile >> ibuf;
	if( ibuf <= 0 ){
		OutputFiles::m_logFile << "Error : Number horizontal coordinate for output must be greater than zero!! nCalcPointW = " << ibuf << std::endl;
		exit(1);
	}

	const int nCalcPointW = ibuf;
	double* CalcPointW = new double[nCalcPointW];
	for( int j = 0; j < nCalcPointW; ++j ){
		double buf(0.0);
		inFile >> buf;
		CalcPointW[j] = buf * 1000.0;
	}

	inFile.close();

	//----- Output EM field and response
	fprintf( OutputFiles::m_csvFileFor2DFwd, "Frequency, %lf\n", freq );
	switch (m_polarization){
		case CommonParameters::EX_POLARIZATION:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Ex-Polarization,\n" );
			break;
		case CommonParameters::EY_POLARIZATION:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Ey-Polarization,\n" );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong polarization!! m_polarization = " << m_polarization << std::endl;
			exit(1);
			break;
	}
	switch (m_planeID){
		case MeshData::YZMinus:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Y-Z(Minus),\n" );
			break;
		case MeshData::YZPlus:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Y-Z(Plus),\n" );
			break;
		case MeshData::ZXMinus:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Z-X(Minus),\n" );
			break;
		case MeshData::ZXPlus:
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Z-X(Plus),\n" );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong plane ID !! m_planeID = " << m_planeID << std::endl;
			exit(1);
			break;
	}
	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		fprintf( OutputFiles::m_csvFileFor2DFwd, "     Y[km],      Z[km],        real(Ex),        imag(Ex),        real(Ey),        imag(Ey),        real(Ez),        imag(Ez),        real(Hx),        imag(Hx),        real(Hy),        imag(Hy),        real(Hz),        imag(Hz),         real(Z),         imag(Z),             Rho,           Phase,\n" );
	}else{//XZ Plane
		fprintf( OutputFiles::m_csvFileFor2DFwd, "     X[km],      Z[km],        real(Ex),        imag(Ex),        real(Ey),        imag(Ey),        real(Ez),        imag(Ez),        real(Hx),        imag(Hx),        real(Hy),        imag(Hy),        real(Hz),        imag(Hz),         real(Z),         imag(Z),             Rho,           Phase,\n" );
	}

	const int numH = nCalcPointH;
	const int numW = nCalcPointW;
	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );

	for( int i = 0; i < numH; ++i ){
		for( int j = 0; j < numW; ++j ){
			int elemID2D(0);
			double xi(0.0);
			double eta(0.0);
			if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){// Y-Z plane
				elemID2D = pMeshData->findElementIncludingPointOnYZPlaneAndReturnElemID2D( m_planeID, CalcPointW[j], CalcPointH[i], xi, eta );
			}else{// Z-X plane
				elemID2D = pMeshData->findElementIncludingPointOnZXPlaneAndReturnElemID2D( m_planeID, CalcPointW[j], CalcPointH[i], xi, eta );
			}
			if( elemID2D < 0 ){
				continue;
			}

			std::complex<double> Ex = calcEx( elemID2D, xi, eta, pMeshData ); 
			std::complex<double> Ey = calcEy( elemID2D, xi, eta, pMeshData ); 
			std::complex<double> Ez = calcEz( elemID2D, xi, eta, pMeshData ); 
			std::complex<double> Hx = calcHx( freq, elemID2D, xi, eta, pMeshData ); 
			std::complex<double> Hy = calcHy( freq, elemID2D, xi, eta, pMeshData ); 
			std::complex<double> Hz = calcHz( freq, elemID2D, xi, eta, pMeshData ); 

			std::complex<double> Z(0.0,0.0);
			if(  m_planeID == MeshData::ZXMinus || m_planeID == MeshData::ZXPlus ){//ZX Plane
				if( imode == CommonParameters::TM_MODE ){
					Ey *= std::complex<double>(-1.0, 0.0);
					Z = Ex / Hy;
				}else{
					Hx *= std::complex<double>(-1.0, 0.0);
					Z = Ey / Hx;
				}
			}else{//YZ Plane
				if( imode == CommonParameters::TM_MODE ){
					Z = Ey / Hx;
				}else{
					Z = Ex / Hy;
				}
			}

			const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
			double rho = static_cast<double>( std::pow(std::abs(Z), 2.0) )/ ( CommonParameters::mu * omega );
			double phase = atan2( imag(Z), real(Z) ) * CommonParameters::rad2deg;

			fprintf( OutputFiles::m_csvFileFor2DFwd, "%10.3lf, %10.3lf, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
				CalcPointW[j]/1000.0, CalcPointH[i]/1000.0, real(Ex), imag(Ex), real(Ey), imag(Ey), real(Ez), imag(Ez), real(Hx), imag(Hx), real(Hy), imag(Hy), real(Hz), imag(Hz), real(Z), imag(Z), rho, phase );
		}
	}


}
