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
#include "Forward3D.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "MeshDataBrickElement.h"
#include "ResistivityBlock.h"
#include <assert.h>
#include <algorithm>

#ifdef _USE_OMP
#include <omp.h>
#endif

Forward3D::Forward3D():
	m_numOfEquation(0),
	m_numOfEquationDegenerated(0),
	m_IDsLocal2Global(NULL),
	m_hasSetIDsLocal2Global(false),
	m_hasMatrixStructureSetAndAnalyzed(false),
	m_matrix3DAnalysis(),
	m_solution(NULL),
	m_polarizationCurrent(CommonParameters::EX_POLARIZATION),
	m_frequencyCurrent(0.0),
	m_orderOfFiniteElement(0)
{

	for( int iPol = 0; iPol < 2; ++iPol ){
		m_IDsGlobal2AfterDegenerated[iPol] = NULL;
		m_hasIDsGlobal2AfterDegenerated[iPol] = false;
	}

}

//Destructer
Forward3D::~Forward3D(){

	if( m_IDsLocal2Global != NULL ){
		const int num = sizeof( m_IDsLocal2Global ) / sizeof( m_IDsLocal2Global[0] );
		for( int i = 0; i < num; ++i ){
			delete [] m_IDsLocal2Global[i];
			m_IDsLocal2Global[i] = NULL;
		}
		delete [] m_IDsLocal2Global;
		m_IDsLocal2Global = NULL;
	}

	for( int iPol = 0; iPol < 2; ++iPol ){
		if( m_IDsGlobal2AfterDegenerated[iPol] != NULL ){
			delete [] m_IDsGlobal2AfterDegenerated[iPol];
			m_IDsGlobal2AfterDegenerated[iPol] = NULL;
		}
	}
	
	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}

	m_matrix3DAnalysis.releaseMemoryMatrixSolver();

}

//Run 3D forward calculation
void Forward3D::forwardCalculation( const double freq, const int iPol ){
	OutputFiles::m_logFile << "Error : forwardCalculation is not implemented in the class Forward3D." << std::endl;
	exit(1);
}

// Set polarization at present for which forward analysis is executed 
void Forward3D::setPolarizationCurrent( const int iPol ){
	m_polarizationCurrent = iPol;
}

// Set frequency at present for which forward analysis is executed 
void Forward3D::setFrequencyCurrent( const double freq ){
	m_frequencyCurrent = freq;
}

// Set order of finite element
void Forward3D::setOrderOfFiniteElement( const int order ){
	m_orderOfFiniteElement = order;
}

//// Get columns number of right-hand sides matrix consisting of interpolator vectors
//int Forward3D::getColumnsNumberRhsMatrixConsistingOfInterpolatorVectors() const{
//	return m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors;
//}
//
//// Increment columns number of right-hand sides matrix consisting of interpolator vectors
//void Forward3D::incrementColumnsNumberRhsMatrixConsistingOfInterpolatorVectors(){
//	++m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors;
//}

// Add values to right-hand sides matrix consisting of interpolator vectors
void Forward3D::addValuesToRhsVectors( const int irow, const int irhs, const std::complex<double>& val ){

	//if( irow >= m_numOfEquationDegenerated || irow < 0 ){
	//	OutputFiles::m_logFile << "Error : Row number of right-hand sides matrix is out of range. irow = " << irow << std::endl;
	//	exit(1);
	//}
	//if( irhs >= m_matrix3DAnalysis.getNumRightHandSideVectors() || irhs < 0 ){
	//	OutputFiles::m_logFile << "Error : Column number of right-hand sides matrix is out of range. irhs = " << irhs << std::endl;
	//	exit(1);
	//}

	//m_rhsMatrixConsistingOfInterpolatorVectors[ irow + m_numOfEquationDegenerated * icolumn ] += val;
	m_matrix3DAnalysis.addRightHandSideVector( irow, val, irhs );

}

//// Allocate memory and initialize right-hand sides matrix consisting of interpolator vectors
//void Forward3D::initializeRhsMatrixConsistingOfInterpolatorVectors( const int nrhs ){
// Initialize right-hand side vectors
void Forward3D::initializeRhsVectors( const int nrhs ){

	//if( nrhs <= 0){
	//	OutputFiles::m_logFile << "Error : Number of interpolator vectors is less than or equals to zero !! nrhs = " << nrhs << std::endl;
	//	exit(1);
	//}

	//m_numColumnsRhsMatrixConsistingOfInterpolatorVectors = nrhs;
	//m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors = 0;

	//if( m_rhsMatrixConsistingOfInterpolatorVectors != NULL ){
	//	delete [] m_rhsMatrixConsistingOfInterpolatorVectors;
	//	m_rhsMatrixConsistingOfInterpolatorVectors = NULL;
	//}

	//m_rhsMatrixConsistingOfInterpolatorVectors = new std::complex<double>[ m_numColumnsRhsMatrixConsistingOfInterpolatorVectors * m_numOfEquationDegenerated ];

	//for( int i = 0; i < m_numColumnsRhsMatrixConsistingOfInterpolatorVectors * m_numOfEquationDegenerated; ++i ){
	//	m_rhsMatrixConsistingOfInterpolatorVectors[i] = std::complex<double>( 0.0, 0.0 );
	//}

	m_matrix3DAnalysis.reallocateMemoryForRightHandSideVectors( nrhs );

}

// Perform solve phase for right-hand sides consisting of interpolator vectors
void Forward3D::solvePhaseForRhsConsistingInterpolatorVectors( const int numInterpolatorVectors, std::complex<double>* solutionForInterpolatorVectors ){

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	m_matrix3DAnalysis.debugWriteMatrix();
	m_matrix3DAnalysis.debugWriteNonZeroRightHandSide();
#endif
	//----- debug <<<<<

	if( numInterpolatorVectors != m_matrix3DAnalysis.getNumRightHandSideVectors() ){
		OutputFiles::m_logFile << "Error : Number of interpolator vectors inputed is not equal to number of right-hand side vectors!!. numInterpolatorVectors = " << numInterpolatorVectors << ", number of right-hand side vectors = "  << m_matrix3DAnalysis.getNumRightHandSideVectors() << std::endl;
		exit(1);
	}

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	int numDivRhs = ptrAnalysisControl->getDivisionNumberOfMultipleRHSInForward();

	assert( numDivRhs > 0 ); 

	if( numDivRhs > numInterpolatorVectors ){
		OutputFiles::m_logFile << "Warning : Division number of right-hand sides ( " << numDivRhs << " ) is greater than the total number of right-hand sides ( " << numInterpolatorVectors <<  " )." << std::endl;
		OutputFiles::m_logFile << "          Thus, the division number is set the total number of right-hand sides." << std::endl;
		numDivRhs = numInterpolatorVectors;
	}

	const int numOfEquationFinallySolved = getNumOfEquationFinallySolved();
	const int numRHSDividedWithoutOdds = numInterpolatorVectors / numDivRhs;
	const int numAdds = numInterpolatorVectors % numDivRhs;
	long long iRhsStart = 0;
	for( int iDiv = 0; iDiv < numDivRhs; ++iDiv ){
		const int numRHSDividedWithout = iDiv < numAdds ? numRHSDividedWithoutOdds + 1 : numRHSDividedWithoutOdds;
		OutputFiles::m_logFile << "# Solve phase is performed simultaneously for " << numRHSDividedWithout  << " right-hand sides" << ptrAnalysisControl->outputElapsedTime() << std::endl;
		const int long long index = static_cast<long long>(numOfEquationFinallySolved) * iRhsStart;
		m_matrix3DAnalysis.solvePhaseMatrixSolver( &solutionForInterpolatorVectors[index], iRhsStart, numRHSDividedWithout );
		iRhsStart += static_cast<long long>(numRHSDividedWithout);
	}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	std::cout << "solution for right-hand sides consisting of interpolator vectors" << std::endl;
	for( int j = 0; j < numInterpolatorVectors; ++j ){
		for( int i = 0; i < numOfEquationFinallySolved; ++i ){
			std::cout << "row irhs " << i << " " << j << " " << solutionForInterpolatorVectors[i+numOfEquationFinallySolved*j] << std::endl;
		}
	}
#endif
	//----- debug <<<<<

	// Reallocate memory for right-hand side vector
	initializeRhsVectors(1);

}

// Calculate  derivative of EM field
void Forward3D::calculateDerivativesOfEMField( const int numInterpolatorVectors,
	const std::complex<double>* const solutionForInterpolatorVectors, std::complex<double>* const derivatives ){

	const int numOfEquationFinallySolved = getNumOfEquationFinallySolved();

	const int iPol = getPolarizationCurrent();
	std::complex<double>* solutionAfterDegenerated = new std::complex<double>[ numOfEquationFinallySolved ];
	copySolutionVectorDegenerated(iPol, solutionAfterDegenerated);

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0 ; i < numOfEquationFinallySolved; ++i ){
		std::cout << "solutionAfterDegenerated[" << i << "]] = " << solutionAfterDegenerated[i] << std::endl;
	}
	for( int i = 0 ; i < numOfEquationFinallySolved; ++i ){
		std::cout << "solutionForInterpolatorVectors[" << i << "]] = " << solutionForInterpolatorVectors[i] << std::endl;
	}
#endif
	//----- debug <<<<<

	const ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();
	const int nBlkNotFixed = ptrResistivityBlock->getNumResistivityBlockNotFixed();
	const int nBlkTotal = ptrResistivityBlock->getNumResistivityBlockTotal();

	int numThreads = 1;
#ifdef _USE_OMP
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
#endif

	//--------------------------------------
	// Allocate variables for each thread
	//--------------------------------------
	std::complex<double>** nonZeroValues = new std::complex<double>*[numThreads];
	for( int i = 0; i < numThreads; ++i ){
		nonZeroValues[i] = new std::complex<double>[numOfEquationFinallySolved];
	}
	std::vector<int>* nonZeroComps = new std::vector<int>[numThreads];
	//--------------------------------------
	//--------------------
	// Private variables
	//--------------------
	int iThread = 0;
	int iblk = 0;
	std::vector<int>::iterator itrEnd;
	std::vector<int>::iterator itr;
	int imdl = 0;
	int ivec(0);
	std::complex<double> work(0.0,0.0);
	int i(0);
	int j(0);
	long long index(0);
#ifdef _USE_OMP
	//--------------------- Start of omp parallel >>>>>
	#pragma omp parallel private( iThread, iblk, itrEnd, itr, imdl, ivec, work, i, j, index )
#endif
	{ 

#ifdef _USE_OMP
	//--------------------- Start of omp for >>>>>
	iThread = omp_get_thread_num();

	#pragma omp single
	{
		if( numThreads != omp_get_num_threads() ){
			OutputFiles::m_logFile << "Error : Number of threads is different from the one obtained previously." << std::endl;
			exit(1);
		}
	}

	#pragma omp for
#endif
	for( iblk = 0; iblk < nBlkTotal; ++iblk ){

		if( ptrResistivityBlock->isFixedResistivityValue( iblk ) ){
			continue;
		}

		//----------------------
		//--- Initialization ---
		//----------------------
		nonZeroComps[iThread].clear();
		for( i = 0; i < numOfEquationFinallySolved; ++i ){
			nonZeroValues[iThread][i] = std::complex<double>(0.0,0.0);
		}
		//----------------------

		calVectorXOfReciprocityAlgorithm( solutionAfterDegenerated, iblk, nonZeroValues[iThread], nonZeroComps[iThread] );

//#ifdef _DEBUG_WRITE
//		for( i = 0; i < numOfEquationFinallySolved; ++i ){
//			std::cout << "iThread iblk i nonZeroValues " << iThread << " " << iblk << " " << i << " " << nonZeroValues[iThread][i] << std::endl;
//		}
//		for( itr = nonZeroComps[iThread].begin(); itr != nonZeroComps[iThread].end(); ++itr ){// Inner product
//			std::cout << "iThread iblk nonZeroComps " << iThread << " " << iblk << " " << *itr << std::endl;
//		}
//#endif

		imdl = ptrResistivityBlock->getModelIDFromBlockID(iblk);
		for( ivec = 0; ivec < numInterpolatorVectors; ++ivec ){
			work = std::complex<double>(0.0,0.0);
			//icount = 0;
			itrEnd = nonZeroComps[iThread].end();
			for( itr = nonZeroComps[iThread].begin(); itr != itrEnd; ++itr ){// Inner product
				index = static_cast<long long>(ivec) * static_cast<long long>(numOfEquationFinallySolved) + static_cast<long long>(*itr);
				work += solutionForInterpolatorVectors[index] * nonZeroValues[iThread][*itr];
			}
			index = static_cast<long long>(imdl) + static_cast<long long>(ivec) * static_cast<long long>(nBlkNotFixed);
			derivatives[index] = work;
		}

	}
	//--------------------- End of omp for >>>>>

	}
	//--------------------- End of omp parallel >>>>>

	delete [] solutionAfterDegenerated;
	for( int iThread = 0; iThread < numThreads; ++iThread ){
		delete [] nonZeroValues[iThread];
	}
	delete [] nonZeroValues;
	delete [] nonZeroComps;
	
	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int ivec = 0; ivec < numInterpolatorVectors; ++ivec ){
		for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
			std::cout << "imdl ivec iblk+ivec*nModel derivatives " << imdl << " " << ivec << " " << imdl + ivec * nBlkNotFixed << " " << derivatives[ imdl + ivec * nBlkNotFixed ] << std::endl;
		}
	}
	std::cout << "---" << std::endl;
#endif
	//----- debug <<<<<

}

//// Allocate memory for derivatives of interpolator vectors
//void Forward3D::allcateMemoryForDerivativeOfInterpolatorVectors( const int numInterpolatorVectors ){
//
//	if( m_derivativeOfInterpolatorVectors != NULL ){
//		delete [] m_derivativeOfInterpolatorVectors;
//		m_derivativeOfInterpolatorVectors = NULL;
//	}
//
//	m_derivativeOfInterpolatorVectors = new ComplexSparseSquareMatrix( numInterpolatorVectors * m_numOfEquationDegenerated );
//
//}

//// Release memory of right-hand sides matrix consisting of interpolator vectors
//void Forward3D::releaseMemoryOfRhsVectors(){
//
//	//m_numColumnsRhsMatrixConsistingOfInterpolatorVectors = 0;
//	//m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors = 0;
//
//	if( m_rhsMatrixConsistingOfInterpolatorVectors != NULL ){
//		delete [] m_rhsMatrixConsistingOfInterpolatorVectors;
//		m_rhsMatrixConsistingOfInterpolatorVectors = NULL;
//	}
//
//}

//// Check number of columns of right-hand sides matrix consisting of interpolator vectors
//void Forward3D::checkNumberOfColumnsRhsMatrixConsistingOfInterpolatorVectors(){
//
//	if( m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors != m_numColumnsRhsMatrixConsistingOfInterpolatorVectors ){
//		OutputFiles::m_logFile << "Error : Number of counter does not match the number of columns !! m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors = " << m_counterOfColumnsNumberRhsMatrixConsistingOfInterpolatorVectors << " , m_numColumnsRhsMatrixConsistingOfInterpolatorVectors = " << m_numColumnsRhsMatrixConsistingOfInterpolatorVectors << m_numOfEquationDegenerated << std::endl;
//		exit(1);
//	}
//
//}

// Copy solution vector degenerated
void Forward3D::copySolutionVectorDegenerated( const int iPol, std::complex<double>* solutionVector ) const{

	for( int i = 0 ; i < m_numOfEquation; ++i ){
		if( m_IDsGlobal2AfterDegenerated[iPol][i] < 0 ){
			continue;
		}
		solutionVector[ m_IDsGlobal2AfterDegenerated[iPol][i] ] = m_solution[i];
	}

}

// Get polarization at present for which forward analysis is executed 
int Forward3D::getPolarizationCurrent() const{
	return m_polarizationCurrent;
}

// Get frequency at present for which forward analysis is executed 
double Forward3D::getFrequencyCurrent() const{
	return m_frequencyCurrent;
}

// Get order of finite element
int Forward3D::getOrderOfFiniteElement() const{
	return m_orderOfFiniteElement;
}

// Get total number of equations after degeneration
int Forward3D::getNumOfEquationDegenerated() const{
	return m_numOfEquationDegenerated;
}

// Get total number of equations finally solved
int Forward3D::getNumOfEquationFinallySolved() const{
	return getNumOfEquationDegenerated();
}

// Release memory of coefficient matrix and sparse solver
void Forward3D::releaseMemoryOfMatrixAndSolver(){
	m_matrix3DAnalysis.releaseMemory();
	m_hasMatrixStructureSetAndAnalyzed = false;
}

// Initialize sparse solver
void Forward3D::initializeSparseSolver(){

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	const int myPE = pAnalysisControl->getMyPE();
	const int imode = pAnalysisControl->getModeOfPARDISO();
	std::ostringstream oocHeaderName;
	oocHeaderName << "ooc_temp_3D_PE" << myPE;
	m_matrix3DAnalysis.initializeMatrixSolver( oocHeaderName.str(), imode );

}
