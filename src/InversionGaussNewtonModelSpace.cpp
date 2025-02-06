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
#include "ObservedData.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "ResistivityBlock.h"
#include "InversionGaussNewtonModelSpace.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mkl_lapacke.h"
#include "mpi.h"

// Default constructer
InversionGaussNewtonModelSpace::InversionGaussNewtonModelSpace():
	Inversion()
{}

// Constructer
InversionGaussNewtonModelSpace::InversionGaussNewtonModelSpace( const int nModel, const int nData ):
	Inversion(nModel, nData)
{}

// Destructer
InversionGaussNewtonModelSpace::~InversionGaussNewtonModelSpace(){
}

// Perform inversion
void InversionGaussNewtonModelSpace::inversionCalculation(){

	// Get process ID and total process number
	//int myProcessID(0);
	//MPI_Comm_rank ( MPI_COMM_WORLD, &myProcessID );
	//int numProcessTotal(0);
	//MPI_Comm_size ( MPI_COMM_WORLD, &numProcessTotal );
	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const int myProcessID = ptrAnalysisControl->getMyPE();
	const int numProcessTotal = ptrAnalysisControl->getTotalPE();

	ObservedData* const ptrObservedData = ObservedData::getInstance();

	//int numDataTotalThisPE = ptrObservedData->getNumObservedDataThisPETotal();
	int numDataThisPE = ptrObservedData->getNumObservedDataThisPETotal();
	int* numDataLocal= new int[numProcessTotal];

#ifdef _DEBUG_WRITE
	for( int i = 0; i < numProcessTotal; ++i ){
		std::cout << "PE numDataThisPE : " << myProcessID << " " << numDataThisPE << std::endl;
	}
#endif

	MPI_Allgather( &numDataThisPE, 1, MPI_INT, numDataLocal, 1, MPI_INT, MPI_COMM_WORLD );

#ifdef _DEBUG_WRITE
	for( int i = 0; i < numProcessTotal; ++i ){
		std::cout << "PE i numDataLocal[i] : " << myProcessID << " " << i << " " << numDataLocal[i] << std::endl;
	}
#endif

	//MPI_Reduce( &numDataThisPE, &numDataTotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	int* displacements = new int[ numProcessTotal + 1 ];
	displacements[0] = 0;
	for( int i = 0; i < numProcessTotal; ++i ){
		displacements[i+1] = displacements[i] + numDataLocal[i];
	}
	const int numDataTotal = displacements[numProcessTotal];

#ifdef _DEBUG_WRITE
	for( int i = 0; i < numProcessTotal + 1; ++i ){
		std::cout << "PE i displacements[i] : " << myProcessID << " " << i << " " << displacements[i] << std::endl;
	}
	std::cout << "PE numDataTotal : " << myProcessID << " " << numDataTotal << std::endl;
#endif

	//-----------------------------------------------------
	// Calculate residual vector
	//-----------------------------------------------------

	// Calculate residual vector of this PE
	double* dataVectorThisPE = new double[numDataThisPE];
	for( int i = 0; i < numDataThisPE; ++i ){
		dataVectorThisPE[i] = 0.0; //Initialize
	}
	ptrObservedData->calculateResidualVectorOfDataThisPE( dataVectorThisPE );

#ifdef _DEBUG_WRITE
	std::cout << "Residual vector of data of this PE" << std::endl;
	for( int i = 0; i < numDataThisPE; ++i ){
		std::cout << "PE i dataVectorThisPE[i] " << myProcessID << " " << i << " " << dataVectorThisPE[i] << std::endl;
	}
#endif

	double* dataVectorTotal(NULL);
	if( myProcessID == 0 ){// If this PE number is zero
		dataVectorTotal = new double[numDataTotal];
	}
	MPI_Gatherv( dataVectorThisPE, numDataThisPE, MPI_DOUBLE, dataVectorTotal, numDataLocal, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD );

#ifdef _DEBUG_WRITE
	if( myProcessID == 0 ){// If this PE number is zero
		std::cout << "Residual vector of total data" << std::endl;
		for( int i = 0; i < numDataTotal; ++i ){
			std::cout << "PE  i dataVectorTotal[i] " << myProcessID << " " << i << " " << dataVectorTotal[i] << std::endl;
		}
	}
#endif

	delete [] numDataLocal;
	delete [] displacements;
	delete [] dataVectorThisPE;

	ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();
	const int nBlkNotFixed = ptrResistivityBlock->getNumResistivityBlockNotFixed();
	const int numModel = getNumberOfModel();

	double* rhsVector = new double[ numModel ];

	const long long int numModel_64 = static_cast<long long int>(numModel);
	const long long int numDataTotal_64 = static_cast<long long int>(numDataTotal);

	// If this PE number is zero -------------------------------------------------------
	if( myProcessID == 0 ){

		//const int nFreq = ptrObservedData->getNumOfFrequenciesCalculatedByThisPE();
		const int nFreq = ptrObservedData->getTotalNumberOfDifferenetFrequencies();

		//=================================================================
		// Construct matriecs and vectors
		//=================================================================
		OutputFiles::m_logFile << "# Number of model : " << numModel << std::endl;
		OutputFiles::m_logFile << "# Number of data  : " << numDataTotal << std::endl;
		OutputFiles::m_logFile << "# Construct matrix and vector of the equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

		double* sensitivityMatrixBuf = new double[ numDataTotal_64 * numModel_64 ];

		for( long long int i = 0; i < numDataTotal_64 * numModel_64; ++i ){
			sensitivityMatrixBuf[i] =0.0;// Initialize
		}
		
		long long int numDataAccumulated_64(0);
		for( int iFreq = 0; iFreq < nFreq; ++iFreq ){

			//const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);
			const int freqID = iFreq;

			std::ostringstream fileName;
			if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
				fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
				fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
			}
			fileName << "sensMatFreq" << freqID;
			FILE* fp = fopen( fileName.str().c_str(), "rb" );
			if( fp == NULL ){
				OutputFiles::m_logFile << "File open error !! : " << fileName.str() << std::endl;
				exit(1);
			}

			//const int numDataThisFreq = ptrObservedData->getNumObservedDataThisPE( iFreq );
			//const int numDataAccumulated = ptrObservedData->getNumObservedDataThisPEAccumulated( iFreq );

			int numDataThisFreq(0);
			int numModelTemp(0);
			fread( &numDataThisFreq, sizeof(int), 1, fp );
			fread( &numModelTemp, sizeof(int), 1, fp );
			if( numModel != numModelTemp ){
				OutputFiles::m_logFile << "Error : numModel is not equal to numModelTemp. numModel = " <<  numModel << ", numModelTemp = " << numModelTemp << std::endl;
				exit(1);
			}
			//if( numDataThisFreq != numDataTemp ){
			//	OutputFiles::m_logFile << "Error : numDataThisFreq is not equal to numDataTemp. numDataThisFreq = " <<  numDataThisFreq << ", numDataTemp = " << numDataTemp << std::endl;
			//	exit(1);
			//}

			const long long int numDataThisFreq_64 = static_cast<long long int>(numDataThisFreq);
			fread( &sensitivityMatrixBuf[ numDataAccumulated_64 * numModel_64 ], sizeof(double), numDataThisFreq_64 * numModel_64, fp );
	
			fclose( fp );

			numDataAccumulated_64 += numDataThisFreq_64;// Add data number
		}

		double* sensitivityMatrix = new double[ numDataTotal_64 * numModel_64 ];

		OutputFiles::m_logFile << "# Convert to column-major. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
		//--------------------------------------------------
		// Transpose
		{
			long long int iMdl(0);
			long long int iDat(0);
			long long int offset(0);
#ifdef _USE_OMP
			#pragma omp parallel for default(shared) private( iMdl, iDat, offset )
#endif
			for( iMdl = 0; iMdl < numModel_64; ++iMdl ){
				offset = iMdl * numDataTotal_64;
				for( iDat = 0; iDat < numDataTotal_64; ++iDat ){
					sensitivityMatrix[ offset + iDat ] = sensitivityMatrixBuf[ iMdl + iDat * numModel_64 ];
				}
			}
		}
		//--------------------------------------------------

		delete [] sensitivityMatrixBuf;
	
		const long long int numElemsOfCoefficientMatrix = numModel_64 * ( numModel_64 + 1 ) / 2;
		OutputFiles::m_logFile << "# Total number of elements in coefficient matrix : " << numElemsOfCoefficientMatrix << std::endl;
		double* matrixToBeInverted = new double[numElemsOfCoefficientMatrix];

		OutputFiles::m_logFile << "# Calculate coefficient matrix. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
		//--------------------------------------------------
		// Matrix-Matrix product
		{
			long long int iMdl1(0);
			long long int iMdl2(0);
			long long int numLeft(0);
			long long int offset1(0);
			long long int offset2(0);
			double work(0.0);
			long long int iDat(0);
#ifdef _USE_OMP
			#pragma omp parallel for default(shared) private( iMdl1, iMdl2, numLeft, offset1, offset2, work, iDat ) schedule(static,1)
#endif
			for( iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				offset1 = iMdl1 * numDataTotal_64;
				numLeft = ( iMdl1 + 1 ) * iMdl1 / 2;
				for( iMdl2 = 0; iMdl2 <= iMdl1; ++iMdl2 ){// Left side
					offset2 = iMdl2 * numDataTotal_64;
					work = 0.0;
					for( iDat = 0; iDat < numDataTotal_64; ++iDat ){
						work += sensitivityMatrix[ offset1 + iDat ] * sensitivityMatrix[ offset2 + iDat ];
					}
					matrixToBeInverted[ iMdl2 + numLeft ] = work;
				}
			}
		}
		//--------------------------------------------------

		//---------------------------------
		// Calculate constraining matrix
		//---------------------------------
		RougheningSquareMatrix constrainingMatrix;
		calcConstrainingMatrix( constrainingMatrix );
		
#ifdef _DEBUG_WRITE
		{
			double* temp = new double[numModel_64*numModel_64];
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				for( long long int iMdl2 = 0; iMdl2 < numModel_64; ++iMdl2 ){// Left side
					temp[ iMdl2 + iMdl1 * numModel_64 ] = 0.0;
				}
			}
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				const long long int numLeft = ( iMdl1 + 1 ) * iMdl1 / 2;
				for( long long int iMdl2 = 0; iMdl2 <= iMdl1; ++iMdl2 ){// Left side
					temp[ iMdl2 + iMdl1 * numModel_64 ] = matrixToBeInverted[ iMdl2 + numLeft ];
				}
			}
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				for( long long int iMdl2 = 0; iMdl2 < numModel_64; ++iMdl2 ){// Left side
					std::cout << " " << temp[ iMdl2 + iMdl1 * numModel_64 ];
				}
				std::cout << std::endl;
			}
			delete [] temp;
		}
#endif
		//Make [R]T[R] matrix, where [R] is a constraining matrix
		DoubleSparseSquareSymmetricMatrix RTRMatrix;
		constrainingMatrix.makeRTRMatrix( RTRMatrix );
		const int nRow = RTRMatrix.getNumRows();
		for( int iRow = 0; iRow < nRow; ++iRow ){
			const int nonZeroEnd = RTRMatrix.getRowIndexCRS(iRow+1);
			for( int iNonZero = RTRMatrix.getRowIndexCRS(iRow); iNonZero < nonZeroEnd; ++iNonZero ){
				const int iCol = RTRMatrix.getColumnsCRS(iNonZero);
				const int numLeft = ( iCol + 1 ) * iCol / 2;
				const double value = RTRMatrix.getValueCRS(iNonZero);
				matrixToBeInverted[ iRow + numLeft ] += value;
			}
		}
#ifdef _DEBUG_WRITE
		{
			double* temp2 = new double[numModel_64*numModel_64];
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				for( long long int iMdl2 = 0; iMdl2 < numModel_64; ++iMdl2 ){// Left side
					temp2[ iMdl2 + iMdl1 * numModel_64 ] = 0.0;
				}
			}
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				const long long int numLeft = ( iMdl1 + 1 ) * iMdl1 / 2;
				for( long long int iMdl2 = 0; iMdl2 <= iMdl1; ++iMdl2 ){// Left side
					temp2[ iMdl2 + iMdl1 * numModel_64 ] = matrixToBeInverted[ iMdl2 + numLeft ];
				}
			}
			for( long long int iMdl1 = 0; iMdl1 < numModel_64; ++iMdl1 ){// Right side
				for( long long int iMdl2 = 0; iMdl2 < numModel_64; ++iMdl2 ){// Left side
					std::cout << " " << temp2[ iMdl2 + iMdl1 * numModel_64 ];
				}
				std::cout << std::endl;
			}
			delete [] temp2;
		}
#endif

		//=================================================
		// Calculate right-hand side vector
		//=================================================
		OutputFiles::m_logFile << "# Calculate right-hand side vector. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

		//----------------------------------------------------------------------------------------
		// Calculate matrix vector product of constraining matrix and vector of model roughness 
		//----------------------------------------------------------------------------------------

		ptrResistivityBlock->copyResistivityValuesNotFixedToVectorLog10( rhsVector );
		ptrObservedData->copyDistortionParamsNotFixedToVector( &rhsVector[nBlkNotFixed] );

#ifdef _DEBUG_WRITE
		for( int i = 0; i < numModel; ++i ){
			std::cout << rhsVector[i] << std::endl;
		}
		std::cout << "---" << std::endl;
#endif
		double* workVector = new double[numModel];
		constrainingMatrix.calcVectorOfModelRoughness( rhsVector, workVector );
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numModel; ++i ){
			std::cout << workVector[i] << std::endl;
		}
		std::cout << "---" << std::endl;
#endif
		constrainingMatrix.calcMatrixVectorProductUsingTransposedMatrix( workVector, rhsVector );
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numModel; ++i ){
			std::cout << rhsVector[i] << std::endl;
		}
		std::cout << "---" << std::endl;
#endif
		delete [] workVector;
		//----------------------------------------------------------------------------------------
		// Calculate matrix vector product of transposed sensitivity matrix and residual vector
		//----------------------------------------------------------------------------------------

		//--------------------------------------------------
		// Matrix-vector product
		{
			long long int iMdl(0);
			long long int offset(0);
			double work(0.0);
			long long int iDat(0);
#ifdef _USE_OMP
			#pragma omp parallel for default(shared) private( iMdl, offset, work, iDat )
#endif
			for( iMdl = 0; iMdl < numModel; ++iMdl ){
				offset = iMdl * numDataTotal_64;
				work = 0.0;
				for( iDat = 0; iDat < numDataTotal_64; ++iDat ){
					 work += sensitivityMatrix[ offset + iDat ] * dataVectorTotal[iDat];
				}
				rhsVector[iMdl] += work;
			}
		}
		//--------------------------------------------------

		//=================================================================
		// Solver linear equation with lapack
		//=================================================================

		// Numerical factorization with lapack
		OutputFiles::m_logFile << "# Start numerical factorization for transformed normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

		const bool positiveDefinite = ( AnalysisControl::getInstance() )->getPositiveDefiniteNormalEqMatrix();
		long long int* ipiv = NULL;
		if( !positiveDefinite ){
			ipiv = new long long int[numModel];
		}

		long long int ierr(0);
		if( positiveDefinite ){
			ierr = LAPACKE_dpptrf( LAPACK_COL_MAJOR, 'U', numModel_64, matrixToBeInverted );
		}
		else{
			ierr = LAPACKE_dsptrf( LAPACK_COL_MAJOR, 'U', numModel_64, matrixToBeInverted, ipiv );
		}

		if( ierr > 0 ) {
			OutputFiles::m_logFile << "Error : Matrix is singular. ierr = " << ierr << std::endl;
			exit(1);
		}else if( ierr < 0 ){
			OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
			exit(1);
		}

		// Solver linear equation with lapack
		OutputFiles::m_logFile << "# Start solve phase for transformed normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
		const long long int nrhs = 1;
		const long long int ldb = numModel_64;
		if( positiveDefinite ){
			ierr = LAPACKE_dpptrs( LAPACK_COL_MAJOR, 'U', numModel_64, nrhs, matrixToBeInverted, rhsVector, ldb );
		}
		else{
			ierr = LAPACKE_dsptrs( LAPACK_COL_MAJOR, 'U', numModel_64, nrhs, matrixToBeInverted, ipiv, rhsVector, ldb );
		}

		if( ierr < 0 ){
			OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
			exit(1);
		}

		if( !positiveDefinite ){
			delete [] ipiv;
		}

		// Refine the solution and estimate the error
		//----- modified by Y.Usui 2013.12.21 : Not remove for future use >>>>>
		//OutputFiles::m_logFile << "# Refine the solution. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
		//const int idx = numModel;
		//double* ferr = new double[nrhs];
		//double* berr = new double[nrhs];
		//ierr = LAPACKE_dsprfs( LAPACK_COL_MAJOR, 'U', numModel, nrhs,  matrixToBeInvertedOrg, matrixToBeInverted, ipiv, rhsVectorOrg, ldb, rhsVector, idx, ferr, berr );
		//if( ierr < 0 ){
		//	OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		//	exit(1);
		//}

//#ifdef _DEBUG_WRITE
//		std::cout << "Component-wise forward and backward errors" << std::endl;
//
//		for( int i = 0; i < numModel; ++i ){
//			std::cout << " i ferr[i] : " << i << " " << ferr[i] << std::endl;
//		}
//		for( int i = 0; i < numModel; ++i ){
//			std::cout << " i berr[i] : " << i << " " << berr[i] << std::endl;
//		}
//
//#endif
//
//		delete [] ferr;
//		delete [] berr;
//		delete [] ipiv;
		//----- modified by Y.Usui 2013.12.21 : Not remove for future use <<<<<

		//=================================================================
		// Release memory
		//=================================================================
		delete[] sensitivityMatrix;
		delete[] matrixToBeInverted;	
	}
	// If this PE number is zero -------------------------------------------------------

	MPI_Bcast( rhsVector, numModel, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	//=================================================================
	// Update resistivity values
	//=================================================================
	ptrResistivityBlock->calctResistivityUpdatedFullFromLog10ResistivityIncres( rhsVector );
	ptrResistivityBlock->updateResistivityValues();

	//=================================================================
	// Update distortion parameters
	//=================================================================
	ptrObservedData->calcDistortionParamsUpdatedFullFromIncrements( &rhsVector[nBlkNotFixed] );
	ptrObservedData->updateDistortionParams();

	//=================================================================
	// Release memory
	//=================================================================
	if( myProcessID == 0 ){// If this PE number is zero
		delete [] dataVectorTotal;
	}

	delete[] rhsVector;

	//=================================================================
	// Synchronize
	//=================================================================
	MPI_Barrier( MPI_COMM_WORLD );

	//----- Not delete for future use >>>>>
//	double** sensitivityMatrix = new double*[ nFreq ];
//	double** sensitivityMatrixMultiplied = new double*[ nFreq ];
//	int* numDataAccumulated = new int[ nFreq + 1 ];
//	numDataAccumulated[0] = 0;
//
//	for( int i = 0; i < nFreq; ++i ){
//		sensitivityMatrix[i] = NULL;// Initialize
//		sensitivityMatrixMultiplied[i] = NULL;// Initialize
//	}
//		
//	for( int iFreq = 0; iFreq < nFreq; ++iFreq ){
//
//		//const int freqID = ptrAnalysisControl->getIDsOfFrequenciesCalculatedByThisPE(iFreq);
//		const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);
//
//		std::ostringstream fileName;
//		fileName << "sensMatFreq" << freqID;
//		FILE* fp = fopen( fileName.str().c_str(), "rb" );
//		if( fp == NULL ){
//			OutputFiles::m_logFile << "File open error !! : " << fileName.str() << std::endl;
//			exit(1);
//		}
//
//		int numData(0);
//		int ntemp(0);
//		fread( &numData, sizeof(int), 1, fp );
//		fread( &ntemp, sizeof(int), 1, fp );
//		if( numModel != ntemp ){
//			OutputFiles::m_logFile << "Error : numModel is not equal to ntemp. numModel = " <<  numModel << ", is ntemp = " << ntemp << std::endl;
//			exit(1);
//		}
//
//		sensitivityMatrix[iFreq] = new double[ numData * numModel ];
//		sensitivityMatrixMultiplied[iFreq] = new double[ numData * numModel ];
//		fread( sensitivityMatrix[iFreq], sizeof(double), numData * numModel, fp );
//		fclose( fp );
//
//		memcpy( sensitivityMatrixMultiplied[iFreq], sensitivityMatrix[iFreq], sizeof(double)*(numData * numModel) );
//		numDataAccumulated[iFreq + 1] = numData + numDataAccumulated[iFreq];
//
//#ifdef _DEBUG_WRITE
//		std::cout << "read sensitivity matrix at inversionCalculation. iFreq = " << iFreq << std::endl;
//		std::cout << "sensitivityMatrix" << std::endl;
//		for( int idat = 0; idat < numData; ++idat ){
//			for( int imdl = 0; imdl < numModel; ++imdl ){
//				std::cout << "idat, imdl, val " << idat << " " <<  imdl << " " << sensitivityMatrix[iFreq][imdl+idat*numModel] << std::endl;
//			}
//		}
//		std::cout << "sensitivityMatrixMultiplied" << std::endl;
//		for( int idat = 0; idat < numData; ++idat ){
//			for( int imdl = 0; imdl < numModel; ++imdl ){
//				std::cout << "idat, imdl, val " << idat << " " <<  imdl << " " << sensitivityMatrixMultiplied[iFreq][imdl+idat*numModel] << std::endl;
//			}
//		}
//#endif
//
//		m_rougheningMatrix.multiplyInverseOfTransposedDifferenceFilter( numData, sensitivityMatrixMultiplied[iFreq] );
//
//#ifdef _DEBUG_WRITE
//		std::cout << "read after multiplinginverse matrix of transpose of difference filter . iFreq = " << iFreq << std::endl;
//		for( int idat = 0; idat < numData; ++idat ){
//			for( int imdl = 0; imdl < numModel; ++imdl ){
//				std::cout << "idat, imdl, val " << idat << " " <<  imdl << " " << sensitivityMatrixMultiplied[iFreq][imdl+idat*numModel] << std::endl;
//			}
//		}
//#endif
//
//		m_rougheningMatrix.multiplyInverseOfDifferenceFilter( numData, sensitivityMatrixMultiplied[iFreq] );
//
//#ifdef _DEBUG_WRITE
//		std::cout << "read after multipling inverse matrix of difference filter . iFreq = " << iFreq << std::endl;
//		for( int idat = 0; idat < numData; ++idat ){
//			for( int imdl = 0; imdl < numModel; ++imdl ){
//				std::cout << "idat, imdl, val " << idat << " " <<  imdl << " " << sensitivityMatrixMultiplied[iFreq][imdl+idat*numModel] << std::endl;
//			}
//		}
//#endif
//
//	}
//
//#ifdef _DEBUG_WRITE
//	for( int iFreq = 0; iFreq < nFreq + 1; ++iFreq ){
//		std::cout << "numDataAccumulated[" << iFreq << "] = " << numDataAccumulated[iFreq] << std::endl;
//	}
//#endif
//
//	const int numDataTotalThisPE = numDataAccumulated[nFreq];
//	
//	double* const matrixToBeInverted = new double[ numDataTotalThisPE * ( numDataTotalThisPE + 1 ) / 2 ];
//
//	for( int iFreq1 = 0; iFreq1 < nFreq; ++iFreq1 ){
//		
//		const int nDat1 = numDataAccumulated[iFreq1+1] - numDataAccumulated[iFreq1];
//
//		// Different frequencies
//		for( int iFreq2 = 0; iFreq2 < iFreq1; ++iFreq2 ){
//
//			const int nDat2 = numDataAccumulated[iFreq2+1] - numDataAccumulated[iFreq2];
//
//			for( int iDat1 = 0; iDat1 < nDat1; ++iDat1 ){
//				for( int iDat2 = 0; iDat2 < nDat2; ++iDat2 ){
//
//					const int offset1 = iDat1 * numModel;
//					const int offset2 = iDat2 * numModel;
//
//					double work(0.0);
//					for( int iMdl = 0; iMdl < numModel; ++iMdl ){
//						work += sensitivityMatrixMultiplied[iFreq1][iMdl+offset1] * sensitivityMatrix[iFreq2][iMdl+offset2];
//					}
//
//					const int col = numDataAccumulated[iFreq1] + iDat1;
//					const int row = numDataAccumulated[iFreq2] + iDat2; 
//					matrixToBeInverted[ row + ( col + 1 ) * col / 2 ] = work;
//
//				}
//			}
//
//		}
//
//		// Same frequencies ( diagonal part )
//		for( int iDat1 = 0; iDat1 < nDat1; ++iDat1 ){
//			for( int iDat2 = 0; iDat2 <= iDat1; ++iDat2 ){
//
//				const int offset1 = iDat1 * numModel;
//				const int offset2 = iDat2 * numModel;
//
//				double work(0.0);
//				for( int iMdl = 0; iMdl < numModel; ++iMdl ){
//					work += sensitivityMatrixMultiplied[iFreq1][iMdl+offset1] * sensitivityMatrix[iFreq1][iMdl+offset2];
//				}
//
//				const int col = numDataAccumulated[iFreq1] + iDat1;
//				const int row = numDataAccumulated[iFreq1] + iDat2; 
//				matrixToBeInverted[ row + ( col + 1 ) * col / 2 ] = work;
//								
//			}
//		}
//
//	}
//
//#ifdef _DEBUG_WRITE
//	for( int icol = 0; icol < numDataTotalThisPE; ++icol ){		
//		for( int irow = 0; irow <= icol; ++irow ){
//			std::cout << "row col val " << irow << " " << icol << " " << matrixToBeInverted[ irow + ( icol + 1 ) * icol / 2 ] << std::endl;
//		}
//	}
//#endif
//
//	//=================================================
//	// Calculate data vector of this PE
//	//=================================================
//	double* dataVectorThisPE = new double[numDataTotalThisPE];
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		dataVectorThisPE[i] = 0.0; //Initialize
//	}
//	// Calculate residual vector of data
//	ptrObservedData->calculateResidualVectorOfDataThisPE( dataVectorThisPE );
//
//#ifdef _DEBUG_WRITE
//	std::cout << "Residual vector of datar of this PE" << std::endl;
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << "i dataVectorThisPE[i] " << i << " " << dataVectorThisPE[i] << std::endl;
//	}
//#endif
//
//	// Matrix-vector product
//	for( int iFreq = 0; iFreq < nFreq; ++iFreq ){
//
//		const int nDat = numDataAccumulated[iFreq+1] - numDataAccumulated[iFreq];
//
//		for( int iDat = 0; iDat < nDat; ++iDat ){
//
//			const int offset = iDat * numModel;
//
//			double work(0.0);
//			for( int iMdl = 0; iMdl < numModel; ++iMdl ){
//				work += sensitivityMatrix[iFreq][iMdl+offset] * log( ptrResistivityBlock->getResistivityValuesFromBlockID( ptrResistivityBlock->getBlockIDFromModelID( iMdl ) ) );
//			}
//
//#ifdef _DEBUG_WRITE
//			std::cout << "iDat+numDataAccumulated[iFreq] work " << iDat+numDataAccumulated[iFreq] << " " << work << std::endl;
//#endif
//			
//			dataVectorThisPE[iDat+numDataAccumulated[iFreq]] += work;
//
//		}
//
//	}
//
//
//#ifdef _DEBUG_WRITE
//	std::cout << "Data vector of this PE" << std::endl;
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << "i dataVectorThisPE[i] " << i << " " << dataVectorThisPE[i] << std::endl;
//	}
//#endif
//
//	// Calculate right-hand side vector
//	double* rhsVectorThisPE = new double[numDataTotalThisPE];
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		rhsVectorThisPE[i] = 0.0; //Initialize
//	}
//
//	// Matrix-vector product
//	for( int iRow = 0; iRow < numDataTotalThisPE; ++iRow ){
//
//		const int offset = ( iRow + 1 ) * iRow / 2;
//		double work(0.0);
//
//		// Lower triangle part of the transposed matrix
//		for( int iCol = 0; iCol < iRow; ++iCol ){
//			work += matrixToBeInverted[ offset + iCol ] * dataVectorThisPE[iCol];
//		}
//
//		// Upper triangle part of the transposed matrix
//		for( int iCol = iRow; iCol < numDataTotalThisPE; ++iCol ){
//			work += matrixToBeInverted[ iRow + ( iCol + 1 ) * iCol / 2 ] * dataVectorThisPE[iCol];
//		}
//
//		rhsVectorThisPE[iRow] = work;
//	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "Right-hand side vector of this PE" << std::endl;
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << "i rhsVectorThisPE[i] " << i << " " << rhsVectorThisPE[i] << std::endl;
//	}
//#endif
//
//	// Add square of damping factor to diagonal componets
//	const double dampingFactor = ptrAnalysisControl->getTradeOffParameterForResistivityValue();
//	const double dampingFactorSquared = dampingFactor * dampingFactor;
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		matrixToBeInverted[ ( i + 3 ) * i / 2 ] += dampingFactorSquared;
//	}
//
//
//#ifdef _DEBUG_WRITE
//	std::cout << "After adding square of damping factor to diagonal componets" << std::endl;
//	for( int icol = 0; icol < numDataTotalThisPE; ++icol ){		
//		for( int irow = 0; irow <= icol; ++irow ){
//			std::cout << "row col val " << irow << " " << icol << " " << matrixToBeInverted[ irow + ( icol + 1 ) * icol / 2 ] << std::endl;
//		}
//	}
//#endif
//
//
//	//=================================================================
//	// Solver linear equation with lapack
//	//=================================================================
////#ifdef _DEBUG_WRITE
////	double a[15] = {3.14, 0.17, 0.79, -0.9, 0.83, 4.53, 1.65, -0.65, -3.7, 5.32, -0.72, 0.28, 1.6, -1.37, 1.98};
////	const int num = 5;
////	const int itmp = LAPACKE_dpptrf( LAPACK_COL_MAJOR, 'U', num, a );
////	if( itmp > 0 ) {
////		OutputFiles::m_logFile << "Error : Test matrix is not positive definite. info = " << itmp << std::endl;
////		exit(1);
////	}else if( itmp < 0 ){
////		OutputFiles::m_logFile << "Error : " << -itmp << "-th parameter has illegal value." << std::endl;
////		exit(1);
////	}
////	for( int i = 0; i < 15; ++i ){
////		std::cout << " i val : " << i << " " << a[i] << std::endl;
////	}
////#endif
//
//	// Numerical factorization with lapack
//	OutputFiles::m_logFile << "# Start numerical factorization for transformed normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
//
//	//const int info = LAPACKE_dpptrf( LAPACK_COL_MAJOR, 'U', numDataTotalThisPE, matrixToBeInverted );
//	//if( info > 0 ) {
//	//	OutputFiles::m_logFile << "Error : Matrix is not positive definite. info = " << info << std::endl;
//	//	exit(1);
//	//}else if( info < 0 ){
//	//	OutputFiles::m_logFile << "Error : " << -info << "-th parameter has illegal value." << std::endl;
//	//	exit(1);
//	//}
//
//	// Copy matrix and vector for refinement of the solution and error estimatation
//	double* matrixToBeInvertedOrg = new double[ numDataTotalThisPE * ( numDataTotalThisPE + 1 ) / 2 ];
//	memcpy( matrixToBeInvertedOrg, matrixToBeInverted, sizeof(double)*(numDataTotalThisPE * ( numDataTotalThisPE + 1 ) / 2 ) );
//	double* rhsVectorThisPEOrg = new double[ numDataTotalThisPE ];
//	memcpy( rhsVectorThisPEOrg, rhsVectorThisPE, sizeof(double)*(numDataTotalThisPE) );
//
//	int* ipiv = new int[numDataTotalThisPE];
//	int ierr(0);
//	ierr = LAPACKE_dsptrf( LAPACK_COL_MAJOR, 'U', numDataTotalThisPE, matrixToBeInverted, ipiv );
//	if( ierr > 0 ) {
//		OutputFiles::m_logFile << "Error : Matrix is singular. ierr = " << ierr << std::endl;
//		exit(1);
//	}else if( ierr < 0 ){
//		OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
//		exit(1);
//	}
//
//	// Solver linear equation with lapack
//	OutputFiles::m_logFile << "# Start solve phase for transformed normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
//	const int nrhs = 1;
//	const int ldb = numDataTotalThisPE;
//	ierr = LAPACKE_dsptrs( LAPACK_COL_MAJOR, 'U', numDataTotalThisPE, nrhs, matrixToBeInverted, ipiv, rhsVectorThisPE, ldb );
//	if( ierr < 0 ){
//		OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
//		exit(1);
//	}
//
//	// Refine the solution and estimate the error
//	OutputFiles::m_logFile << "# Refine the solution. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
//	const int idx = numDataTotalThisPE;
//	double* ferr = new double[nrhs];
//	double* berr = new double[nrhs];
//	ierr = LAPACKE_dsprfs( LAPACK_COL_MAJOR, 'U', numDataTotalThisPE, nrhs,  matrixToBeInvertedOrg, matrixToBeInverted, ipiv, rhsVectorThisPEOrg, ldb, rhsVectorThisPE, idx, ferr, berr );
//	if( ierr < 0 ){
//		OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
//		exit(1);
//	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "Component-wise forward and backward errors" << std::endl;
//
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << " i ferr[i] : " << i << " " << ferr[i] << std::endl;
//	}
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << " i berr[i] : " << i << " " << berr[i] << std::endl;
//	}
//
//#endif
//
//	delete [] ferr;
//	delete [] berr;
//	delete [] ipiv;
//
//#ifdef _DEBUG_WRITE
//	std::cout << "Solution vector for inversion" << std::endl;
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		std::cout << " i rhsVectorThisPE[i] : " << i << " " << rhsVectorThisPE[i] << std::endl;
//	}
//
//	std::cout << "Check solution vector for inversion" << std::endl;
//	// Matrix-vector product
//	for( int iRow = 0; iRow < numDataTotalThisPE; ++iRow ){
//
//		const int offset = ( iRow + 1 ) * iRow / 2;
//		double work(0.0);
//
//		// Lower triangle part of the transposed matrix
//		for( int iCol = 0; iCol < iRow; ++iCol ){
//			work += matrixToBeInvertedOrg[ offset + iCol ] * rhsVectorThisPE[iCol];
//		}
//
//		// Upper triangle part of the transposed matrix
//		for( int iCol = iRow; iCol < numDataTotalThisPE; ++iCol ){
//			work += matrixToBeInvertedOrg[ iRow + ( iCol + 1 ) * iCol / 2 ] * rhsVectorThisPE[iCol];
//		}
//
//		std::cout << " iRow work : " << iRow << " " << work << std::endl;
//	}
//#endif
//
//	delete [] matrixToBeInvertedOrg;
//	delete [] rhsVectorThisPEOrg;
//
//	//=================================================================
//	// Calculate new model
//	//=================================================================
//	OutputFiles::m_logFile << "# Calculate new model. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
//	double* newModelValues = new double[numModel];
//	//const double dampingFactor = ptrAnalysisControl->getTradeOffParameterForResistivityValue();
//	const double multipliedFactor = 1.0/dampingFactorSquared;
//
//	for( int i = 0; i < numDataTotalThisPE; ++i ){
//		// Subtract solution vector
//		dataVectorThisPE[i] -= rhsVectorThisPE[i];
//	}
//
//	for( int iFreq = 0; iFreq < nFreq; ++iFreq ){
//		
//		const int nDat = numDataAccumulated[iFreq+1] - numDataAccumulated[iFreq];
//
//		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
//
//			double work(0.0);
//			for( int iDat = 0; iDat < nDat; ++iDat ){
//				work += sensitivityMatrixMultiplied[iFreq][ iMdl + iDat * numModel ] * dataVectorThisPE[ numDataAccumulated[iFreq] + iDat ];
//			}
//
//			//newModelValues[iMdl] = work * multipliedFactor;
//			newModelValues[iMdl] = pow( 10.0, work * multipliedFactor );
//		}
//
//	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "New model" << std::endl;
//	for( int iMdl = 0; iMdl < numModel; ++iMdl ){
//		std::cout << " iMdl newModelValues[iMdl] : " << iMdl << " " << newModelValues[iMdl] << std::endl;
//	}
//#endif
//
//	//=================================================================
//	// Update resistivity values
//	//=================================================================
//	ptrResistivityBlock->setResistivityValuesIncreFull( newModelValues );
//
//
//	//=================================================================
//	// Release memory
//	//=================================================================
//	for( int i = 0; i < nFreq; ++i ){
//		delete [] sensitivityMatrix[i];
//	}
//	delete[] sensitivityMatrix;
//
//	for( int i = 0; i < nFreq; ++i ){
//		delete [] sensitivityMatrixMultiplied[i];
//	}
//	delete[] sensitivityMatrixMultiplied;
//
//	delete [] numDataAccumulated;
//
//	delete [] matrixToBeInverted;
//
//	delete [] dataVectorThisPE;
//
//	delete [] rhsVectorThisPE;
//
//	delete [] newModelValues;
//
	//----- Not delete for future use <<<<<

}
