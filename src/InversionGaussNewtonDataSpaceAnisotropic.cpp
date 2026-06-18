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
#include "InversionGaussNewtonDataSpaceAnisotropic.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include "mpi.h"

// Default constructer
InversionGaussNewtonDataSpaceAnisotropic::InversionGaussNewtonDataSpaceAnisotropic():
	Inversion()
{}

// Constructer
InversionGaussNewtonDataSpaceAnisotropic::InversionGaussNewtonDataSpaceAnisotropic( const int nModel, const int nData ):
	Inversion(nModel, nData)
{}

// Destructer
InversionGaussNewtonDataSpaceAnisotropic::~InversionGaussNewtonDataSpaceAnisotropic(){
}

// Perform inversion
void InversionGaussNewtonDataSpaceAnisotropic::inversionCalculation(){

	// Get process ID and total process number
	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const int myProcessID = ptrAnalysisControl->getMyPE();
	const int numProcessTotal = ptrAnalysisControl->getTotalPE();

#ifdef _DEBUG_WRITE
	std::ostringstream oss;
	oss << "debug_" << myProcessID << ".txt";
	std::ofstream fout;
	fout.open( oss.str().c_str() );
#endif

	ObservedData* const ptrObservedData = ObservedData::getInstance();
	ResistivityBlockAnisotropic* const ptrResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	const int nBlkNotFixed = ptrResistivityBlock->getNumberOfUnfixedResistivityParameters();
	const int numModel = getNumberOfModel();
	const int nFreqThisPE = ptrObservedData->getNumOfFrequenciesCalculatedByThisPE();

	OutputFiles::m_logFile << "# Number of model : " << numModel << std::endl;
	OutputFiles::m_logFile << "# Number of frequencies of this PE : " << nFreqThisPE << std::endl;

	//---------------------------------
	// Calculate constraining matrix
	//---------------------------------
	OutputFiles::m_logFile << "# Calculate constraining matrix. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
	RougheningMatrix constrainingMatrixFull;
	RougheningMatrix constrainingMatrixWithoutCrossProductTerm;
	RougheningMatrix derivativeMatrixOfCrossProductTerm;
	std::vector<double> crossProductVector;
	calcConstrainingMatrices(constrainingMatrixFull, constrainingMatrixWithoutCrossProductTerm, derivativeMatrixOfCrossProductTerm, crossProductVector);

	//-----------------------------------------------
	// Calculate vector of model roughness
	//-----------------------------------------------
	// modelVector : m
	// dVector : RTRm = [R]T*[b-[R]*m] - [deriv_G]T*G
	//-----------------------------------------------
	double* dVector(NULL);
	if( myProcessID == 0 ){
		double* modelVector = new double[numModel];
		ptrResistivityBlock->copyUnfixedResistivityParametersPreToVector(modelVector);
		ptrObservedData->copyDistortionParamsNotFixedPreToVector( &modelVector[nBlkNotFixed] );
		const int numRows = constrainingMatrixWithoutCrossProductTerm.getNumRows();
		double* vectorRxm = new double[numRows];
		constrainingMatrixWithoutCrossProductTerm.calcResidualVector( modelVector, vectorRxm );
		delete [] modelVector;
		dVector = new double[numModel];
		constrainingMatrixWithoutCrossProductTerm.calcMatrixVectorProductUsingTransposedMatrix( vectorRxm, dVector );
		delete [] vectorRxm;
		double* dVectorCrossProduct = new double[numModel];
		derivativeMatrixOfCrossProductTerm.calcMatrixVectorProductUsingTransposedMatrix(crossProductVector, dVectorCrossProduct);
		for (int iMdl = 0; iMdl < numModel; ++iMdl) {
			dVector[iMdl] -= dVectorCrossProduct[iMdl];
		}
		delete[] dVectorCrossProduct;
	}

	//----------------------------------------------------------------------------------------
	// Make [R]T[R] matrix, where [R] is a constraining matrix
	//----------------------------------------------------------------------------------------
	OutputFiles::m_logFile << "# Make [R]T[R] matrix." << ptrAnalysisControl->outputElapsedTime() << std::endl;
	DoubleSparseSquareSymmetricMatrix RTRMatrix;
	if( ptrAnalysisControl->isSmallValueToRougheningMatrixDiagonals() ){
		// Add small value to diagonals of  [R]T[R] matrix
		constrainingMatrixFull.makeRTRMatrix( RTRMatrix, ptrAnalysisControl->getSmallValueAddedToDiagonals() );
	}else{
		constrainingMatrixFull.makeRTRMatrix( RTRMatrix );
	}
	constrainingMatrixFull.releaseMemory();

	//----------------------------------
	// Initialization
	//----------------------------------
	std::ostringstream oocHeaderName;
	oocHeaderName << "ooc_temp_inv_3D_PE" << myProcessID;
	RTRMatrix.initializeMatrixSolver( oocHeaderName.str(), ptrAnalysisControl->getModeOfPARDISO() );

	//----------------------------------
	// Analysis
	//----------------------------------
	OutputFiles::m_logFile << "# Analyse constraining matrix." << ptrAnalysisControl->outputElapsedTime() << std::endl;
	RTRMatrix.analysisPhaseMatrixSolver();

	//----------------------------------
	// Factorization
	//----------------------------------
	OutputFiles::m_logFile << "# Factorize constraining matrix." << ptrAnalysisControl->outputElapsedTime() << std::endl;
	RTRMatrix.factorizationPhaseMatrixSolver();

	//------------------------------------------------------------------------------------------------------
	// Calculate residual vector and sensitivity matrix multiplied by inverse of the constraing matrix
	//------------------------------------------------------------------------------------------------------
	int numDataThisPE = ptrObservedData->getNumObservedDataThisPETotal();
	OutputFiles::m_logFile << "# Number of data of this PE : " << numDataThisPE << std::endl;
#ifdef _DEBUG_WRITE
	for( int i = 0; i < numProcessTotal; ++i ){
		fout << "PE numDataThisPE : " << myProcessID << " " << numDataThisPE << std::endl;
	}
#endif

	double* residualVectorThisPE = new double[numDataThisPE];
	for (int i = 0; i < numDataThisPE; ++i) {
		residualVectorThisPE[i] = ptrObservedData->getValueOfResidualVectorComponentPre(i);
	}
#ifdef _DEBUG_WRITE
	for( int i = 0; i < numDataThisPE; ++i ){
		fout << "PE i residualVectorThisPE[i] " << myProcessID << " " << i << " " << residualVectorThisPE[i] << std::endl;
	}
#endif

	double* vectorJTxResidialLocal = new double[numModel];
	for( int iMdl = 0; iMdl < numModel; ++iMdl ){
		vectorJTxResidialLocal[iMdl] = 0.0;// Zero clear
	}
	int offset(0);
	for( int iFreq = 0; iFreq < nFreqThisPE; ++iFreq ){
		//assert( offset == ptrObservedData->getNumObservedDataThisPEAccumulated( iFreq ) );
		const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);

		std::ostringstream fileName;
		if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
		}
		fileName << "sensMatFreq" << freqID;

		OutputFiles::m_logFile << "# Read sensitivity matrix from "<< fileName.str() << "." << ptrAnalysisControl->outputElapsedTime() << std::endl;
		int numDataThisFreq(0);
		int numModelTemp(0);
		double* sensitivityMatrix = NULL;
		readSensitivityMatrix( fileName.str(), numDataThisFreq, numModelTemp, sensitivityMatrix );
		if( numDataThisFreq != ptrObservedData->getNumObservedDataThisPE( iFreq ) ){
			OutputFiles::m_logFile << "Error : Number of data written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}
		if( numModel != numModelTemp ){
			OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}

#ifdef _DEBUG_WRITE
		fout << "read sensitivity matrix at inversionCalculation. iFreq = " << iFreq << std::endl;
		fout << "sensitivityMatrix" << std::endl;
		for( int idat = 0; idat < numDataThisFreq; ++idat ){
			for( int imdl = 0; imdl < numModel; ++imdl ){
				fout << "idat, imdl, val " << idat + offset << " " <<  imdl << " " << sensitivityMatrix[imdl+idat*numModel] << std::endl;
			}
		}
#endif
		//-----------------------------------------------	
		// Matrix-vector product
		//-----------------------------------------------	
		if (CommonParameters::IS_BLAS_USED) {
			CBLAS_ORDER order = CblasRowMajor;
			CBLAS_TRANSPOSE trans = CblasTrans;
			MKL_INT m = static_cast<MKL_INT>(numDataThisFreq);
			MKL_INT n = static_cast<MKL_INT>(numModel);
			double alpha = 1.0;
			double beta = 1.0;		
			MKL_INT lda = n;
			MKL_INT incx = 1;
			MKL_INT incy = 1;
			cblas_dgemv(order, trans, m, n, alpha, sensitivityMatrix, lda, &residualVectorThisPE[offset], incx, beta, vectorJTxResidialLocal, incy);
		}else{
			const long long numModel_64 = static_cast<long long>(numModel);
			const long long offset_64 = static_cast<long long>(offset);
			const long long nDataEnd = offset_64 + static_cast<long long>(numDataThisFreq);
			for( long long iDat = offset_64; iDat < nDataEnd; ++iDat ){
				const double val = residualVectorThisPE[iDat];
				const long long offsetTemp = ( iDat - offset_64 ) * numModel_64;
				for( long long iMdl = 0; iMdl < numModel_64; ++iMdl ){
					vectorJTxResidialLocal[iMdl] += sensitivityMatrix[ iMdl + offsetTemp ] * val;
				}
			}
		}
		//-----------------------------------------------

		delete[] sensitivityMatrix;
		offset += numDataThisFreq;// Add data number
	}

	delete[] residualVectorThisPE;
	if( myProcessID != 0 ){
		// PE = 0 uses constraining matrix later
		RTRMatrix.releaseMemory();
	}

	//-------------------------------------------------
	// Merge vectors
	//-------------------------------------------------
	double* vectorJTxResidialGlobal = NULL;
	if( myProcessID == 0 ){
		vectorJTxResidialGlobal = new double[numModel];
	}
	MPI_Reduce( vectorJTxResidialLocal, vectorJTxResidialGlobal, numModel, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	delete [] vectorJTxResidialLocal;

#ifdef _DEBUG_WRITE
	if( myProcessID == 0 ){
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			fout << "iMdl, vectorJTxResidialGlobal : " << iMdl << " " << vectorJTxResidialGlobal[iMdl] << std::endl;
		}
	}
#endif

	//-------------------------------------------------
	// Calculate Inv([R]T*[R])*d vector
	//-------------------------------------------------
	double* vectorInvRTRd = new double[numModel];
	if( myProcessID == 0 ){//----- Processs ID = 0 Only ----->>>>>
		//------------------------------------------
		// vectorJTxResidialGlobal : [J]T*rd
		// dVector : -RTRm-[deriv_G]T*G => -RTRm-[deriv_G]T*G+[J]T*rd
		//------------------------------------------
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			dVector[iMdl] += vectorJTxResidialGlobal[iMdl];
		}
		delete [] vectorJTxResidialGlobal;

#ifdef _DEBUG_WRITE
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			fout << "iMdl, dVector : " << iMdl << " " << dVector[iMdl] << std::endl;
		}
#endif
		
		//-------------------------------------------------
		// dVector : -RTRm-[deriv_G]T*G+[J]T*rd
		// vectorInvRTRd : inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
		//-------------------------------------------------
		RTRMatrix.solvePhaseMatrixSolver( 1, dVector, vectorInvRTRd );

#ifdef _DEBUG_WRITE
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			fout << "iMdl, vectorInvRTRd : " << iMdl << " " << vectorInvRTRd[iMdl] << std::endl;
		}
#endif
	}//---- Processs ID = 0 Only -----<<<<<
	MPI_Bcast( vectorInvRTRd, numModel, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	//-------------------------------------------------
	// Calculate right-hand side vector
	//-------------------------------------------------
	OutputFiles::m_logFile << "# Calculate right-hand side vector. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
	double* rhsVectorLocal = new double[numDataThisPE];
	for( int i = 0; i < numDataThisPE; ++i ){
		rhsVectorLocal[i] = 0.0;// Zero clear
	}
	offset = 0;
	for( int iFreq = 0; iFreq < nFreqThisPE; ++iFreq ){
		//assert( offset == ptrObservedData->getNumObservedDataThisPEAccumulated( iFreq ) );
		const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);

		std::ostringstream fileName;
		if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
		}
		fileName << "sensMatFreq" << freqID;

		OutputFiles::m_logFile << "# Read sensitivity matrix from "<< fileName.str() << "." << ptrAnalysisControl->outputElapsedTime() << std::endl;
		int numDataThisFreq(0);
		int numModelTemp(0);
		double* sensitivityMatrix(NULL);
		readSensitivityMatrix( fileName.str(), numDataThisFreq, numModelTemp, sensitivityMatrix );

		if( numDataThisFreq != ptrObservedData->getNumObservedDataThisPE( iFreq ) ){
			OutputFiles::m_logFile << "Error : Number of data written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}
		if( numModel != numModelTemp ){
			OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}

		if (CommonParameters::IS_BLAS_USED) {
			CBLAS_ORDER order = CblasRowMajor;
			CBLAS_TRANSPOSE trans = CblasNoTrans;
		    MKL_INT m = static_cast<MKL_INT>(numDataThisFreq);
			MKL_INT n = static_cast<MKL_INT>(numModel);
			double alpha = 1.0;
			double beta = 0.0;
			MKL_INT lda = n;
			MKL_INT incx = 1;
			MKL_INT incy = 1;
			cblas_dgemv(order, trans, m, n, alpha, sensitivityMatrix, lda, vectorInvRTRd, incx, beta, &rhsVectorLocal[offset], incy);
		}else{
			long long iDat(0);
			long long offsetTemp(0);
			double work(0.0);
			long long iMdl(0);
			const long long numDataThisFreq_64 = static_cast<long long>(numDataThisFreq);
			const long long numModel_64 = static_cast<long long>(numModel);
			const long long offset_64 = static_cast<long long>(offset);
#ifdef _USE_OMP
			#pragma omp parallel for default(shared) private( iDat, offsetTemp, work, iMdl )
#endif
			for( iDat = 0; iDat < numDataThisFreq_64; ++iDat ){
				offsetTemp = iDat * numModel_64;
				work = 0.0;
				for( iMdl = 0; iMdl < numModel_64; ++iMdl ){
					work += sensitivityMatrix[ iMdl + offsetTemp ] * vectorInvRTRd[iMdl];
				}
				rhsVectorLocal[offset_64 + iDat] = work;
			}
		}

		delete [] sensitivityMatrix;

		offset += numDataThisFreq;// Add data number
	}
	delete [] vectorInvRTRd;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < numDataThisPE; ++i ){
		fout << "i rhsVectorLocal[i] " << i << " " << rhsVectorLocal[i] << std::endl;
	}
#endif

	//------------------------------------------------------------
	// Make array for garhering
	//------------------------------------------------------------
	int* numDataLocal(NULL);
	if( myProcessID == 0 ){
		numDataLocal= new int[numProcessTotal];
	}

	MPI_Gather( &numDataThisPE, 1, MPI_INT, numDataLocal, 1, MPI_INT, 0, MPI_COMM_WORLD );

	int* numDataAccumulated(NULL);
	int numDataTotal(-1);
	if( myProcessID == 0 ){
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numProcessTotal; ++i ){
			fout  << "PE i numDataLocal[i] : " << myProcessID << " " << i << " " << numDataLocal[i] << std::endl;
		}
#endif
		numDataAccumulated = new int[ numProcessTotal + 1 ];
		numDataAccumulated[0] = 0;
		for( int i = 0; i < numProcessTotal; ++i ){
			numDataAccumulated[i+1] = numDataAccumulated[i] + numDataLocal[i];
		}

		numDataTotal = numDataAccumulated[numProcessTotal];
		OutputFiles::m_logFile << "# Number of total data  : " << numDataTotal << std::endl;

#ifdef _DEBUG_WRITE
		for( int i = 0; i < numProcessTotal + 1; ++i ){
			fout << "PE i numDataAccumulated[i] : " << myProcessID << " " << i << " " << numDataAccumulated[i] << std::endl;
		}
		fout << "PE numDataTotal : " << myProcessID << " " << numDataTotal << std::endl;
#endif
	}

	//------------------------------------------------------------
	// Merge right-hands-side vector
	//------------------------------------------------------------
	double* rhsVectorGlobal(NULL);
	if( myProcessID == 0 ){
		rhsVectorGlobal = new double[numDataTotal];
	}
	MPI_Gatherv( rhsVectorLocal, numDataThisPE, MPI_DOUBLE, rhsVectorGlobal, numDataLocal, numDataAccumulated, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	delete [] rhsVectorLocal;
	
#ifdef _DEBUG_WRITE
	if( myProcessID == 0 ){
		fout << "RHS vector" << std::endl;
		for( int i = 0; i < numDataTotal; ++i ){
			fout << "i rhsVectorGlobal[i] " << i << " " << rhsVectorGlobal[i] << std::endl;
		}
	}
#endif

	//------------------------------------------------------------
	// Calculate coefficient matrix
	//------------------------------------------------------------
	double* matrixToBeInverted(NULL);
	double denominatorOfGCV(0.0);
	if( myProcessID == 0 ){//---- Processs ID = 0 Only ----->>>>>
		OutputFiles::m_logFile << "# Calculate coefficient matrix. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

		const long long int numDataTotal_64 = static_cast<long long int>(numDataTotal);
		double* matrix = new double[numDataTotal_64 * numDataTotal_64];
#ifdef _DEBUG_WRITE
		for( long long i = 0; i < static_cast<long long>(numDataTotal) * static_cast<long long>(numDataTotal); ++i ){
			matrix[i] = -9999.999;
		}
#endif

		//--------------------------------------------------------------
		// Read sensitivity matrix and calculate matrix-matrix product
		//--------------------------------------------------------------
		const int nFreq = ptrObservedData->getTotalNumberOfDifferenetFrequencies();
		int offsetRows(0);
		for( int ifreqLeft = 0 ;ifreqLeft < nFreq; ++ifreqLeft ){
			std::ostringstream fileNameLeft;
			if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
				fileNameLeft << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
				fileNameLeft << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
			}
			fileNameLeft << "sensMatFreq" << ifreqLeft;
			int numDataThisFreqLeft(0);
			int numModelTempLeft(0);
			double* sensitivityMatrixLeft(NULL);
			readSensitivityMatrixMod(RTRMatrix, fileNameLeft.str(), numDataThisFreqLeft, numModelTempLeft, sensitivityMatrixLeft);
			if (numModel != numModelTempLeft) {
				OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
				exit(1);
			}

			int offsetCols = offsetRows;
			for( int ifreqRight = ifreqLeft ;ifreqRight < nFreq; ++ifreqRight ){
				std::ostringstream fileNameRight;
				if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
					fileNameRight << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
					fileNameRight << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
				}
				fileNameRight << "sensMatFreq" << ifreqRight;
				int numDataThisFreqRight(0);
				int numModelTempRight(0);
				double* sensitivityMatrixRight(NULL);
				OutputFiles::m_logFile << "# Read sensitivity matrix from " << fileNameRight.str() << "." << ptrAnalysisControl->outputElapsedTime() << std::endl;
				readSensitivityMatrix(fileNameRight.str(), numDataThisFreqRight, numModelTempRight, sensitivityMatrixRight);
				if (numModel != numModelTempRight) {
					OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
					exit(1);
				}

				//-----------------------------
				//--- Matrix-Matrix product ---
				//-----------------------------
				if (CommonParameters::IS_BLAS_USED) {
					double* result = new double [ static_cast<long long>(numDataThisFreqLeft) * static_cast<long long>(numDataThisFreqRight) ];
					MKL_INT m = static_cast<MKL_INT>(numDataThisFreqLeft);
					MKL_INT n = static_cast<MKL_INT>(numDataThisFreqRight);
					MKL_INT k = static_cast<MKL_INT>(numModel);
					MKL_INT lda = k;
					MKL_INT ldb = k;
					MKL_INT ldc = n;
					double alpha = 1.0;
					double beta = 0.0;
					CBLAS_ORDER  order = CblasRowMajor;
					CBLAS_TRANSPOSE transA = CblasNoTrans;
					CBLAS_TRANSPOSE transB = CblasTrans;
					cblas_dgemm(order, transA, transB, m, n, k, alpha, sensitivityMatrixLeft, lda, sensitivityMatrixRight, ldb, beta, result, ldc);
					for( long long irow = 0; irow < numDataThisFreqLeft; ++irow ){
						const long long row = irow + static_cast<long long>(offsetRows);
						for( long long icol = 0; icol < numDataThisFreqRight; ++icol ){
							const long long col = icol + static_cast<long long>(offsetCols);
							matrix[ row * numDataTotal_64 + col ] = result[ irow * static_cast<long long>(numDataThisFreqRight) + icol ];
						}
					}
					delete [] result;
				}else{
					long long iDatLeft(0);
					long long iDatRight(0);
					long long offsetLeft(0);
					long long offsetRight(0);
					long long index(0);
					double work(0.0);
					long long iMdl(0);
					const long long numDataThisFreqLeft_64 = static_cast<long long>(numDataThisFreqLeft);
					const long long numDataThisFreqRight_64 = static_cast<long long>(numDataThisFreqRight);
					const long long numModel_64 = static_cast<long long>(numModel);
					const long long offsetRows_64 = static_cast<long long>(offsetRows);
					const long long offsetCols_64 = static_cast<long long>(offsetCols);
#ifdef _USE_OMP
					#pragma omp parallel for default(shared) private( iDatLeft, iDatRight, offsetLeft, offsetRight, index, work, iMdl )
#endif
					for( iDatLeft = 0; iDatLeft < numDataThisFreqLeft_64; ++iDatLeft ){
						offsetLeft = iDatLeft * numModel_64;
						for( iDatRight = 0; iDatRight < numDataThisFreqRight_64; ++iDatRight ){
							offsetRight = iDatRight * numModel_64;
							work = 0.0;
							for( iMdl = 0; iMdl < numModel_64; ++iMdl ){
								work += sensitivityMatrixLeft[offsetLeft+iMdl] * sensitivityMatrixRight[offsetRight+iMdl];
							}
							index = (iDatLeft + offsetRows_64) * numDataTotal_64 + (iDatRight + offsetCols_64);
							matrix[index] = work;
						}
					}

				}
				//-----------------------------
				delete [] sensitivityMatrixRight;
				offsetCols += numDataThisFreqRight;// Increment offset of column number
			}
			delete [] sensitivityMatrixLeft;
			offsetRows += numDataThisFreqLeft;// Increment offset of row number
		}

		//----------------------------------------------------------------------
		// Copy this matrix for GCV calculation before the unit matrix is added
		//----------------------------------------------------------------------
		double* matrixK = new double[numDataTotal_64 * numDataTotal_64];
		for (long long int row = 0; row < numDataTotal_64; ++row) {
			matrixK[row * numDataTotal_64 + row] = matrix[row * numDataTotal_64 + row];
			for (long long int col = row + 1; col < numDataTotal_64; ++col) {
				matrixK[row * numDataTotal_64 + col] = matrix[row * numDataTotal_64 + col];
				matrixK[col * numDataTotal_64 + row] = matrix[row * numDataTotal_64 + col];
			}
		}

		//----------------------------------------------
		// Add unit matrix
		//----------------------------------------------
		for( long long int row = 0; row < numDataTotal_64; ++row ){
			const long long int col = row;
			matrix[ row * numDataTotal_64 + col ] += 1.0;
		}

		const long long numElemsOfCoefficientMatrixTotal = numDataTotal_64 * ( numDataTotal_64 + 1 ) / 2;
		OutputFiles::m_logFile << "# Total number of elements in coefficient matrix : " << numElemsOfCoefficientMatrixTotal << std::endl;
		matrixToBeInverted = new double [numElemsOfCoefficientMatrixTotal];

		//----------------------------------------------
		// Copy only upper triangle components
		//----------------------------------------------
		long long int index(0);
		for( long long int row = 0; row < numDataTotal_64; ++row ){
			for(long long int col = row; col < numDataTotal_64; ++col ){
				matrixToBeInverted[index] = matrix[ row * numDataTotal_64 + col ];
				++index;
			}
		}
#ifdef _DEBUG_WRITE
		for( int i = 0; i < numElemsOfCoefficientMatrixTotal; ++i ){
			fout << "i, matrixToBeInverted : " << i << " " << matrixToBeInverted[i] << std::endl;
		}
#endif
		delete [] matrix;

		//----------------------------------------------
		// Numerical factorization with lapack
		//----------------------------------------------
		OutputFiles::m_logFile << "# Start numerical factorization for normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

		const long long int numModel_64 = static_cast<long long int>(numModel);
		const bool positiveDefinite = ptrAnalysisControl->getPositiveDefiniteNormalEqMatrix();
		long long int* ipiv = NULL;
		if( !positiveDefinite ){
			ipiv = new long long int[numDataTotal_64];
		}

		long long int ierr(0);
		if( positiveDefinite ){
			ierr = LAPACKE_dpptrf( LAPACK_COL_MAJOR, 'L', numDataTotal_64, matrixToBeInverted );
		}
		else{
			ierr = LAPACKE_dsptrf( LAPACK_COL_MAJOR, 'L', numDataTotal_64, matrixToBeInverted, ipiv );
		}

		if( ierr > 0 ) {
			OutputFiles::m_logFile << "Error : Matrix is singular. ierr = " << ierr << std::endl;
			exit(1);
		}else if( ierr < 0 ){
			OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
			exit(1);
		}

		//----------------------------------------------
		// Solver linear equation with lapack
		//----------------------------------------------
		OutputFiles::m_logFile << "# Start solve phase for normal equation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
		const long long int nrhs = 1;
		const long long int ldb = numDataTotal_64;
		if( positiveDefinite ){
			ierr = LAPACKE_dpptrs( LAPACK_COL_MAJOR, 'L', numDataTotal_64, nrhs, matrixToBeInverted, rhsVectorGlobal, ldb );
		}
		else{
			ierr = LAPACKE_dsptrs( LAPACK_COL_MAJOR, 'L', numDataTotal_64, nrhs, matrixToBeInverted, ipiv, rhsVectorGlobal, ldb );
		}

		if( ierr < 0 ){
			OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
			exit(1);
		}

		//----------------------------------------------
		// Caluculate the hat matrix for GCV
		//----------------------------------------------
		if (ptrAnalysisControl->needsGCVCalculation()) {
			double* matrixKMod = new double[numDataTotal_64 * numDataTotal_64];
			memcpy(matrixKMod, matrixK, sizeof(double) * (numDataTotal_64 * numDataTotal_64));
			// Solver linear equation with lapack
			OutputFiles::m_logFile << "# Start solve phase for the hat matrix calculation. " << ptrAnalysisControl->outputElapsedTime() << std::endl;
			const long long int nrhsForGCV = numDataTotal_64;
			if (positiveDefinite) {
				ierr = LAPACKE_dpptrs(LAPACK_COL_MAJOR, 'L', numDataTotal_64, nrhsForGCV, matrixToBeInverted, matrixKMod, ldb);
			}
			else {
				ierr = LAPACKE_dsptrs(LAPACK_COL_MAJOR, 'L', numDataTotal_64, nrhsForGCV, matrixToBeInverted, ipiv, matrixKMod, ldb);
			}
			if (ierr < 0) {
				OutputFiles::m_logFile << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
				exit(1);
			}
#ifdef _DEBUG_WRITE
			fout << "matrixK" << std::endl;
			for (int icol = 0; icol < numDataTotal_64; ++icol) {
				for (long long irow = 0; irow < numDataTotal_64; ++irow) {
					fout << matrixK[icol * numDataTotal_64 + irow] << ",";
				}
				fout << std::endl;
			}
			fout << "---";
			fout << "matrixKMod" << std::endl;
			for (int icol = 0; icol < numDataTotal_64; ++icol) {
				for (long long irow = 0; irow < numDataTotal_64; ++irow) {
					fout << matrixKMod[icol * numDataTotal_64 + irow] << ",";
				}
				fout << std::endl;
			}
			fout << "---";
#endif
			// Matrix-Matrix product
			m_traceOfHatMatrixWithoutDampingFactor = 0.0;
			if (CommonParameters::IS_BLAS_USED) {
				double* result = new double[numDataTotal_64 * numDataTotal_64];
				MKL_INT m = static_cast<MKL_INT>(numDataTotal_64);
				MKL_INT n = static_cast<MKL_INT>(numDataTotal_64);
				MKL_INT k = static_cast<MKL_INT>(numDataTotal_64);
				MKL_INT lda = k;
				MKL_INT ldb = k;
				MKL_INT ldc = n;
				double alpha = 1.0;
				double beta = 0.0;
				CBLAS_ORDER  order = CblasColMajor;
				CBLAS_TRANSPOSE transA = CblasNoTrans;
				CBLAS_TRANSPOSE transB = CblasNoTrans;
				cblas_dgemm(order, transA, transB, m, n, k, alpha, matrixK, lda, matrixKMod, ldb, beta, result, ldc);
#ifdef _DEBUG_WRITE
				fout << "hat matrix" << std::endl;
				for (int icol = 0; icol < numDataTotal_64; ++icol) {
					for (long long irow = 0; irow < numDataTotal_64; ++irow) {
						fout << matrixK[icol * numDataTotal_64 + irow] - result[icol * numDataTotal_64 + irow] << ",";
					}
					fout << std::endl;
				}
				fout << "---";
#endif
				for (long long int row = 0; row < numDataTotal_64; ++row) {
					m_traceOfHatMatrixWithoutDampingFactor += matrixK[row * numDataTotal_64 + row] - result[row * numDataTotal_64 + row];
				}
				delete[] result;
			}
			else {
				for (long long iDatLeft = 0; iDatLeft < numDataTotal_64; ++iDatLeft) {
					const long long offsetLeft = iDatLeft * numDataTotal_64;
					double work(0.0);
					for (long long iDatRight = 0; iDatRight < numDataTotal_64; ++iDatRight) {
						// @note: This equation is corredt only if K matrix is symmetric
						work += matrixK[offsetLeft + iDatRight] * matrixKMod[offsetLeft + iDatRight];
					}
					m_traceOfHatMatrixWithoutDampingFactor += matrixK[iDatLeft * numDataTotal_64 + iDatLeft] - work;
				}
#ifdef _DEBUG_WRITE
				fout << "hat matrix" << std::endl;
				for (long long iDatLeft = 0; iDatLeft < numDataTotal_64; ++iDatLeft) {
					for (long long iDatRight = 0; iDatRight < numDataTotal_64; ++iDatRight) {
						double work(0.0);
						for (long long i = 0; i < numDataTotal_64; ++i) {
							work += matrixK[iDatLeft * numDataTotal_64 + i] * matrixKMod[iDatRight * numDataTotal_64 + i];
						}
						fout << matrixK[iDatLeft * numDataTotal_64 + iDatRight] - work << ",";
					}
					fout << std::endl;
				}
				fout << "---";
#endif
			}
			OutputFiles::m_logFile << "# Trace of hat matrix without damping factor : " << m_traceOfHatMatrixWithoutDampingFactor << " ." << ptrAnalysisControl->outputElapsedTime() << std::endl;
			delete[] matrixKMod;
		}
		//-----------------------------

		if( matrixToBeInverted != NULL ){
			delete [] matrixToBeInverted;
			matrixToBeInverted = NULL;
		}
		if( !positiveDefinite ){
			delete [] ipiv;
			ipiv = NULL;
		}
		delete[] matrixK;

#ifdef _DEBUG_WRITE
		for( int iDat = 0; iDat < numDataTotal; ++iDat ){
			fout << "iDat, rhsVectorGlobal : " << iDat << " " << rhsVectorGlobal[iDat] << std::endl;
		}
#endif

	}//----- Treadted by only PE 0 ------------------<<<<<<<<<<<<<<<<<<<<<<<<
	
	//-----------------------------------------------------------------
	// Scatter result vector
	//-----------------------------------------------------------------
	double* dataVectorThisPE = new double [numDataThisPE];
	MPI_Scatterv( rhsVectorGlobal, numDataLocal, numDataAccumulated, MPI_DOUBLE, dataVectorThisPE, numDataThisPE, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if( rhsVectorGlobal != NULL ){
		delete [] rhsVectorGlobal;
		rhsVectorGlobal= NULL;
	}
	if (numDataLocal != NULL) {
		delete[] numDataLocal;
		numDataLocal = NULL;
	}
	if (numDataAccumulated != NULL) {
		delete[] numDataAccumulated;
		numDataAccumulated = NULL;
	}

#ifdef _DEBUG_WRITE
	for( int iDat = 0; iDat < numDataThisPE; ++iDat ){
		fout << "iDat, dataVectorThisPE : " << iDat << " " << dataVectorThisPE[iDat] << std::endl;
	}
#endif

	//----------------------------------------------------------------------------------------------------------------
	// Matrix-vector product
	//----------------------------------------------------------------------------------------------------------------
	// dataVectorThisPE : inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
	// workVector : [J]T*inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
	//----------------------------------------------------------------------------------------------------------------
	double* workVector = new double[numModel];
	for( int iMdl = 0; iMdl < numModel; ++iMdl ){
		workVector[iMdl] = 0.0;// Zero clear
	}
	offset = 0;
	for( int iFreq = 0; iFreq < nFreqThisPE; ++iFreq ){
		assert( offset == ptrObservedData->getNumObservedDataThisPEAccumulated( iFreq ) );
		const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);
		std::ostringstream fileName;
		if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
		}
		fileName << "sensMatFreq" << freqID;
		OutputFiles::m_logFile << "# Read sensitivity matrix from "<< fileName.str() << "." << ptrAnalysisControl->outputElapsedTime() << std::endl;
		int numDataThisFreq(0);
		int numModelTemp(0);
		double* sensitivityMatrix(NULL);
		readSensitivityMatrix( fileName.str(), numDataThisFreq, numModelTemp, sensitivityMatrix );
		if( numDataThisFreq != ptrObservedData->getNumObservedDataThisPE( iFreq ) ){
			OutputFiles::m_logFile << "Error : Number of data written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}
		if( numModel != numModelTemp ){
			OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}	
		//-----------------------------------------------	
		// Matrix-vector product
		//-----------------------------------------------
		if (CommonParameters::IS_BLAS_USED) {
			CBLAS_ORDER order = CblasRowMajor;
			CBLAS_TRANSPOSE trans = CblasTrans;
		    MKL_INT m = static_cast<MKL_INT>(numDataThisFreq);
			MKL_INT n = static_cast<MKL_INT>(numModel);
			double alpha = 1.0;
			double beta = 1.0;
			MKL_INT lda = n;
			MKL_INT incx = 1;
			MKL_INT incy = 1;
			cblas_dgemv(order, trans, m, n, alpha, sensitivityMatrix, lda, &dataVectorThisPE[offset], incx, beta, workVector, incy);
		}else{
			const long long numModel_64 = static_cast<long long>(numModel);
			const long long offset_64 = static_cast<long long>(offset);
			const long long nDataEnd = offset_64 + static_cast<long long>(numDataThisFreq);
			for( long long iDat = offset_64; iDat < nDataEnd; ++iDat ){
				const long long offsetTemp = (iDat - offset_64) * numModel_64;
				const double val = dataVectorThisPE[iDat];
				for( long long iMdl = 0; iMdl < numModel_64; ++iMdl ){
					workVector[iMdl] += sensitivityMatrix[ iMdl + offsetTemp ] * val;
				}
			}
		}
		//-----------------------------------------------
		delete [] sensitivityMatrix;
		offset += numDataThisFreq;// Add data number
	}

	if( dataVectorThisPE != NULL ){
		delete [] dataVectorThisPE;
		dataVectorThisPE = NULL;
	}

#ifdef _DEBUG_WRITE
	for( int iMdl = 0; iMdl < numModel; ++iMdl ){
		fout << "iMdl, workVector : " << iMdl << " " << workVector[iMdl] << std::endl;
	}
#endif

	double* workVectorSum(NULL);
	if( myProcessID == 0 ){
		workVectorSum = new double [numModel];
	}
	MPI_Reduce( workVector, workVectorSum, numModel, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	if( workVector != NULL ){
		delete [] workVector ;
		workVector = NULL;
	}

#ifdef _DEBUG_WRITE
	if( myProcessID == 0 ){
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			fout << "iMdl, workVectorSum : " << iMdl << " " << workVectorSum[iMdl] << std::endl;
		}
	}
#endif

	double* resultVector = new double[numModel];
	if( myProcessID == 0 ){
		//-------------------------------------------------------------------------------------------------------------------
		// workVectorSum : [J]T*inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
		// dVector : -RTRm-[deriv_G]T*G+[J]T*rd
		//   => [-RTRm-[deriv_G]T*G+[J]T*rd] - [J]T*inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
		//-------------------------------------------------------------------------------------------------------------------
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			dVector[iMdl] -= workVectorSum[iMdl];
		}
		delete [] workVectorSum;
		//-------------------------------------------------------------------------------------------------------------------

#ifdef _DEBUG_WRITE
		for( int iMdl = 0; iMdl < numModel; ++iMdl ){
			fout << "iMdl, dVector : " << iMdl << " " << dVector[iMdl] << std::endl;
		}
#endif
		//-------------------------------------------------------------------------------------------------------------------
		// dVector : [-RTRm-[deriv_G]T*G+[J]T*rd] - [J]T*inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
		// resultVector : inv([R]T*[R]) * [-RTRm-[deriv_G]T*G+[J]T*rd] - [J]T*inv([I]+[J]*inv([R]T*[R])*[J]T)*[J]*inv([R]T*[R])*[-RTRm-[deriv_G]T*G+[J]T*rd]
		//-------------------------------------------------------------------------------------------------------------------
		RTRMatrix.solvePhaseMatrixSolver( 1, dVector, resultVector );
		RTRMatrix.releaseMemory();
		delete [] dVector;
		//-------------------------------------------------------------------------------------------------------------------
	}

	MPI_Bcast( resultVector, numModel, MPI_DOUBLE, 0, MPI_COMM_WORLD );

#ifdef _DEBUG_WRITE
	for( int iMdl = 0; iMdl < numModel; ++iMdl ){
		fout << "iMdl, resultVector : " << iMdl << " " << resultVector[iMdl] << std::endl;
	}
#endif

	//=================================================================
	// Calculate predicted vector of weighted observed data for each PE
	//=================================================================
	if (ptrAnalysisControl->needsGCVCalculation()) {
		OutputFiles::m_logFile << "# Calculate predicted vector of weighted observed data." << ptrAnalysisControl->outputElapsedTime() << std::endl;
		if (m_predictedObservedDataVectorThisPE == NULL){
			m_predictedObservedDataVectorThisPE = new double[numDataThisPE];
		}
		calculatePredictedVectorOfWeightedObservedDataLocal(numDataThisPE, nFreqThisPE, numModel, CommonParameters::IS_BLAS_USED, resultVector, m_predictedObservedDataVectorThisPE);
	}

	//=================================================================
	// Calculate model increments
	//=================================================================
	OutputFiles::m_logFile << "# Calculate model increments. " << ptrAnalysisControl->outputElapsedTime() << std::endl;

	//=================================================================
	// Update resistivity values
	//=================================================================
	ptrResistivityBlock->calcUnfixedAnisotropicResistivityParametersUpdatedFull( resultVector );
	ptrResistivityBlock->updateResistivityValues();

	//=================================================================
	// Update distortion parameters
	//=================================================================
	ptrObservedData->calcDistortionParamsUpdatedFullFromIncrements( &resultVector[nBlkNotFixed] );
	ptrObservedData->updateDistortionParams();

	delete [] resultVector;

	//=================================================================
	// Synchronize
	//=================================================================
	MPI_Barrier( MPI_COMM_WORLD );

}

// Calculate constraining matrices
void InversionGaussNewtonDataSpaceAnisotropic::calcConstrainingMatrices(DoubleSparseMatrix& constrainingMatrixFull, DoubleSparseMatrix& constrainingMatrixWithoutCrossProductTerm,
	DoubleSparseMatrix& derivativeMatrixOfCrossProductTerm, std::vector<double>& crossProductVector) const{

	const int numDistortionParamsNotFixed = (ObservedData::getInstance())->getNumDistortionParamsNotFixed();
	const ResistivityBlockAnisotropic* const ptrResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	const int numberOfUnfixedResistivityParameters = ptrResistivityBlock->getNumberOfUnfixedResistivityParameters();
	const int numModel = numberOfUnfixedResistivityParameters + numDistortionParamsNotFixed;

	ptrResistivityBlock->confirmNoElementHavingGeneralAnisotropy();
	
	// Difference filters
	std::vector< std::vector<int> > nonZeroColsWithoutCrossProduct;
	std::vector< std::vector<double> > matValuesWithoutCrossProduct;
	std::vector<double> rhsValuesWithoutCrossProduct;
	const int nElem = (AnalysisControl::getInstance())->getPointerOfMeshData()->getNumElemTotal();
	nonZeroColsWithoutCrossProduct.reserve(nElem);
	matValuesWithoutCrossProduct.reserve(nElem);
	rhsValuesWithoutCrossProduct.reserve(nElem);
	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const double factor1 = ptrAnalysisControl->getTradeOffParameterForResistivityValue();
	ptrResistivityBlock->calcSpatialDifferenceFilterForTransverselyIsotropicConductivity(factor1, nonZeroColsWithoutCrossProduct, matValuesWithoutCrossProduct, rhsValuesWithoutCrossProduct);
	const double factor2 = ptrAnalysisControl->getTradeOffParameterForDegreeOfAnisotropy();
	ptrResistivityBlock->calcDifferenceFilterBetweenComponentsForTransverselyIsotropicConductivity(factor2, nonZeroColsWithoutCrossProduct, matValuesWithoutCrossProduct, rhsValuesWithoutCrossProduct);
	addGalvanicDistortionTermToConstrainingMatrix(numberOfUnfixedResistivityParameters, nonZeroColsWithoutCrossProduct, matValuesWithoutCrossProduct, rhsValuesWithoutCrossProduct);
	addValuesToConstrainingMatrix(numModel, nonZeroColsWithoutCrossProduct, matValuesWithoutCrossProduct, rhsValuesWithoutCrossProduct, constrainingMatrixWithoutCrossProductTerm);

	// Cross-gradient term
	const double factor3 = ptrAnalysisControl->getTradeOffParameterForAnisotropyDirection();
	ptrResistivityBlock->calcCrossProductForTransverselyIsotropicConductivity(factor3, crossProductVector);
	std::vector< std::vector<int> > nonZeroColsCrossProduct;
	std::vector< std::vector<double> > matValuesCrossProduct;
	std::vector<double> rhsValuesCrossProduct;
	nonZeroColsCrossProduct.reserve(nElem);
	matValuesCrossProduct.reserve(nElem);
	rhsValuesCrossProduct.reserve(nElem);
	ptrResistivityBlock->calcDerivativeOfCrossProductForTransverselyIsotropicConductivity(factor3, nonZeroColsCrossProduct, matValuesCrossProduct, rhsValuesCrossProduct);
	addValuesToConstrainingMatrix(numModel, nonZeroColsCrossProduct, matValuesCrossProduct, rhsValuesCrossProduct, derivativeMatrixOfCrossProductTerm);
#ifdef _DEBUG_WRITE
	int count(0);
	for (std::vector<double>::const_iterator itr = crossProductVector.begin(); itr != crossProductVector.end(); ++itr, ++count){
		std::cout << "i " << count << " " << *itr << std::endl;
	}
#endif

	// Total constraining matrix
	std::vector< std::vector<int> > nonZeroColsFull;
	std::vector< std::vector<double> > matValuesFull;
	std::vector<double> rhsValuesFull;
	const int numRows = static_cast<int>(nonZeroColsWithoutCrossProduct.size()) + static_cast<int>(nonZeroColsCrossProduct.size());
	nonZeroColsFull.reserve(numRows);
	matValuesFull.reserve(numRows);
	rhsValuesFull.reserve(numRows);
	nonZeroColsFull = nonZeroColsWithoutCrossProduct;
	matValuesFull = matValuesWithoutCrossProduct;
	rhsValuesFull = rhsValuesWithoutCrossProduct;
	nonZeroColsFull.insert(nonZeroColsFull.end(), nonZeroColsCrossProduct.begin(), nonZeroColsCrossProduct.end());
	matValuesFull.insert(matValuesFull.end(), matValuesCrossProduct.begin(), matValuesCrossProduct.end());
	rhsValuesFull.insert(rhsValuesFull.end(), rhsValuesCrossProduct.begin(), rhsValuesCrossProduct.end());
	addValuesToConstrainingMatrix(numModel, nonZeroColsFull, matValuesFull, rhsValuesFull, constrainingMatrixFull);

}

// Add the galvanic distortion term to constraining matrix
void InversionGaussNewtonDataSpaceAnisotropic::addGalvanicDistortionTermToConstrainingMatrix(const int numberOfUnfixedResistivityParameters, std::vector< std::vector<int> >& nonZeroCols,
	std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const {

	const ObservedData* const ptrObservedData = ObservedData::getInstance();
	const int numDistortionParamsNotFixed = ptrObservedData->getNumDistortionParamsNotFixed();

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	int iMdl = numberOfUnfixedResistivityParameters;
	if (ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		// For distortion matrix components
		const double factor = ptrAnalysisControl->getTradeOffParameterForDistortionMatrixComplexity();
		for (int iParamsNotFixed = 0; iParamsNotFixed < numDistortionParamsNotFixed; ++iParamsNotFixed) {
			std::vector<int> cols;
			cols.push_back(iMdl);
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			vals.push_back(factor);
			matValues.push_back(vals);
			rhsValues.push_back(0.0);
			++iMdl;
		}
	}
	else if (ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
		// For gains and rotations of distortion matrix
		const double factor2 = ptrAnalysisControl->getTradeOffParameterForGainsOfDistortionMatrix();
		const double factor3 = ptrAnalysisControl->getTradeOffParameterForRotationsOfDistortionMatrix();
		for (int iParamsNotFixed = 0; iParamsNotFixed < numDistortionParamsNotFixed; ++iParamsNotFixed) {
			std::vector<int> cols;
			cols.push_back(iMdl);
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			if (ptrObservedData->getTypesOfDistortionParamsNotFixed(iParamsNotFixed) == ObservedDataStationMT::EX_GAIN ||
				ptrObservedData->getTypesOfDistortionParamsNotFixed(iParamsNotFixed) == ObservedDataStationMT::EY_GAIN) {
				// Gains
				vals.push_back(factor2);
			}
			else {
				// Rotations
				vals.push_back(factor3);
			}
			matValues.push_back(vals);
			rhsValues.push_back(0.0);
			++iMdl;
		}
	}
	else if (ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		// For gains of distortion matrix
		const double factor2 = ptrAnalysisControl->getTradeOffParameterForGainsOfDistortionMatrix();
		for (int iParamsNotFixed = 0; iParamsNotFixed < numDistortionParamsNotFixed; ++iParamsNotFixed) {
			std::vector<int> cols;
			cols.push_back(iMdl);
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			vals.push_back(factor2);
			matValues.push_back(vals);
			rhsValues.push_back(0.0);
			++iMdl;
		}
	}

}

// Add values to constraing matrix
void InversionGaussNewtonDataSpaceAnisotropic::addValuesToConstrainingMatrix(const int numModel, const std::vector< std::vector<int> >& nonZeroCols,
	const std::vector< std::vector<double> >& matValues, const std::vector<double>& rhsValues, DoubleSparseMatrix& constrainingMatrix) const{

	assert(nonZeroCols.size() == matValues.size());
	assert(nonZeroCols.size() == rhsValues.size());
	const int numRows = static_cast<int>(nonZeroCols.size());
	constrainingMatrix.setNumRowsAndColumns(numRows, numModel);

	for (int iRow = 0; iRow < numRows; ++iRow) {
		int count(0);
		for (std::vector<int>::const_iterator itr = nonZeroCols[iRow].begin(); itr != nonZeroCols[iRow].end(); ++itr, ++count) {
			constrainingMatrix.setStructureAndAddValueByTripletFormat(iRow, *itr, matValues[iRow][count]);
		}
		constrainingMatrix.addRightHandSideVector(iRow, rhsValues[iRow]);
	}
	constrainingMatrix.convertToCRSFormat();

}

// Calculate predicted vector of weighted observed data for each PE
void InversionGaussNewtonDataSpaceAnisotropic::calculatePredictedVectorOfWeightedObservedDataLocal(	const int numDataThisPE, const int nFreqThisPE, 
	const int numModel, const bool useBLAS, const double* const resultVector, double* predictedVectorLocal) const {

	int offset(0);
	for (int i = 0; i < numDataThisPE; ++i) {
		predictedVectorLocal[i] = 0.0;// Zero clear
	}
	const ObservedData* const ptrObservedData = ObservedData::getInstance();
	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	for (int iFreq = 0; iFreq < nFreqThisPE; ++iFreq) {
		const int freqID = ptrObservedData->getIDsOfFrequenciesCalculatedByThisPE(iFreq);
		std::ostringstream fileName;
		if (!ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
#ifdef _LINUX
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\/";
#else
			fileName << ptrAnalysisControl->getDirectoryOfOutOfCoreFilesForSensitivityMatrix() + "\\";
#endif
		}
		fileName << "sensMatFreq" << freqID;
		OutputFiles::m_logFile << "# Read sensitivity matrix from " << fileName.str() << "." << ptrAnalysisControl->outputElapsedTime() << std::endl;
		int numDataThisFreq(0);
		int numModelTemp(0);
		double* sensitivityMatrix(NULL);
		readSensitivityMatrix(fileName.str(), numDataThisFreq, numModelTemp, sensitivityMatrix);
		if (numDataThisFreq != ptrObservedData->getNumObservedDataThisPE(iFreq)) {
			OutputFiles::m_logFile << "Error : Number of data written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}
		if (numModel != numModelTemp) {
			OutputFiles::m_logFile << "Error : Number of model written in out-of-core file is unequal to the internal one !!" << std::endl;
			exit(1);
		}
		if (CommonParameters::IS_BLAS_USED) {
			CBLAS_ORDER order = CblasRowMajor;
			CBLAS_TRANSPOSE trans = CblasNoTrans;
			MKL_INT m = static_cast<MKL_INT>(numDataThisFreq);
			MKL_INT n = static_cast<MKL_INT>(numModel);
			double alpha = 1.0;
			double beta = 0.0;
			MKL_INT lda = n;
			MKL_INT incx = 1;
			MKL_INT incy = 1;
			cblas_dgemv(order, trans, m, n, alpha, sensitivityMatrix, lda, resultVector, incx, beta, &predictedVectorLocal[offset], incy);
		}
		else {
			long long iDat(0);
			long long offsetTemp(0);
			double work(0.0);
			long long iMdl(0);
			const long long numDataThisFreq_64 = static_cast<long long>(numDataThisFreq);
			const long long numModel_64 = static_cast<long long>(numModel);
			const long long offset_64 = static_cast<long long>(offset);
#ifdef _USE_OMP
#pragma omp parallel for default(shared) private( iDat, offsetTemp, work, iMdl )
#endif
			for (iDat = 0; iDat < numDataThisFreq_64; ++iDat) {
				offsetTemp = iDat * numModel_64;
				work = 0.0;
				for (iMdl = 0; iMdl < numModel_64; ++iMdl) {
					work += sensitivityMatrix[iMdl + offsetTemp] * resultVector[iMdl];
				}
				predictedVectorLocal[offset_64 + iDat] = work;
			}
		}
		delete[] sensitivityMatrix;
		offset += numDataThisFreq;// Add data number
	}
#ifdef _DEBUG_WRITE
	for (int i = 0; i < numDataThisPE; ++i) {
		std::cout << "i predictedVectorGlobal[i] bfr" << i << " " << predictedVectorLocal[i] << std::endl;
	}
#endif

}
