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
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>

#include "AnalysisControl.h"
#include "ResistivityBlock.h"
#include "MeshData.h"
#include "MeshDataBrickElement.h"
#include "MeshDataNonConformingHexaElement.h"
#include "mpi.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "ObservedData.h"
#include "OutputFiles.h"
#include "InversionGaussNewtonModelSpace.h"
#include "InversionGaussNewtonDataSpace.h"
#include "Forward3DBrickElement0thOrder.h"
#include "Forward3DTetraElement0thOrder.h"
#include "mkl.h"
#include <assert.h>

#ifdef _USE_OMP
#include <omp.h>
#endif

// Return the the instance of the class
AnalysisControl* AnalysisControl::getInstance(){
   	static AnalysisControl instance;// The only instance
  	return &instance;
}

// Constructer
AnalysisControl::AnalysisControl():
	m_myPE(NULL),
	m_totalPE(NULL),
	m_numThreads(1),
	m_startTime(NULL),
	m_boundaryConditionBottom(AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR),
	m_orderOfFiniteElement(0),
	m_modeOfPARDISO(PARDISOSolver::INCORE_MODE),
	m_numberingMethod(AnalysisControl::NOT_ASSIGNED),
	m_isOutput2DResult(false),
	m_tradeOffParameterForResistivityValue(0.0),
	m_tradeOffParameterForDistortionMatrixComplexity(0.0),
	m_tradeOffParameterForDistortionGain(0.0),
	m_tradeOffParameterForDistortionRotation(0.0),
	m_iterationNumInit(0),
	m_iterationNumMax(0),
	m_thresholdValueForDecreasing(0.001),
	m_decreaseRatioForConvegence(1.0),
	m_stepLengthDampingFactorCur(0.5),
	m_stepLengthDampingFactorMin(0.1),
	m_stepLengthDampingFactorMax(1.0),
	m_numOfIterIncreaseStepLength(3),
	m_factorDecreasingStepLength(0.50),
	m_factorIncreasingStepLength(1.25),
	m_numCutbackMax(5),
	m_holdMemoryForwardSolver(false),
	m_ptrForward3DBrickElement0thOrder(NULL),
	m_ptrForward3DTetraElement0thOrder(NULL),
	m_ptrInversion(NULL),
	m_objectFunctionalPre(0.0),
	m_dataMisfitPre(0.0),
	m_modelRoughnessPre(0.0),
	m_normOfDistortionMatrixDifferencesPre(0.0),
	m_normOfGainsPre(0.0),
	m_normOfRotationsPre(0.0),
	m_numConsecutiveIterFunctionalDecreasing(0),
	m_continueWithoutCutback(false),
	m_maxMemoryPARDISO(3000),
	m_typeOfMesh(MeshData::HEXA),
	m_typeOfRoughningMatrix(AnalysisControl::USE_ELEMENTS_SHARE_FACES),
	m_typeOfElectricField(AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD),
	m_isTypeOfElectricFieldSetIndivisually(false),
	m_typeOfOwnerElement(AnalysisControl::USE_LOWER_ELEMENT),
	m_isTypeOfOwnerElementSetIndivisually(false),
	m_divisionNumberOfMultipleRHSInForward(1),
	m_divisionNumberOfMultipleRHSInInversion(1),
	m_binaryOutput(true),
	m_positiveDefiniteNormalEqMatrix(false),
	m_typeOfDistortion(AnalysisControl::DISTORTION_TYPE_UNDEFINED),
	m_inversionMethod(Inversion::GAUSS_NEWTON_MODEL_SPECE),
	m_isObsLocMovedToCenter(false),
	m_apparentResistivityAndPhaseTreatmentOption(NO_SPECIAL_TREATMENT_APP_AND_PHASE),
	m_isRougheningMatrixOutputted(false),
#ifdef _ANISOTOROPY
	m_typeOfDataSpaceAlgorithm(NEW_DATA_SPACE_ALGORITHM),
	m_typeOfAnisotropy(NO_ANISOTROPY)
#else
	m_typeOfDataSpaceAlgorithm(NEW_DATA_SPACE_ALGORITHM),
	m_useDifferenceFilter(false),
	m_degreeOfLpOptimization(2),
	m_lowerLimitOfDifflog10RhoForLpOptimization(0.01),
	m_upperLimitOfDifflog10RhoForLpOptimization(2.0),
	m_maxIterationIRWLSForLpOptimization(3),
	m_thresholdIRWLSForLpOptimization(1.0)
#endif
{
	for( int iDir = 0; iDir < 3; ++iDir ){
		m_alphaWeight[iDir] = 1.0;
	}

#ifdef _USE_OMP
	m_numThreads = omp_get_num_threads();
#endif

	// Measure the start time
	time(&m_startTime);

	// Get process ID and total process number
	MPI_Comm_rank ( MPI_COMM_WORLD, &m_myPE );
	MPI_Comm_size ( MPI_COMM_WORLD, &m_totalPE );
	
	// Assign flag specifing use backward element to the member variable specifing
	// which backward or forward element is used for calculating EM field
	m_useBackwardOrForwardElement.directionX = AnalysisControl::BACKWARD_ELEMENT;
	m_useBackwardOrForwardElement.directionY = AnalysisControl::BACKWARD_ELEMENT;

}

//Destructer
AnalysisControl::~AnalysisControl(){

	if( !m_outputParametersForVis.empty() ){
		m_outputParametersForVis.clear();
	}

	if( m_ptrForward3DBrickElement0thOrder != NULL ){
		delete m_ptrForward3DBrickElement0thOrder;
		m_ptrForward3DBrickElement0thOrder = NULL;
	}

	if( m_ptrForward3DTetraElement0thOrder != NULL ){
		delete m_ptrForward3DTetraElement0thOrder;
		m_ptrForward3DTetraElement0thOrder = NULL;
	}

	if( m_ptrInversion != NULL ){
		delete m_ptrInversion;
		m_ptrInversion = NULL;
	}

}

void AnalysisControl::run(){

#ifdef _LINUX
	OutputFiles::m_logFile << "# Start " << CommonParameters::programName << " Linux Version " << CommonParameters::versionID << " " << outputElapsedTime() << std::endl;
#else
	OutputFiles::m_logFile << "# Start " << CommonParameters::programName << " Windows Version " << CommonParameters::versionID << " " << outputElapsedTime() << std::endl;
#endif

	// Get process ID 
	const int myProcessID = getMyPE();

	//---------------------------------------------------
	//--- Read analysis control data from control.dat ---
	//---------------------------------------------------
	OutputFiles::m_logFile << "# Read analysis control data from control.dat." << outputElapsedTime() << std::endl;
	inputControlData();

	//-------------------------------------------------------
	//--- Create object of Forward analysis and inversion ---
	//-------------------------------------------------------
	if( m_typeOfMesh == MeshData::HEXA ){
		m_ptrForward3DBrickElement0thOrder = new Forward3DBrickElement0thOrder();
	}else if( m_typeOfMesh == MeshData::TETRA ){
		m_ptrForward3DTetraElement0thOrder = new Forward3DTetraElement0thOrder();
	}else if( m_typeOfMesh == MeshData::NONCONFORMING_HEXA ){
		m_ptrForward3DNonConformingHexaElement0thOrder = new Forward3DNonConformingHexaElement0thOrder();
	}else{
		OutputFiles::m_logFile << "Error : Type of mesh is wrong !! : " << m_typeOfMesh << "." << std::endl;
		exit(1);
	}

	switch(getInversionMethod()){
		case Inversion::GAUSS_NEWTON_MODEL_SPECE:
			m_ptrInversion = new InversionGaussNewtonModelSpace();
			break;
		case Inversion::GAUSS_NEWTON_DATA_SPECE:
			m_ptrInversion = new InversionGaussNewtonDataSpace();
			break;
		default:
			OutputFiles::m_logFile << "Error : Type of inversion method is wrong  !! : " << getInversionMethod() << std::endl;
			exit(1);
			break;
	}

	//------------------------------------
	//--- Read mesh data from mesh.dat ---
	//------------------------------------
	OutputFiles::m_logFile << "# Read mesh data from mesh.dat ." << outputElapsedTime() << std::endl;
	getPointerOfForward3D()->callInputMeshData();

	//---------------------------------------------
	//--- Read data of resisitivity block model ---
	//---------------------------------------------
	OutputFiles::m_logFile << "# Read data of resisitivity block model ." << outputElapsedTime() << std::endl;
	ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();	
	pResistivityBlock->inputResisitivityBlock();

	//-------------------------------------------
	//--- Read observed data from observe.dat ---
	//-------------------------------------------
	OutputFiles::m_logFile << "# Read observed data ." << outputElapsedTime() << std::endl;
	ObservedData* pObservedData = ObservedData::getInstance();	
	pObservedData->inputObservedData();
	pObservedData->calcFrequenciesCalculatedByThisPE();

	const int nfreq = pObservedData->getTotalNumberOfDifferenetFrequencies();
	if( nfreq <= 0 ){
		OutputFiles::m_logFile << "Error : Total number of frequencies is less than zero !!" << std::endl;
		exit(1);
	}

	//-----------------------------------
	//--- Read distortion matrix data ---
	//-----------------------------------
	if( estimateDistortionMatrix() ){
		OutputFiles::m_logFile << "# Read distortion matrix data ." << outputElapsedTime() << std::endl;
		//pObservedData->inputStaticShiftData();
		pObservedData->inputDistortionMatrixData();
	}

	//------------------------------------------------------------------------------------------
	//--- Allocate memory for the calculated values and errors of all stations               ---
	//--- after setting up frequencies calculated by this PE, at which observed value exists ---
	//------------------------------------------------------------------------------------------
	pObservedData->allocateMemoryForCalculatedValuesOfAllStations();

	//-------------------------------------------
	//--- Find element including each station ---
	//-------------------------------------------
	OutputFiles::m_logFile << "# Find element including each station ." << outputElapsedTime() << std::endl;
	pObservedData->findElementIncludingEachStation();

	//------------------------------------------------------------------------
	//--- Output information of locations of observed stations to vtk file ---
	//------------------------------------------------------------------------
	OutputFiles* const ptrOutputFiles = OutputFiles::getInstance();
	if( myProcessID == 0 ){// If this PE number is zero
		ptrOutputFiles->openVTKFileForObservedStation();
		pObservedData->outputLocationsOfObservedStationsToVtk();
	}

	//----------------------------------------------------------------
	//--- Initialize response functions and errors of all stations ---
	//----------------------------------------------------------------
	pObservedData->initializeResponseFunctionsAndErrorsOfAllStations();

	//-----------------------------------------------------
	//--- Output number of model parameters to log file ---
	//-----------------------------------------------------
	m_ptrInversion->outputNumberOfModel();

	//-----------------------------------
	//--- Calculate Roughening Matrix ---
	//-----------------------------------
	OutputFiles::m_logFile << "# Calculate Roughening Matrix ." << outputElapsedTime() << std::endl;
	pResistivityBlock->calcRougheningMatrix();

	//-------------------------------------------
	//--- Output geometory file and case file ---
	//-------------------------------------------
	if( !m_outputParametersForVis.empty() && writeBinaryFormat() && myProcessID == 0 ){// Write to BINARY file
		ptrOutputFiles->outputCaseFile();
		getPointerOfMeshData()->outputMeshDataToBinary();
		pResistivityBlock->outputResistivityDataToBinary();
	}

	//---------------------
	//--- Open cnv file ---
	//---------------------
	ptrOutputFiles->openCnvFile();

	AnalysisControl::ConvergenceBehaviors convergenceFlag = AnalysisControl::DURING_RETRIALS;

	for( int iter = m_iterationNumInit; iter <= m_iterationNumMax; ++iter ){

		m_iterationNumCurrent = iter;

		// Open csv file in which the results of 2D forward computations is written
		if( m_isOutput2DResult ){
			ptrOutputFiles->openCsvFileFor2DFwd( iter );
		}

		if( !m_outputParametersForVis.empty() && !writeBinaryFormat() ){
			ptrOutputFiles->openVTKFile( iter );
			pResistivityBlock->outputResistivityDataToVTK();
		}

		if( m_iterationNumMax > iter && doesOutputToVTK( AnalysisControl::OUTPUT_SENSITIVITY) ){// if output sensitivity
			m_ptrInversion->allocateMemoryForSensitivityScalarValues();
		}

		int iCutBack = 0;
		for( ; iCutBack <= m_numCutbackMax; ++iCutBack ){

			OutputFiles::m_logFile << "###############################################################################" << std::endl;
			OutputFiles::m_logFile << "# Start Forward Computation.  Iteration : " << iter << ",  Retrial : " << iCutBack << std::endl;
			OutputFiles::m_logFile << "###############################################################################" << std::endl;

			//---------------------------
			//--- Forward computation ---
			//---------------------------
			calcForwardComputation( iter );

			//--------------------------------------------
			//--- Output information about convergence ---
			//--------------------------------------------
			convergenceFlag = adjustStepLengthDampingFactor( iter, iCutBack );
			if( convergenceFlag == AnalysisControl::GO_TO_NEXT_ITERATION || convergenceFlag == AnalysisControl::INVERSIN_CONVERGED ){
				break;// Go out of the loop
			}

			//-----------------------------------------------------------
			//--- Change resistivity values and distortion parameters ---
			//-----------------------------------------------------------
			pResistivityBlock->updateResistivityValues();
			pObservedData->updateDistortionParams();
		}

		if( iCutBack > m_numCutbackMax ){
			OutputFiles::m_logFile << "# Reach maximum retrial number." << std::endl;
			break;
		}

		// Output induction arrows to vtk file
		pObservedData->outputInductionArrowToVtk( iter );

		// Open csv file in which the results of 3D forward computations is written
		ptrOutputFiles->openCsvFileFor3DFwd( iter );
		// Output results
		pObservedData->outputCalculatedValuesOfAllStations();

		// Output resistivity model
		if( writeBinaryFormat() ){// Write to BINARY file
			if( myProcessID == 0 ){
				pResistivityBlock->outputResistivityValuesToBinary( iter );
			}
		}
		else{// Write to ASCII file
			pResistivityBlock->outputResistivityValuesToVTK();
		}

		if( myProcessID == 0 && iter > m_iterationNumInit ){// If this PE number is zero and iteration number is not the first one
			pResistivityBlock->outputResisitivityBlock( iter );
			if( estimateDistortionMatrix() ){
				pObservedData->outputDistortionParams( iter );
			}
		}

		// Output sensitivity
		if( doesCalculateSensitivity(iter) && doesOutputToVTK( AnalysisControl::OUTPUT_SENSITIVITY) ){// if output sensitivity			
			if( writeBinaryFormat() ){// Write to BINARY file
				m_ptrInversion->outputSensitivityScalarValuesToBinary( iter );
			}
			else{// Write to ASCII file
				m_ptrInversion->outputSensitivityScalarValuesToVtk();
			}
			m_ptrInversion->releaseMemoryOfSensitivityScalarValues();
		}

		//-----------------
		//--- Inversion ---
		//-----------------
		if( convergenceFlag == AnalysisControl::INVERSIN_CONVERGED ){
			OutputFiles::m_logFile << "# Converged." << std::endl;
			break;
		}

		if( iter >= m_iterationNumMax ){
			OutputFiles::m_logFile << "# Reach maximum iteration number." << std::endl;
			break;
		}

		OutputFiles::m_logFile << "###############################################################################" << std::endl;
		OutputFiles::m_logFile << "# Start Inversion.  Iteration : " << iter << std::endl;
		OutputFiles::m_logFile << "###############################################################################" << std::endl;

		ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();
		ObservedData* const ptrObservedData = ObservedData::getInstance();
		ptrResistivityBlock->copyResistivityValuesNotFixedCurToPre();
		ptrObservedData->copyDistortionParamsCurToPre();
		if( useDifferenceFilter() ){
			const int maxIter = getMaxIterationIRWLSForLpOptimization();
			double modelRoughnessPre = ptrResistivityBlock->calcModelRoughnessForDifferenceFilter();
			for( int iter = 0; iter < maxIter; ++iter ){
				OutputFiles::m_logFile << "# Iteration number of reweighted iterative algorithm for Lp optimization : " << iter+1 << std::endl;
				m_ptrInversion->inversionCalculation();
				const double modelRoughness = ptrResistivityBlock->calcModelRoughnessForDifferenceFilter();
				OutputFiles::m_logFile << "# Model-roughness is changed from " << modelRoughnessPre << " to " << modelRoughness << std::endl;
				if( modelRoughnessPre > CommonParameters::EPS &&
					fabs( modelRoughness - modelRoughnessPre ) / modelRoughnessPre < getThresholdIRWLSForLpOptimization() * 0.01 ){
					OutputFiles::m_logFile << "# Reweighted iterative algorithm for Lp optimization is converged" << std::endl;
					break;
				}
				modelRoughnessPre = modelRoughness;
			}
			// Delete old out-of-core files
			m_ptrInversion->deleteOutOfCoreFileAll();
		}else{
			m_ptrInversion->inversionCalculation();
		}
	}

	if( m_ptrForward3DBrickElement0thOrder != NULL ){
		delete m_ptrForward3DBrickElement0thOrder;
		m_ptrForward3DBrickElement0thOrder = NULL;
	}

	if( m_ptrForward3DTetraElement0thOrder != NULL ){
		delete m_ptrForward3DTetraElement0thOrder;
		m_ptrForward3DTetraElement0thOrder = NULL;
	}

	if( m_ptrInversion != NULL ){
		delete m_ptrInversion;
		m_ptrInversion = NULL;
	}

	m_ptrInversion->deleteOutOfCoreFileAll();

	OutputFiles::m_logFile << "# End " << CommonParameters::programName << " " << outputElapsedTime() << std::endl;

}

// Read analysis control data from "control.dat"
void AnalysisControl::inputControlData(){

	std::ifstream inFile( "control.dat", std::ios::in );
	if( inFile.fail() )
	{
		//std::cerr << "File open error : control.dat !!" << std::endl;
		OutputFiles::m_logFile << "File open error : control.dat !!" << std::endl;
		exit(1);
	}

	// Flag specifing whether each parameter has already read from control.dat
	bool hasAlreadyRead[numParamWrittenInControlFile];
	for( int i = 0; i < numParamWrittenInControlFile; ++i ){
		hasAlreadyRead[i] = false;
	}


	ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();

	while(!inFile.eof())
	{
		std::string line;
		inFile >> line;

#ifdef _DEBUG_WRITE
		std::cout << "line : " << line << std::endl;
#endif

		double dbuf(0.0);
		int ibuf(0);
		if( line.substr(0,25).compare("BOUNDARY_CONDITION_BOTTOM") == 0 ){// Read the type of boundary condition at the bottom of the model
			const int paramID= AnalysisControl::BOUNDARY_CONDITION_BOTTOM;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : BOUNDARY_CONDITION_BOTTOM" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			if( ibuf != AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL &&
				ibuf != AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
				OutputFiles::m_logFile << "Error : Wrong type of boundary condition at the bottom of the model !! " << ibuf << "." << std::endl;
				exit(1);
			}
			m_boundaryConditionBottom = ibuf;
			hasAlreadyRead[paramID] = true;
		//}else if( line.substr(0,5).compare("ORDER") == 0 ){// Read the order of finite element
		//	const int paramID= AnalysisControl::ORDER;
		//	if( hasAlreadyRead[paramID] == true ){
		//		OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ORDER" << std::endl;
		//		exit(1);
		//	}
		//	inFile >> ibuf;
		//	if( ibuf != 0 && ibuf != 1 ){
		//		OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! " << ibuf << "." << std::endl;
		//		exit(1);
		//	}
		//	m_orderOfFiniteElement = ibuf;
		//	hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,9).compare("MESH_TYPE") == 0 ){// Type of mesh
			const int paramID= AnalysisControl::MESH_TYPE;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : MESH_TYPE" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			if( ibuf != MeshData::HEXA && ibuf != MeshData::TETRA && ibuf != MeshData::NONCONFORMING_HEXA ){
				OutputFiles::m_logFile << "Error : The number following MESH_TYPE must be " << MeshData::HEXA << ", " << MeshData::TETRA << " or " << MeshData::NONCONFORMING_HEXA << " !!" << std::endl;
				exit(1);
			}
			m_typeOfMesh = ibuf;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,11).compare("NUM_THREADS") == 0 ){// Read total number of threads
			const int paramID= AnalysisControl::NUM_THREADS;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : NUM_THREADS" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			if( ibuf < 0 ){
				OutputFiles::m_logFile << "Error : Number of threads must be greater than or equals to 1 !! " << std::endl;
				exit(1);
			}
			m_numThreads = ibuf;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,10).compare("FWD_SOLVER") == 0 ){
			const int paramID= AnalysisControl::FWD_SOLVER;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : FWD_SOLVER" << std::endl;
				exit(1);
			}			
			inFile >> ibuf;
			if( ibuf != PARDISOSolver::INCORE_MODE && ibuf != PARDISOSolver::SELECT_MODE_AUTOMATICALLY && ibuf !=PARDISOSolver::OUT_OF_CORE_MODE ){
				OutputFiles::m_logFile << "Error : Parameter specifing the mode of forward solver must be 0, 1 or 2 !! " << ibuf << "." << std::endl;
				exit(1);
			}
			m_modeOfPARDISO = ibuf;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,9).compare("MEM_LIMIT") == 0 ){
			const int paramID= AnalysisControl::MEM_LIMIT;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : MEM_LIMIT" << std::endl;
				exit(1);
			}
			inFile >> dbuf;
			m_maxMemoryPARDISO = static_cast<int>(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,16).compare("NUMBERING_METHOD") == 0 ){
			const int paramID= AnalysisControl::NUMBERING_METHOD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : NUMBERING_METHOD" << std::endl;
				exit(1);
			}			
			inFile >> ibuf;
			if( ibuf != AnalysisControl::NOT_ASSIGNED && ibuf != AnalysisControl::XYZ && ibuf != AnalysisControl::YZX && ibuf != AnalysisControl::ZXY ){
				OutputFiles::m_logFile << "Error : Number of parameter specifing the way numbering must be -1, 0, 1 or 2 !!" << std::endl;
				exit(1);
			}
			m_numberingMethod = ibuf;
			hasAlreadyRead[paramID] = true;	
		}else if( line.substr(0,12).compare("OUTPUT_PARAM") == 0 ){
			const int paramID= AnalysisControl::OUTPUT_PARAM;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OUTPUT_PARAM" << std::endl;
				exit(1);
			}			
			int num(0);
			inFile >> num;
			if( num < 0 ){
				OutputFiles::m_logFile << "Error : Number of parameter to be outputed to VTK is less than 0 !!" << std::endl;
				exit(1);
			}else if( num > 0 ){
				for( int i = 0; i < num; ++i ){
					inFile >> ibuf;
					m_outputParametersForVis.insert(ibuf);
				}
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,16).compare("OUTPUT_OPTION") == 0 ){
			const int paramID= AnalysisControl::OUTPUT_OPTION;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OUTPUT_OPTION" << std::endl;
				exit(1);
			}			
			int ibufX(0);
			int ibufY(0);
			inFile >> ibufX >> ibufY;
			if( ibufX == 0 ){
				m_useBackwardOrForwardElement.directionX = AnalysisControl::BACKWARD_ELEMENT;
			}else{
				m_useBackwardOrForwardElement.directionX = AnalysisControl::FORWARD_ELEMENT;
			}
			if( ibufY == 0 ){
				m_useBackwardOrForwardElement.directionY = AnalysisControl::BACKWARD_ELEMENT;
			}else{
				m_useBackwardOrForwardElement.directionY = AnalysisControl::FORWARD_ELEMENT;
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,17).compare("OUTPUT_2D_RESULTS") == 0 ){
			const int paramID= AnalysisControl::OUTPUT_2D_RESULTS;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OUTPUT_2D_RESULTS" << std::endl;
				exit(1);
			}	
			m_isOutput2DResult = true;
			hasAlreadyRead[paramID] = true;
		//}else if( line.substr(0,12).compare("FORWARD_ONLY") == 0 ){
		//	const int paramID= AnalysisControl::FORWARD_ONLY;
		//	if( hasAlreadyRead[paramID] == true ){
		//		OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : FORWARD_ONLY" << std::endl;
		//		exit(1);
		//	}	
		//	m_performForwardOnly = true;
		//	hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,15).compare("TRADE_OFF_PARAM") == 0 ){

			if( !hasAlreadyRead[AnalysisControl::DISTORTION] ){
				OutputFiles::m_logFile << "Error : You must write DISTORTION data above TRADE_OFF_PARAM" << std::endl;
				exit(1);
			}

			const int paramID= AnalysisControl::TRADE_OFF_PARAM;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : TRADE_OFF_PARAM" << std::endl;
				exit(1);
			}	

			switch(m_typeOfDistortion){
				case AnalysisControl::NO_DISTORTION:
					inFile >> m_tradeOffParameterForResistivityValue;			
					break;
				case AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE:
					inFile >> m_tradeOffParameterForResistivityValue >> m_tradeOffParameterForDistortionMatrixComplexity;			
					break;
				case AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS:
					inFile >> m_tradeOffParameterForResistivityValue >> m_tradeOffParameterForDistortionGain >> m_tradeOffParameterForDistortionRotation;			
					break;
				case AnalysisControl::ESTIMATE_GAINS_ONLY:
					inFile >> m_tradeOffParameterForResistivityValue >> m_tradeOffParameterForDistortionGain;			
					break;
				default:
					OutputFiles::m_logFile << "Error : Wrong type of distortion : " << ibuf << std::endl;
					exit(1);
					break;
			}

			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,9).compare("ITERATION") == 0 ){
			const int paramID= AnalysisControl::ITERATION;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ITERATION" << std::endl;
				exit(1);
			}	
			inFile >> m_iterationNumInit >> m_iterationNumMax;			
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,18).compare("DECREASE_THRESHOLD") == 0 ){
			const int paramID= AnalysisControl::DECREASE_THRESHOLD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : DECREASE_THRESHOLD" << std::endl;
				exit(1);
			}	
			inFile >> m_thresholdValueForDecreasing;
			if( m_thresholdValueForDecreasing < 0 ){
				OutputFiles::m_logFile << "Error : Threshold value for determining if objective functional decrease must be positive." << std::endl;
				exit(1);
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,8).compare("CONVERGE") == 0 ){
			const int paramID= AnalysisControl::CONVERGE;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : CONVERGE" << std::endl;
				exit(1);
			}	
			inFile >> m_decreaseRatioForConvegence;		
			if( m_decreaseRatioForConvegence < 0 ){
				OutputFiles::m_logFile << "Error : Criterion for convergence must be positive." << std::endl;
				exit(1);
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,7).compare("RETRIAL") == 0 ){
			const int paramID= AnalysisControl::RETRIAL;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : RETRIAL" << std::endl;
				exit(1);
			}	
			inFile >> m_numCutbackMax;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,11).compare("STEP_LENGTH") == 0 ){
			const int paramID= AnalysisControl::STEP_LENGTH;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : STEP_LENGTH" << std::endl;
				exit(1);
			}	
			inFile >> m_stepLengthDampingFactorCur >> m_stepLengthDampingFactorMin
					>> m_stepLengthDampingFactorMax >> m_numOfIterIncreaseStepLength
					>> m_factorDecreasingStepLength >> m_factorIncreasingStepLength;			
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,10).compare("DISTORTION") == 0 ){
			const int paramID= AnalysisControl::DISTORTION;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : DISTORTION" << std::endl;
				exit(1);
			}	
			inFile >> ibuf;

			if( ibuf != AnalysisControl::NO_DISTORTION && 
				ibuf != AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE &&
				ibuf != AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS &&
				ibuf != AnalysisControl::ESTIMATE_GAINS_ONLY ){
				OutputFiles::m_logFile << "Error : Wrong type ID is specified below DISTORTION : " << ibuf << std::endl;
				exit(1);
				break;
			}

			m_typeOfDistortion = ibuf;

			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,12).compare("ROUGH_MATRIX") == 0 ){
			const int paramID= AnalysisControl::ROUGH_MATRIX;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ROUGH_MATRIX" << std::endl;
				exit(1);
			}	

			inFile >> ibuf;
			if( ibuf >= AnalysisControl::EndOfTypeOfRoughningMatrix ){
				OutputFiles::m_logFile << "Error : Inputted parameter specifing the way of creating roughning matrix is wrong !! : " << ibuf << std::endl;
				exit(1);
			}else{
				m_typeOfRoughningMatrix = ibuf;
			}

			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,10).compare("ELEC_FIELD") == 0 ){
			const int paramID= AnalysisControl::ELEC_FIELD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ELEC_FIELD" << std::endl;
				exit(1);
			}	

			inFile >> ibuf;
			if( ibuf < 0 ){
				m_isTypeOfElectricFieldSetIndivisually = true;
			}
			else if(ibuf != AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD &&
					ibuf != AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD ){
				OutputFiles::m_logFile << "Error : Unknown type of the electric field is specified in ELEC_FIELD : " << ibuf << std::endl;
				exit(1);
			}
			m_typeOfElectricField = ibuf;

			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,15).compare("DIV_NUM_RHS_FWD") == 0 ){
			const int paramID= AnalysisControl::DIV_NUM_RHS_FWD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : DIV_NUM_RHS_FWD" << std::endl;
				exit(1);
			}	

			inFile >> m_divisionNumberOfMultipleRHSInForward;

			if( m_divisionNumberOfMultipleRHSInForward < 1 ){
				OutputFiles::m_logFile << "Error : Division number of right-hand sides must be greater than zero !! Specified number is " << m_divisionNumberOfMultipleRHSInForward << "." << std::endl;
				exit(1);
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,15).compare("DIV_NUM_RHS_INV") == 0 ){
			const int paramID= AnalysisControl::DIV_NUM_RHS_INV;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : DIV_NUM_RHS_INV" << std::endl;
				exit(1);
			}	

			inFile >> m_divisionNumberOfMultipleRHSInInversion;

			if( m_divisionNumberOfMultipleRHSInInversion < 1 ){
				OutputFiles::m_logFile << "Error : Division number of right-hand sides must be greater than zero !! Specified number is " << m_divisionNumberOfMultipleRHSInInversion << "." << std::endl;
				exit(1);
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,18).compare("RESISTIVITY_BOUNDS") == 0 ){
			const int paramID= AnalysisControl::RESISTIVITY_BOUNDS;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : RESISTIVITY_BOUNDS" << std::endl;
				exit(1);
			}	
			inFile >> ibuf;
			ptrResistivityBlock->setTypeBoundConstraints(ibuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,10).compare("OFILE_TYPE") == 0 ){
			const int paramID= AnalysisControl::OFILE_TYPE;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OFILE_TYPE" << std::endl;
				exit(1);
			}	

			inFile >> ibuf;
			
			if( ibuf == 0 ){// ASCII format
				m_binaryOutput = false;
			}
			else{// Binary format
				m_binaryOutput = true;
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,12).compare("HOLD_FWD_MEM") == 0 ){
			const int paramID= AnalysisControl::HOLD_FWD_MEM;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : HOLD_FWD_MEM" << std::endl;
				exit(1);
			}	
			m_holdMemoryForwardSolver = true;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,12).compare("ALPHA_WEIGHT") == 0 ){
			const int paramID= AnalysisControl::ALPHA_WEIGHT;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ALPHA_WEIGHT" << std::endl;
				exit(1);
			}
			for( int iDir = 0; iDir < 3; ++iDir ){
				double dbuf(0.0);
				inFile >> dbuf;
				if(dbuf < 0.0){
					OutputFiles::m_logFile << "Error : Weighting factor of alpha must be positive !! : " << dbuf << std::endl;
					exit(1);
				}
				m_alphaWeight[iDir] = dbuf;
			}

			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,25).compare("INV_MAT_POSITIVE_DEFINITE") == 0 ){
			const int paramID= AnalysisControl::INV_MAT_POSITIVE_DEFINITE;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : INV_MAT_POSITIVE_DEFINITE" << std::endl;
				exit(1);
			}
			m_positiveDefiniteNormalEqMatrix = true;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,18).compare("BOTTOM_RESISTIVITY") == 0 ){
			const int paramID= AnalysisControl::BOTTOM_RESISTIVITY;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : BOTTOM_RESISTIVITY" << std::endl;
				exit(1);
			}

			ptrResistivityBlock->setFlagIncludeBottomResistivity(true);
			double dbuf(0.0);
			inFile >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << "Error : Bottom resistivity is set to be negative !! : " << dbuf << std::endl;
				exit(1);
			}
			ptrResistivityBlock->setBottomResistivity(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,23).compare("BOTTOM_ROUGHNING_FACTOR") == 0 ){
			const int paramID= AnalysisControl::BOTTOM_ROUGHNING_FACTOR;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : BOTTOM_ROUGHNING_FACTOR" << std::endl;
				exit(1);
			}
			inFile >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << "Error : Roughning factor at bottom is set to be negative !! : " << dbuf << std::endl;
				exit(1);
			}
			ptrResistivityBlock->setRoughningFactorAtBottom(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,10).compare("INV_METHOD") == 0 ){
			const int paramID= AnalysisControl::INV_METHOD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : INV_METHOD" << std::endl;
				exit(1);
			}
			inFile >> m_inversionMethod;
			if( m_inversionMethod != Inversion::GAUSS_NEWTON_DATA_SPECE && m_inversionMethod != Inversion::GAUSS_NEWTON_MODEL_SPECE ){
				OutputFiles::m_logFile << "Error : Type of inversion method is wrong  !! : " << m_inversionMethod << std::endl;
				exit(1);
			}
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,16).compare("BOUNDS_DIST_THLD") == 0 ){
			const int paramID= AnalysisControl::BOUNDS_DIST_THLD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : BOUNDS_DIST_THLD" << std::endl;
				exit(1);
			}
			inFile >> dbuf;
			if( dbuf <= 0.0 ){
				OutputFiles::m_logFile << "Error : Minimum distance to resistivity bounds must be positive !!" << std::endl;
				exit(1);
			}
			ptrResistivityBlock->setMinDistanceToBounds(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,3).compare("IDW") == 0 ){
			const int paramID= AnalysisControl::IDW;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : IDW" << std::endl;
				exit(1);
			}
			inFile >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << "Error : Factor of inverse distance weighting must not be negative !!" << std::endl;
				exit(1);
			}
			ptrResistivityBlock->setInverseDistanceWeightingFactor(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,11).compare("SMALL_VALUE") == 0 ){
			const int paramID= AnalysisControl::SMALL_VALUE;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : SMALL_VALUE" << std::endl;
				exit(1);
			}
			inFile >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << "Error : Small value added to the diagonals of roughning matrix must not be negative !!" << std::endl;
				exit(1);
			}
			ptrResistivityBlock->setFlagAddSmallValueToDiagonals(true);
			ptrResistivityBlock->setSmallValueAddedToDiagonals(dbuf);
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,12).compare("MOVE_OBS_LOC") == 0 ){
			const int paramID= AnalysisControl::MOVE_OBS_LOC;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : MOVE_OBS_LOC" << std::endl;
				exit(1);
			}
			m_isObsLocMovedToCenter = true;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,13).compare("OWNER_ELEMENT") == 0 ){
			const int paramID= AnalysisControl::OWNER_ELEMENT;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OWNER_ELEMENT" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			if( ibuf < 0 ){
				m_isTypeOfOwnerElementSetIndivisually = true;
			}
			else if(ibuf != AnalysisControl::USE_LOWER_ELEMENT && ibuf != AnalysisControl::USE_UPPER_ELEMENT ){
				OutputFiles::m_logFile << "Error : Unknown type of owner element is specified in OWNER_ELEMENT : " << ibuf << std::endl;
				exit(1);
			}
			m_typeOfOwnerElement = ibuf;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,14).compare("APP_PHS_OPTION") == 0 ){
			const int paramID= AnalysisControl::APP_PHS_OPTION;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : APP_PHS_OPTION" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			if( ibuf != NO_SPECIAL_TREATMENT_APP_AND_PHASE && ibuf != USE_Z_IF_SIGN_OF_RE_Z_DIFFER ){
				OutputFiles::m_logFile << "Error : Unknown option is specified in APP_PHS_OPTION : " << ibuf << std::endl;
				exit(1);
			}
			m_apparentResistivityAndPhaseTreatmentOption = ibuf;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,19).compare("OUTPUT_ROUGH_MATRIX") == 0 ){
			const int paramID= AnalysisControl::OUTPUT_ROUGH_MATRIX;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : OUTPUT_ROUGH_MATRIX" << std::endl;
				exit(1);
			}
			m_isRougheningMatrixOutputted = true;
			hasAlreadyRead[paramID] = true;
		}else if( line.substr(0,17).compare("DATA_SPACE_METHOD") == 0 ){
			const int paramID= AnalysisControl::DATA_SPACE_METHOD;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : DATA_SPACE_METHOD" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			m_typeOfDataSpaceAlgorithm = ibuf;
			hasAlreadyRead[paramID] = true;
#ifdef _ANISOTOROPY
		}else if( line.substr(0,10).compare("ANISOTROPY") == 0 ){
			const int paramID= AnalysisControl::ANISOTROPY;
			if( hasAlreadyRead[paramID] == true ){
				OutputFiles::m_logFile << "Error : Already read the data from control.dat !! : ANISOTROPY" << std::endl;
				exit(1);
			}
			inFile >> ibuf;
			m_typeOfAnisotropy = ibuf;
			hasAlreadyRead[paramID] = true;
#endif
		}else if( line.substr(0,11).compare("DIFF_FILTER") == 0 ){
			m_useDifferenceFilter = true;
			inFile >> ibuf;
			m_degreeOfLpOptimization = ibuf;
			inFile >> dbuf;
			m_lowerLimitOfDifflog10RhoForLpOptimization = dbuf;
			inFile >> dbuf;
			m_upperLimitOfDifflog10RhoForLpOptimization = dbuf;
			inFile >> ibuf;
			m_maxIterationIRWLSForLpOptimization = ibuf;
			inFile >> dbuf;
			m_thresholdIRWLSForLpOptimization = dbuf;
		}else if (line.substr(0,19).compare("SENSE_MAT_DIRECTORY") == 0) {
			inFile >> m_directoryOfOutOfCoreFilesForSensitivityMatrix;
		}else if( line.substr(0,3).compare("END") == 0 ){
			break;
		}else{
			OutputFiles::m_logFile << "Error : Improper data !! " << line << std::endl;
			exit(1);
		}
		
	}
	
	inFile.close();

	if( !hasAlreadyRead[AnalysisControl::DISTORTION] ){
		OutputFiles::m_logFile << "Error : You must write DISTORTION data in control.dat" << std::endl;
		exit(1);
	}

	if( m_iterationNumMax > m_iterationNumInit && !hasAlreadyRead[AnalysisControl::TRADE_OFF_PARAM] ){
		OutputFiles::m_logFile << "Error : Trade-off parameter must be specified for inversion !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "#==============================================================================" << std::endl;
	OutputFiles::m_logFile << "# Summary of control data" << std::endl;
	OutputFiles::m_logFile << "#==============================================================================" << std::endl;

	if( m_boundaryConditionBottom == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){
		OutputFiles::m_logFile << "# 1D boundary condition is specified at the bottom boundary." << std::endl;
	}else if (m_boundaryConditionBottom == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		OutputFiles::m_logFile << "# Condition of perfect electric conductor is specified at the bottom boundary." << std::endl;
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of boundary condition at the bottom of the model !! m_boundaryConditionBottom = " << m_boundaryConditionBottom << "." << std::endl;
		exit(1);
	}

	if( m_typeOfMesh == MeshData::HEXA ){
		OutputFiles::m_logFile << "# Type of mesh : Hexahedron." << std::endl;
	}else if ( m_typeOfMesh == MeshData::TETRA ){
		OutputFiles::m_logFile << "# Type of mesh : Tetrahedron." << std::endl;
	}else if ( m_typeOfMesh == MeshData::NONCONFORMING_HEXA ){
		OutputFiles::m_logFile << "# Type of mesh : Deformed hexahedron." << std::endl;
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh !! m_typeOfMesh = " << m_typeOfMesh << "." << std::endl;
		exit(1);
	}

	//if( m_orderOfFiniteElement == 0 ){
	//	OutputFiles::m_logFile << "# 0th order elements are used for forward analysis." << std::endl;
	//}else if( m_orderOfFiniteElement == 1 ){
	//	OutputFiles::m_logFile << "# 1st order elements are used for forward analysis." << std::endl;
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! m_orderOfFiniteElement = " << m_orderOfFiniteElement << std::endl;
	//	exit(1);
	//}

	OutputFiles::m_logFile << "# Number of threads is specified to be " << m_numThreads << " ." << std::endl;
	// Specifies the number of threads to use
#ifdef _USE_OMP
	omp_set_num_threads( m_numThreads );
#endif
	mkl_set_num_threads( m_numThreads );

	switch( m_modeOfPARDISO ){
		case PARDISOSolver::INCORE_MODE:
			OutputFiles::m_logFile << "# In-core mode is used." << std::endl;
			break;
		case PARDISOSolver::SELECT_MODE_AUTOMATICALLY:
			OutputFiles::m_logFile << "# Either in-core or out-of-core is used depending on the memory required." << std::endl;
			break;
		case PARDISOSolver::OUT_OF_CORE_MODE:
			OutputFiles::m_logFile << "# Out-of-core mode is used." << std::endl;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong value m_modeOfPARDISO !! m_modeOfPARDISO = " << m_modeOfPARDISO << std::endl;
			exit(1);
			break;
	}

	OutputFiles::m_logFile << "# Division number of right-hand sides at solve phase in forward calculation : " << m_divisionNumberOfMultipleRHSInForward << std::endl;
	OutputFiles::m_logFile << "# Division number of right-hand sides at solve phase in inversion : " << m_divisionNumberOfMultipleRHSInInversion << std::endl;

	if( m_modeOfPARDISO == PARDISOSolver::SELECT_MODE_AUTOMATICALLY || m_modeOfPARDISO == PARDISOSolver::OUT_OF_CORE_MODE ){

#ifdef _LINUX
		std::ostringstream strMem;
		strMem << m_maxMemoryPARDISO;
		if( setenv( "MKL_PARDISO_OOC_MAX_CORE_SIZE", strMem.str().c_str(), 1 ) != 0 ){
			OutputFiles::m_logFile << "Error : Environment variable MKL_PARDISO_OOC_MAX_CORE_SIZE was not set correctly ! " << std::endl;
			exit(1);			
		}
#else
		std::ostringstream strEnv;
		strEnv << "MKL_PARDISO_OOC_MAX_CORE_SIZE=" << m_maxMemoryPARDISO;
#ifdef _DEBUG_WRITE
		std::cout << "strEnv " << strEnv.str() << std::endl;// For debug
#endif
		//if( putenv( strEnv.str().c_str() ) != 0 ){
		//	OutputFiles::m_logFile << "Error : Environment variable MKL_PARDISO_OOC_MAX_CORE_SIZE was not set correctly ! " << std::endl;
		//	exit(1);			
		//}
#endif

		OutputFiles::m_logFile << "# Maximum value of the memory used by out-of-core mode of forward solver : " << m_maxMemoryPARDISO << " [MB]" << std::endl;
	}

	if( m_positiveDefiniteNormalEqMatrix ){
		OutputFiles::m_logFile << "# Coefficient matrix of normal equation is assumed to be positive definite." << std::endl;
	}
	else{
		OutputFiles::m_logFile << "# Coefficient matrix of normal equation is assumed to be indefinite." << std::endl;
	}

	switch( m_numberingMethod ){
	case AnalysisControl::NOT_ASSIGNED:
			OutputFiles::m_logFile << "# Renumbering is not performed." << std::endl;
			break;
		case AnalysisControl::XYZ:
			OutputFiles::m_logFile << "# Numbering edges or nodes in the way X => Y => Z ." << std::endl;
			break;
		case AnalysisControl::YZX:
			OutputFiles::m_logFile << "# Numbering edges or nodes in the way Y => Z => X ." << std::endl;
			break;
		case AnalysisControl::ZXY:
			OutputFiles::m_logFile << "# Numbering edges or nodes in the way Z => X => Y ." << std::endl;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong value m_numberingMethod !! m_numberingMethod = " << m_modeOfPARDISO << std::endl;
			exit(1);
			break;
	}

	if( m_typeOfMesh == MeshData::HEXA ){
		if(	m_typeOfElectricField != AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD ){
			OutputFiles::m_logFile << "Warning : Horizontal electric field is used for hexahedral mesh." << std::endl;
			m_typeOfElectricField = AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD;
		}
	}
	if( m_isTypeOfElectricFieldSetIndivisually ){
		OutputFiles::m_logFile << "# Electric field type of each site is specified in observe.dat." << std::endl;
	}
	else{
		switch( m_typeOfElectricField ){
			case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
				OutputFiles::m_logFile << "# Tangential electric field is used for calculating response functions." << std::endl;
				break;
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				OutputFiles::m_logFile << "# Horizontal electric field is used for calculating response functions." << std::endl;
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of the electric field : " << m_typeOfElectricField << std::endl;
				exit(1);
				break;
		}
	}

	if( m_isTypeOfOwnerElementSetIndivisually ){
		OutputFiles::m_logFile << "# Owner element type of each site is specified in observe.dat." << std::endl;
	}
	else{
		switch(m_typeOfOwnerElement){
			case AnalysisControl::USE_LOWER_ELEMENT:
				OutputFiles::m_logFile << "# EM field is interpolated from the values of the edges of the lower element." << std::endl;
				break;
			case AnalysisControl::USE_UPPER_ELEMENT:
				OutputFiles::m_logFile << "# EM field is interpolated from the values of the edges of the upper element." << std::endl;
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of owner element : " << m_typeOfOwnerElement << std::endl;
				exit(1);
				break;
		}
	}

	switch(m_apparentResistivityAndPhaseTreatmentOption){
		case AnalysisControl::NO_SPECIAL_TREATMENT_APP_AND_PHASE:
			break;
		case AnalysisControl::USE_Z_IF_SIGN_OF_RE_Z_DIFFER:
			OutputFiles::m_logFile << "# Impedance tensor is used instead of apparent resistivity and phase if signs of Re(Z) are different between observed and calculated responses." << std::endl;
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown type of owner element : " << m_typeOfOwnerElement << std::endl;
			exit(1);
			break;
	}

	if( m_typeOfMesh == MeshData::HEXA ){
		if( m_useBackwardOrForwardElement.directionX == AnalysisControl::BACKWARD_ELEMENT ){
			OutputFiles::m_logFile << "# Element of -X direction is used for points locating on boudarieds of elements." << std::endl;
		}else{
			OutputFiles::m_logFile << "# Element of +X direction is used for points locating on boudarieds of elements." << std::endl;
		}
	
		if( m_useBackwardOrForwardElement.directionY == AnalysisControl::BACKWARD_ELEMENT ){
			OutputFiles::m_logFile << "# Element of -Y direction is used for points locating on boudarieds of elements." << std::endl;
		}else{
			OutputFiles::m_logFile << "# Element of +Y direction is used for points locating on boudarieds of elements." << std::endl;
		}
	}

	if( m_isObsLocMovedToCenter ){
		OutputFiles::m_logFile << "# Observation point is moved to the horizontal center of the element including it." << std::endl;
	}

	if( m_holdMemoryForwardSolver ){
		OutputFiles::m_logFile << "# Hold memory of coefficient matrix and sparse solver after forward calculation." << std::endl;
	}
	else{
		OutputFiles::m_logFile << "# Release memory of coefficient matrix and sparse solver after forward calculation." << std::endl;
	}

	const int bountConstraingMethod = ( ResistivityBlock::getInstance() )->getTypeBoundConstraints();
	if(  bountConstraingMethod == ResistivityBlock::SIMPLE_BOUND_CONSTRAINING ){
		OutputFiles::m_logFile << "# Type of bound constraints method : Simple bound constraining" << std::endl;
	}
	else if( bountConstraingMethod == ResistivityBlock::TRANSFORMING_METHOD ){
		OutputFiles::m_logFile << "# Type of bound constraints method : Transforming method" << std::endl;
	}
	else{
		OutputFiles::m_logFile << "Error : Wrong type of bound constraining method !! : " << bountConstraingMethod << " ." << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Minimum distance to resistivity bounds in common logarithm scale : " << 
		ptrResistivityBlock->getMinDistanceToBounds() << " ." << std::endl;

	if( ptrResistivityBlock->includeBottomResistivity() ){
		OutputFiles::m_logFile << "# Bottom resistivity : " << ptrResistivityBlock->getBottomResistivity() << " [Ohm-m]" << std::endl;
		OutputFiles::m_logFile << "# Roughning factor at the bottom : " << ptrResistivityBlock->getRoughningFactorAtBottom() << std::endl;
	}
	else if( ptrResistivityBlock->getFlagAddSmallValueToDiagonals() ){
		if( getTypeOfDataSpaceAlgorithm() == AnalysisControl::NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX ){
			OutputFiles::m_logFile << "# Small value added to the diagonals of [R]T*[R] matrix : " << ptrResistivityBlock->getSmallValueAddedToDiagonals() << std::endl;
		}else{
			OutputFiles::m_logFile << "# Small value added to the diagonals of roughning matrix : " << ptrResistivityBlock->getSmallValueAddedToDiagonals() << std::endl;
		}
	}
	else if( getInversionMethod() == Inversion::GAUSS_NEWTON_DATA_SPECE ){
		OutputFiles::m_logFile << "Error : You must give small number added to diagonals of roughning matrix" << std::endl;
		OutputFiles::m_logFile << "        when data space inverson method is selected !!" << std::endl;
#ifdef _DEBUG_WRITE
#else
			exit(1);
#endif
	}

	if( m_useDifferenceFilter ){
		if( getTypeOfDataSpaceAlgorithm() != AnalysisControl::NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX ){
			OutputFiles::m_logFile << "Error : You must select " << AnalysisControl::NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX << " as DATA_SPACE_METHOD when you use Lp optimization" << std::endl;
			exit(1);
		}
		OutputFiles::m_logFile << "# Degree of Lp optimization : " << m_degreeOfLpOptimization << std::endl;
		OutputFiles::m_logFile << "# Range of difference of log10(rho) for Lp optimization : " << m_lowerLimitOfDifflog10RhoForLpOptimization << " - " << m_upperLimitOfDifflog10RhoForLpOptimization << std::endl;
		OutputFiles::m_logFile << "# Maximum iteration number of IRWLS for Lp optimization : " << m_maxIterationIRWLSForLpOptimization << std::endl;
		OutputFiles::m_logFile << "# Convergence criteria of IRWLS for Lp optimization [%] : " << m_thresholdIRWLSForLpOptimization << std::endl;
	}

#ifdef _ANISOTOROPY
	switch (getTypeOfAnisotropy()){
		case AnalysisControl::NO_ANISOTROPY:
			// No anisotropy => Nothing to do
			break;
		case AnalysisControl::AXIAL_ANISOTROPY:
			OutputFiles::m_logFile << "# Axial anisotropy is considered." << std::endl;
			if( m_typeOfMesh != MeshData::TETRA ){
				OutputFiles::m_logFile << "Error : Axial anisotropys is supported only for tetrahedral mesh !!" << std::endl;
				exit(1);
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong type of anisotropy : " << getTypeOfAnisotropy() << std::endl;
			exit(1);
			break;
	}
#endif

	// Open VTK file
	if( !m_outputParametersForVis.empty() ){
		if( writeBinaryFormat() ){
			OutputFiles::m_logFile << "# Following variables are written to BINARY file." << std::endl;
		}
		else{
			OutputFiles::m_logFile << "# Following variables are written to to ASCII file." << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK) ){
			OutputFiles::m_logFile << "#  - Resistivity" << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK) ){
			OutputFiles::m_logFile << "#  - Electric field" << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK) ){
			OutputFiles::m_logFile << "#  - Magnetic field" << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_CURRENT_DENSITY) ){
			OutputFiles::m_logFile << "#  - Current density" << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_SENSITIVITY) ){
			OutputFiles::m_logFile << "#  - Sensitivity" << std::endl;
		}
		if( doesOutputToVTK( AnalysisControl::OUTPUT_SENSITIVITY_DENSITY) ){
			OutputFiles::m_logFile << "#  - Sensitivity density" << std::endl;
		}
	}

	// Open csv file in which the results of 2D forward computations is written
	if( m_isOutput2DResult ){
		OutputFiles::m_logFile << "# Output results of 2D forward computations to csv file." << std::endl;
		//OutputFiles* const ptrOutputFiles = OutputFiles::getInstance();
		//ptrOutputFiles->openCsvFileFor2DFwd();
	}

	OutputFiles::m_logFile << "# Method of inversion : ";
	switch(getInversionMethod()){
		case Inversion::GAUSS_NEWTON_MODEL_SPECE:
			OutputFiles::m_logFile << "Gauss-newton method (Model space)" << std::endl;
			break;
		case Inversion::GAUSS_NEWTON_DATA_SPECE:
			switch(getTypeOfDataSpaceAlgorithm()){
				case AnalysisControl::NEW_DATA_SPACE_ALGORITHM:
					OutputFiles::m_logFile << "Gauss-newton method (Data space)" << std::endl;
					break;
				case AnalysisControl::NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX:
					OutputFiles::m_logFile << "Gauss-newton method (Data space) using inverse of [R]T*[R] matrix" << std::endl;
					break;
				default:
					OutputFiles::m_logFile << "Error : Type of data space inversion algorithm is wrong  !! : " << getTypeOfDataSpaceAlgorithm() << std::endl;
					exit(1);
					break;
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Type of inversion method is wrong  !! : " << getInversionMethod() << std::endl;
			exit(1);
			break;
	}

	if( estimateDistortionMatrix() ){
		if( m_typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			OutputFiles::m_logFile << "# Components of distortion matrices are estimated directly as model parameters." << std::endl;
		}
		else if( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
			OutputFiles::m_logFile << "# Gains and rotations of distortion matrices are estimated as model parameters." << std::endl;
		}
		else if( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY ){
			OutputFiles::m_logFile << "# Gains of distortion matrices are estimated as model parameters." << std::endl;
		}
	}
	else{
		OutputFiles::m_logFile << "# Distortion matrices are NOT estimated as model parameters." << std::endl;
	}

	if( m_typeOfRoughningMatrix == AnalysisControl::USE_ELEMENTS_SHARE_FACES ){
		OutputFiles::m_logFile << "# Roughening matrix is created using shared faces of elements." << std::endl;
	}else if( m_typeOfRoughningMatrix == AnalysisControl::USER_DEFINED_ROUGHNING ){
		OutputFiles::m_logFile << "# Roughening matrix is created from user-defined roughning factor." << std::endl;
	}else if( m_typeOfRoughningMatrix == AnalysisControl::USE_RESISTIVITY_BLOCKS_SHARE_FACES ){
		OutputFiles::m_logFile << "# Roughening matrix is created using shared faces of resistivity blocks." << std::endl;
	}else if( m_typeOfRoughningMatrix == AnalysisControl::USE_ELEMENTS_SHARE_FACES_AREA_VOL_RATIO ){
		OutputFiles::m_logFile << "# Roughening matrix is created using shared faces of elements (weighting by area-volume ratio)." << std::endl;
	}else{
		OutputFiles::m_logFile << "Error : Number of parameter specifing the way of creating roughning matrix must be 0 , 1 or 2 !!" << std::endl;
		exit(1);
	}

	if( m_isRougheningMatrixOutputted ){
		OutputFiles::m_logFile << "# Roughening matrix is outputted." << std::endl;
	}

	if( m_iterationNumMax < m_iterationNumInit ){
		OutputFiles::m_logFile << "# Inital number of iteration must be less than or equal to the maximum number." << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Trade-off parameter for resistivity value : " << m_tradeOffParameterForResistivityValue << " ." << std::endl;
	if( estimateDistortionMatrix() ){
		if( m_typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			OutputFiles::m_logFile << "# Trade-off parameter for distortion strength : " << m_tradeOffParameterForDistortionMatrixComplexity << " ." << std::endl;
		}
		else if( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
			OutputFiles::m_logFile << "# Trade-off parameter for gains of distortion matrix : " << m_tradeOffParameterForDistortionGain << " ." << std::endl;
			OutputFiles::m_logFile << "# Trade-off parameter for rotations of distortion matrix : " << m_tradeOffParameterForDistortionRotation << " ." << std::endl;
		}
		else if( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY ){
			OutputFiles::m_logFile << "# Trade-off parameter for gains of distortion matrix : " << m_tradeOffParameterForDistortionGain << " ." << std::endl;
		}
	}

	OutputFiles::m_logFile << "# Weighting factor of alpha (X,Y,Z) = (" << m_alphaWeight[0] << ", " << m_alphaWeight[1] << ", " << m_alphaWeight[2] << ") ." << std::endl;

	OutputFiles::m_logFile << "# Factor of inverse distance weighting : " << ptrResistivityBlock->getInverseDistanceWeightingFactor() << "." << std::endl;

	OutputFiles::m_logFile << "# Initial iteration number : " << m_iterationNumInit << "." << std::endl;

	OutputFiles::m_logFile << "# Maximun iteration number : " << m_iterationNumMax << "." << std::endl;

	OutputFiles::m_logFile << "# Threshold value for determining if objective functional decrease : " << m_thresholdValueForDecreasing << "." << std::endl;

	OutputFiles::m_logFile << "# Convergence criterion of inversion is that change ratios of objective function and its components are less than " << m_decreaseRatioForConvegence << " [%] ." << std::endl;

	if( m_stepLengthDampingFactorCur < 0.0 || m_stepLengthDampingFactorCur  > 1.0 ){
		OutputFiles::m_logFile << "Error : Initial factor of step-length damping " << m_stepLengthDampingFactorCur << " must be greater than or equal to zero and less than or equal to one." << std::endl;
		exit(1);
	}
	OutputFiles::m_logFile << "# Initial factor of step-length damping : " << m_stepLengthDampingFactorCur << "." << std::endl;

	if( m_stepLengthDampingFactorCur < 0.0 || m_stepLengthDampingFactorCur  > 1.0 ){
		OutputFiles::m_logFile << "Error : Minimum factor of step-length damping " << m_stepLengthDampingFactorCur << " must be greater than or equal to zero and less than or equal to one." << std::endl;
		exit(1);
	}
	if( m_stepLengthDampingFactorCur < m_stepLengthDampingFactorMin ){
		OutputFiles::m_logFile << "Error : Minimum factor of step-length damping must be less than or equal to the initial one." << std::endl;
		exit(1);
	}
	OutputFiles::m_logFile << "# Minimum factor of step-length damping : " << m_stepLengthDampingFactorMin << "." << std::endl;

	if( m_stepLengthDampingFactorCur < 0.0 || m_stepLengthDampingFactorCur  > 1.0 ){
		OutputFiles::m_logFile << "Error : Maximum factor of step-length damping " << m_stepLengthDampingFactorCur << " must be greater than or equal to zero and less than or equal to one." << std::endl;
		exit(1);
	}
	if( m_stepLengthDampingFactorCur > m_stepLengthDampingFactorMax ){
		OutputFiles::m_logFile << "Error : Maximum factor of step-length damping must be greater than or equal to the initial one." << std::endl;
		exit(1);
	}
	OutputFiles::m_logFile << "# Maximum factor of step-length damping : " << m_stepLengthDampingFactorMax << "." << std::endl;

	if( m_factorDecreasingStepLength < 0 || m_factorDecreasingStepLength > 1.0 ){
		OutputFiles::m_logFile << "Error : Factors of step-length damping is  must be less than or equal to the initial one." << std::endl;
		exit(1);
	}
	OutputFiles::m_logFile << "# If residual increase, factor of step-length damping is muliplied by " << m_factorDecreasingStepLength << " times." << std::endl;

	OutputFiles::m_logFile << "# If residual decrease " << m_numOfIterIncreaseStepLength << " times in a row, factor of step-length damping is muliplied by " << m_factorIncreasingStepLength << " times." << std::endl;

	if( m_numCutbackMax < 0 ){
		m_continueWithoutCutback = true;
		m_numCutbackMax = 0;
		OutputFiles::m_logFile << "# Continue iteration without retrials." << std::endl;
	}else{
		m_continueWithoutCutback = false;
		OutputFiles::m_logFile << "# Maximum number of retrials : " << m_numCutbackMax << "." << std::endl;
	}

	if (!getDirectoryOfOutOfCoreFilesForSensitivityMatrix().empty()) {
		OutputFiles::m_logFile << "# Directory of out-of-core files for the sensitivitry matrix: " << getDirectoryOfOutOfCoreFilesForSensitivityMatrix() << std::endl;
	}

	OutputFiles::m_logFile << "#==============================================================================" << std::endl;
}

// Calculate elapsed time
std::string AnalysisControl::outputElapsedTime() const{

	time_t curTime(NULL);
	time(&curTime);

	std::ostringstream output;
	output << "( " << difftime(curTime, m_startTime) << " sec )";

	return output.str();

}

// Get type of boundary condition at the bottom of the model
int AnalysisControl::getBoundaryConditionBottom() const{
	return m_boundaryConditionBottom;
}

// Get order of finite element
int AnalysisControl::getOrderOfFiniteElement() const{
	return m_orderOfFiniteElement;
}

// Get process ID
int AnalysisControl::getMyPE() const{
	return m_myPE;
}

// Get total number of processes
int AnalysisControl::getTotalPE() const{
	return m_totalPE;
}

// Get total number of threads
int AnalysisControl::getNumThreads() const{
	return m_numThreads;
}

// Get flag specifing either incore or out-of-core version of PARDISO is used
int AnalysisControl::getModeOfPARDISO() const{
	return m_modeOfPARDISO;
}

// Get flag specifing the way of numbering of edges or nodess
int AnalysisControl::getNumberingMethod() const{
	return m_numberingMethod;
}

// Get flag specifing whether the results of 2D forward calculations are outputed
bool AnalysisControl::getIsOutput2DResult() const{
	return m_isOutput2DResult;
}

// Get current iteration number
int AnalysisControl::getIterationNumInit() const{
	return m_iterationNumInit;
}

// Get current iteration number
int AnalysisControl::getIterationNumCurrent() const{
	return m_iterationNumCurrent;
}

// Get maximum iteration number
int AnalysisControl::getIterationNumMax() const{
	return m_iterationNumMax;
}

//// Get Total number of frequencies calculated by this PE
//int AnalysisControl::getNumOfFrequenciesCalculatedByThisPE() const{
//	return m_numOfFrequenciesCalculatedByThisPE;
//}
//		
//// Get IDs of Frequencies calculated by this PE
//int AnalysisControl::getIDsOfFrequenciesCalculatedByThisPE( const int num ) const{
//	
//	if( num >= m_numOfFrequenciesCalculatedByThisPE || num < 0 ){
//		OutputFiles::m_logFile << "Error : num is out range of getIDsOfFrequenciesCalculatedByThisPE !!" << std::endl;
//		exit(1);
//	}
//
//	return m_IDsOfFrequenciesCalculatedByThisPE[num];
//
//}
//
//// Get values of Frequencies calculated by this PE
//double AnalysisControl::getValuesOfFrequenciesCalculatedByThisPE( const int num ) const{
//
//	if( num >= m_numOfFrequenciesCalculatedByThisPE || num < 0 ){
//		OutputFiles::m_logFile << "Error : num is out range of m_valuesOfFrequenciesCalculatedByThisPE !!" << std::endl;
//		exit(1);
//	}
//
//	return m_valuesOfFrequenciesCalculatedByThisPE[num];
//}

// Get member variable specifing which backward or forward element is used for calculating EM field
const AnalysisControl::UseBackwardOrForwardElement AnalysisControl::getUseBackwardOrForwardElement() const{
	return m_useBackwardOrForwardElement;
}

// Get whether the specified parameter is outputed to VTK file
bool AnalysisControl::doesOutputToVTK( const int paramID ) const{

	if( m_outputParametersForVis.find( paramID ) == m_outputParametersForVis.end() ){
		return false;
	}else{
		return true;
	}

}

// Get damping factor for resistivity value
double AnalysisControl::getTradeOffParameterForResistivityValue() const{

	return m_tradeOffParameterForResistivityValue;

}

// Get trade-off parameter for distortion matrix complexity
double AnalysisControl::getTradeOffParameterForDistortionMatrixComplexity() const{

	assert( m_typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE );
	return m_tradeOffParameterForDistortionMatrixComplexity; 

}

// Get trade-off parameter for gains of distortion matrix
double AnalysisControl::getTradeOffParameterForGainsOfDistortionMatrix() const{
	assert( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS || m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY );
	return m_tradeOffParameterForDistortionGain;
}

// Get trade-off parameter for rotations of distortion matrix
double AnalysisControl::getTradeOffParameterForRotationsOfDistortionMatrix() const{
	assert( m_typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS );
	return m_tradeOffParameterForDistortionRotation;
}

// Get current factor of step length damping
double AnalysisControl::getStepLengthDampingFactorCur() const{

	return m_stepLengthDampingFactorCur;

}

// Get maximum number of cutbacks.
int AnalysisControl::getNumCutbackMax() const{

	return m_numCutbackMax;

}

// Get flag whether memory of solver is held after forward calculation
bool AnalysisControl::holdMemoryForwardSolver() const{
	
	return m_holdMemoryForwardSolver;

}

// Get type of mesh
int AnalysisControl::getTypeOfMesh() const{

	return m_typeOfMesh;

}

// Get flag specifing whether distortion matrix is estimated or not
bool AnalysisControl::estimateDistortionMatrix() const{

	return ( m_typeOfDistortion != AnalysisControl::NO_DISTORTION );

}

// Get type of galvanic distortion
int AnalysisControl::getTypeOfDistortion() const{

	return m_typeOfDistortion;

}

// Get flag specifing the way of creating roughning matrix
int AnalysisControl::geTypeOfRoughningMatrix() const{

	return m_typeOfRoughningMatrix;

}

// Get type of the electric field used to calculate response functions
int AnalysisControl::getTypeOfElectricField() const{

	return m_typeOfElectricField;

}

// Flag specifing whether type of the electric field of each site is specified indivisually
bool AnalysisControl::isTypeOfElectricFieldSetIndivisually() const{

	return m_isTypeOfElectricFieldSetIndivisually;

}

// Tyep of owner element of observation sites
int AnalysisControl::getTypeOfOwnerElement() const{

	return m_typeOfOwnerElement;

}

// Flag specifing whether the type of owner element of each site is specified indivisually
bool AnalysisControl::isTypeOfOwnerElementSetIndivisually() const{

	return m_isTypeOfOwnerElementSetIndivisually;

}

// Get division number of right-hand sides at solve phase in forward calculation
int AnalysisControl::getDivisionNumberOfMultipleRHSInForward() const{

	return m_divisionNumberOfMultipleRHSInForward;

}

// Get division number of right-hand sides at solve phase in inversion
int AnalysisControl::getDivisionNumberOfMultipleRHSInInversion() const{

	return m_divisionNumberOfMultipleRHSInInversion;

}

// Get weighting factor of alpha
double AnalysisControl::getAlphaWeight(const int iDir) const{
	assert( iDir >= 0 && iDir < 3 );
	return m_alphaWeight[iDir];
}

// Get flag specifing whether the cofficient matrix of the normal equation is positive definite or not
bool AnalysisControl::getPositiveDefiniteNormalEqMatrix() const{
	return m_positiveDefiniteNormalEqMatrix;
}

// Get flag specifing whether output file for paraview is binary or ascii
bool AnalysisControl::writeBinaryFormat() const{

	return m_binaryOutput;

}

// Get inversion method
int AnalysisControl::getInversionMethod() const{
	return m_inversionMethod;
}

// Get flag specifing whether observation point is moved to the horizontal center of the element including it
int AnalysisControl::getIsObsLocMovedToCenter() const{
	return m_isObsLocMovedToCenter;
}

// Get option about treatment of apparent resistivity & phase
int AnalysisControl::getApparentResistivityAndPhaseTreatmentOption() const{
	return m_apparentResistivityAndPhaseTreatmentOption;
}

// Get flag specifing whether roughening matrix is outputed
bool AnalysisControl::getIsRougheningMatrixOutputted() const{
	return m_isRougheningMatrixOutputted;
}

// Get type of data space algorithm
int AnalysisControl::getTypeOfDataSpaceAlgorithm() const{
	return m_typeOfDataSpaceAlgorithm;
}

// Get flag specifing whether Lp optimization with difference filter is used
bool AnalysisControl::useDifferenceFilter() const{
	return m_useDifferenceFilter;
}

// Get degree of Lp optimization
int AnalysisControl::getDegreeOfLpOptimization() const{
	return m_degreeOfLpOptimization;
}

// Get lower limit of the difference of log10(rho) for Lp optimization
double AnalysisControl::getLowerLimitOfDifflog10RhoForLpOptimization() const{
	return m_lowerLimitOfDifflog10RhoForLpOptimization;
}

// Get upper limit of the difference of log10(rho) for Lp optimization
double AnalysisControl::getUpperLimitOfDifflog10RhoForLpOptimization() const{
	return m_upperLimitOfDifflog10RhoForLpOptimization;
}

// Get maximum iteration number of IRWLS for Lp optimization
int AnalysisControl::getMaxIterationIRWLSForLpOptimization() const{
	return m_maxIterationIRWLSForLpOptimization;
}

// Get threshold value for deciding convergence about IRWLS for Lp optimization
double AnalysisControl::getThresholdIRWLSForLpOptimization() const{
	return m_thresholdIRWLSForLpOptimization;
}

// Get directory of out-of-core files for the sensitivitry matrix
std::string AnalysisControl::getDirectoryOfOutOfCoreFilesForSensitivityMatrix() const {
	return m_directoryOfOutOfCoreFilesForSensitivityMatrix;
}

#ifdef _ANISOTOROPY
// Get type of anisotropy
int AnalysisControl::getTypeOfAnisotropy() const{
	return m_typeOfAnisotropy;
}

// Get flag specifing whether anisotropy of conductivity is taken into account
bool AnalysisControl::isAnisotropyConsidered() const{
	if( getTypeOfAnisotropy() == AnalysisControl::NO_ANISOTROPY ){
		return false;
	}else{
		return true;
	}
}
#endif

// Get pointer to the object of class MeshData
const MeshData* AnalysisControl::getPointerOfMeshData() const{

	if( getPointerOfForward3D() == NULL ){
		OutputFiles::m_logFile << "Error : Pointer to the class Forward3D is NULL." << std::endl;
		exit(1);
	}

	return getPointerOfForward3D()->getPointerToMeshData();

}

// Get pointer to the object of class MeshDataBrickElement
const MeshDataBrickElement* AnalysisControl::getPointerOfMeshDataBrickElement() const{

	if( m_ptrForward3DBrickElement0thOrder == NULL ){
		OutputFiles::m_logFile << "Error : m_ptrForward3DBrickElement0thOrder is NULL." << std::endl;
		exit(1);
	}

	return m_ptrForward3DBrickElement0thOrder->getPointerToMeshDataBrickElement();

}

// Get pointer to the object of class MeshDataTetraElement
const MeshDataTetraElement* AnalysisControl::getPointerOfMeshDataTetraElement() const{

	if( m_ptrForward3DTetraElement0thOrder == NULL ){
		OutputFiles::m_logFile << "Error : m_ptrForward3DTetraElement0thOrder is NULL." << std::endl;
		exit(1);
	}

	return m_ptrForward3DTetraElement0thOrder->getPointerToMeshDataTetraElement();

}

// Get pointer to the object of class MeshDataNonConformingHexaElement
const MeshDataNonConformingHexaElement* AnalysisControl::getPointerOfMeshDataNonConformingHexaElement() const{

	if( m_ptrForward3DNonConformingHexaElement0thOrder == NULL ){
		OutputFiles::m_logFile << "Error : m_ptrForward3DNonConformingHexaElement0thOrder is NULL." << std::endl;
		exit(1);
	}

	return m_ptrForward3DNonConformingHexaElement0thOrder->getPointerToMeshDataNonConformingHexaElement();

}

// Calculate forward computation
void AnalysisControl::calcForwardComputation( const int iter ){

	Forward3D* ptrForward3D = getPointerOfForward3D();
	if( ptrForward3D == NULL ){
		OutputFiles::m_logFile << "Error : Pointer to the class Forward3D is NULL." << std::endl;
		exit(1);
	}

	if( m_ptrInversion == NULL ){
		OutputFiles::m_logFile << "Error : m_ptrInversion is NULL." << std::endl;
		exit(1);
	}

	ObservedData* const pObservedData = ObservedData::getInstance();	

	const int numOfFrequenciesCalculatedByThisPE = pObservedData->getNumOfFrequenciesCalculatedByThisPE();
	for( int ifreq = 0; ifreq < numOfFrequenciesCalculatedByThisPE; ++ifreq ){

		//const int ifreq = m_IDsOfFrequenciesCalculatedByThisPE[ifreq];
		const double frquencyValue = pObservedData->getValuesOfFrequenciesCalculatedByThisPE( ifreq );

		for( int iPol = 0; iPol < 2; ++iPol ){

			std::string polarizationName;
			if( iPol == 0 ){
				polarizationName = "Ex-polarization";
			}else{
				polarizationName = "Ey-polarization";
			}

			OutputFiles::m_logFile << "#================================================================================================" << std::endl;
			OutputFiles::m_logFile << "# Start forward calculation. Frequency : " << frquencyValue << " [Hz], Polarization : " << polarizationName << std::endl;
			OutputFiles::m_logFile << "#================================================================================================" << std::endl;

			ptrForward3D->forwardCalculation( frquencyValue, iPol );

			pObservedData->calculateEMFieldOfAllStations( ptrForward3D, frquencyValue, iPol, ifreq );

			if( doesCalculateSensitivity(iter) ){
				m_ptrInversion->calculateDerivativesOfEMField( ptrForward3D, frquencyValue, iPol );
			}

		}

		OutputFiles::m_logFile << "#==============================================================================" << std::endl;
		OutputFiles::m_logFile << "# Calculate response functions. " << outputElapsedTime() << std::endl;
		pObservedData->calculateResponseFunctionOfAllStations( ifreq );
		if( doesCalculateSensitivity(iter) ){
			m_ptrInversion->calculateSensitivityMatrix( ifreq, frquencyValue );
		}
		OutputFiles::m_logFile << "#==============================================================================" << std::endl;

	}

	OutputFiles::m_logFile << "# Release memory of coefficient matrix and sparse solver. " << outputElapsedTime() << std::endl;
	if(!m_holdMemoryForwardSolver){// Release memory of sparse solver
		ptrForward3D->releaseMemoryOfMatrixAndSolver();
	}

}

// Adjust factor of step length damping and output convergence data to cnv file
AnalysisControl::ConvergenceBehaviors AnalysisControl::adjustStepLengthDampingFactor( const int iterCur, const int iCutbackCur ){

	// Get process ID
	//MPI_Comm_rank( MPI_COMM_WORLD, &myProcessID );
	const int myProcessID = getMyPE();

	ObservedData* const pObservedData = ObservedData::getInstance();

	//const double dataMisfit = pObservedData->calculateErrorSumOfSquares();
	double dataMisfitThisPE = pObservedData->calculateErrorSumOfSquaresThisPE();

#ifdef _DEBUG_WRITE
	std::cout << "PE dataMisfitThisPE : " << myProcessID << " " << dataMisfitThisPE << std::endl;// For debug
#endif

	double dataMisfit(0.0);
	MPI_Reduce( &dataMisfitThisPE, &dataMisfit, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

#ifdef _DEBUG_WRITE
	if( myProcessID == 0 ){// Zero process only ---------------------
		std::cout << "dataMisfit = " << dataMisfit << std::endl;// For debug
	}
#endif

	int iynConverged(0);
	int iynGoNextIteration(0);

	int numDataThisPE = pObservedData->getNumObservedDataThisPETotal();
	int numDataTotal(0);
	MPI_Reduce( &numDataThisPE, &numDataTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

	if( myProcessID == 0 ){// Zero process only ---------------------
	
		const double rms = sqrt( dataMisfit / static_cast<double>( numDataTotal ) );

		double modelRoughness(0.0);
		if( useDifferenceFilter() ){
			modelRoughness = ( ResistivityBlock::getInstance() )->calcModelRoughnessForDifferenceFilter();
		}else{
			modelRoughness = ( ResistivityBlock::getInstance() )->calcModelRoughnessForLaplacianFilter();
		}
		const double modelRoughnessMultipliedByAlphaAlpha = modelRoughness * pow(m_tradeOffParameterForResistivityValue, 2);

		double objectFunctionalCur = dataMisfit	+ modelRoughnessMultipliedByAlphaAlpha;

		double distortionMatrixNorm = -1.0;
		double normOfGains = -1.0;
		double normOfRotations = -1.0;

		//----------------------------------------
		// Output convergence data to cnv file
		//----------------------------------------
		if( !OutputFiles::m_cnvFile.is_open() ){
			OutputFiles::m_logFile << "Error : CNV file has not been opened." << std::endl;
			exit(1);
		}

		if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			distortionMatrixNorm = pObservedData->calculateSumSquareOfDistortionMatrixComplexity();
			objectFunctionalCur += m_tradeOffParameterForDistortionMatrixComplexity * m_tradeOffParameterForDistortionMatrixComplexity * distortionMatrixNorm;

			OutputFiles::m_cnvFile.precision(6);
			OutputFiles::m_cnvFile << std::setw(10) << iterCur << std::setw(10) << iCutbackCur
				<< std::setw(15) << std::scientific << m_tradeOffParameterForResistivityValue
				<< std::setw(15) << std::scientific << m_tradeOffParameterForDistortionMatrixComplexity
				<< std::setw(15) << std::scientific << m_stepLengthDampingFactorCur
				<< std::setw(15) << std::scientific << modelRoughness
				<< std::setw(15) << std::scientific << distortionMatrixNorm
				<< std::setw(15) << std::scientific << dataMisfit
				<< std::setw(15) << std::scientific << rms
				<< std::setw(15) << std::scientific << objectFunctionalCur
				<< std::endl;	

		}
		else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
			normOfGains = pObservedData->calculateSumSquareOfDistortionMatrixGains();
			objectFunctionalCur += m_tradeOffParameterForDistortionGain * m_tradeOffParameterForDistortionGain * normOfGains;
			normOfRotations = pObservedData->calculateSumSquareOfDistortionMatrixRotations();
			objectFunctionalCur += m_tradeOffParameterForDistortionRotation * m_tradeOffParameterForDistortionRotation * normOfRotations;

			OutputFiles::m_cnvFile.precision(6);
			OutputFiles::m_cnvFile << std::setw(10) << iterCur << std::setw(10) << iCutbackCur
				<< std::setw(15) << std::scientific << m_tradeOffParameterForResistivityValue
				<< std::setw(15) << std::scientific << m_tradeOffParameterForDistortionGain
				<< std::setw(15) << std::scientific << m_tradeOffParameterForDistortionRotation
				<< std::setw(15) << std::scientific << m_stepLengthDampingFactorCur
				<< std::setw(15) << std::scientific << modelRoughness
				<< std::setw(15) << std::scientific << normOfGains
				<< std::setw(15) << std::scientific << normOfRotations
				<< std::setw(15) << std::scientific << dataMisfit
				<< std::setw(15) << std::scientific << rms
				<< std::setw(15) << std::scientific << objectFunctionalCur
				<< std::endl;	

		}
		else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
			normOfGains = pObservedData->calculateSumSquareOfDistortionMatrixGains();
			objectFunctionalCur += m_tradeOffParameterForDistortionGain * m_tradeOffParameterForDistortionGain * normOfGains;

			OutputFiles::m_cnvFile.precision(6);
			OutputFiles::m_cnvFile << std::setw(10) << iterCur << std::setw(10) << iCutbackCur
				<< std::setw(15) << std::scientific << m_tradeOffParameterForResistivityValue
				<< std::setw(15) << std::scientific << m_tradeOffParameterForDistortionGain
				<< std::setw(15) << std::scientific << m_stepLengthDampingFactorCur
				<< std::setw(15) << std::scientific << modelRoughness
				<< std::setw(15) << std::scientific << normOfGains
				<< std::setw(15) << std::scientific << dataMisfit
				<< std::setw(15) << std::scientific << rms
				<< std::setw(15) << std::scientific << objectFunctionalCur
				<< std::endl;
		}
		else{

			OutputFiles::m_cnvFile.precision(6);
			OutputFiles::m_cnvFile << std::setw(10) << iterCur << std::setw(10) << iCutbackCur
				<< std::setw(15) << std::scientific << m_tradeOffParameterForResistivityValue
				<< std::setw(15) << std::scientific << m_stepLengthDampingFactorCur
				<< std::setw(15) << std::scientific << modelRoughness
				<< std::setw(15) << std::scientific << dataMisfit
				<< std::setw(15) << std::scientific << rms
				<< std::setw(15) << std::scientific << objectFunctionalCur
				<< std::endl;			

		}		

		// Perform convergence test
		if( checkConvergence( objectFunctionalCur ) ){
			iynConverged = 1;
		}

		//----------------------------------------
		// Adjust factor of step length damping
		//----------------------------------------

		if( iterCur <= getIterationNumInit() || m_continueWithoutCutback ){// First iteration or no cutback

			m_objectFunctionalPre = objectFunctionalCur;
			//goNextIteration = true;// Go next iteration
			iynGoNextIteration = 1;// Go next iteration

		}else{

			const double stepLengthDampingFactorPre = m_stepLengthDampingFactorCur;

			if( objectFunctionalCur < m_objectFunctionalPre - m_thresholdValueForDecreasing ){
				// Value of objective functional decrease from the one of previous iteration

				//goNextIteration = true;
				iynGoNextIteration = 1;// Go next iteration

				if( ++m_numConsecutiveIterFunctionalDecreasing >= m_numOfIterIncreaseStepLength ){
					m_stepLengthDampingFactorCur *= m_factorIncreasingStepLength;// Increase factor of step length damping
				}

				m_objectFunctionalPre = objectFunctionalCur;

				m_dataMisfitPre = dataMisfit;
				m_modelRoughnessPre = modelRoughness;
				m_normOfDistortionMatrixDifferencesPre = distortionMatrixNorm;
				m_normOfGainsPre = normOfGains;
				m_normOfRotationsPre = normOfRotations;

			}else{
				// Value of objective functional increase from the one of previous iteration

				//goNextIteration = false;
				iynGoNextIteration = 0;// Not go next iteration

				m_numConsecutiveIterFunctionalDecreasing = 0;// reset value

				m_stepLengthDampingFactorCur *= m_factorDecreasingStepLength;// Decreaes factor of step length damping

			}

			if( m_stepLengthDampingFactorCur < m_stepLengthDampingFactorMin ){// Reach minimum value
				m_stepLengthDampingFactorCur = m_stepLengthDampingFactorMin;
			}

			if( m_stepLengthDampingFactorCur > m_stepLengthDampingFactorMax ){// Reach maximum value
				m_stepLengthDampingFactorCur = m_stepLengthDampingFactorMax;
			}

			const double threshold = 1.0E-8;
			if( fabs( m_stepLengthDampingFactorCur - stepLengthDampingFactorPre ) > threshold ){

				OutputFiles::m_logFile << "# Factor of step length damping change from " << stepLengthDampingFactorPre << " to " << m_stepLengthDampingFactorCur << "." << std::endl;

			}

		}

	}//--------------------------------------------------------------
	
	MPI_Bcast( &iynConverged, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &iynGoNextIteration, 1, MPI_INT, 0, MPI_COMM_WORLD );

	MPI_Bcast( &m_objectFunctionalPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_dataMisfitPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_modelRoughnessPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_normOfDistortionMatrixDifferencesPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_normOfGainsPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_normOfRotationsPre, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	MPI_Bcast( &m_stepLengthDampingFactorCur, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &m_numConsecutiveIterFunctionalDecreasing, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

#ifdef _DEBUG_WRITE
	std::cout << "PE iynConverged : " <<  myProcessID << " " << iynConverged << std::endl;// For debug
	std::cout << "PE iynGoNextIteration : " <<  myProcessID << " " << iynGoNextIteration << std::endl;// For debug
	std::cout << "PE m_objectFunctionalPre : " <<  myProcessID << " " << m_objectFunctionalPre << std::endl;// For debug
	std::cout << "PE m_dataMisfitPre : " <<  myProcessID << " " << m_dataMisfitPre << std::endl;// For debug
	std::cout << "PE m_modelRoughnessPre : " <<  myProcessID << " " << m_modelRoughnessPre << std::endl;// For debug
	std::cout << "PE m_normOfDistortionMatrixDifferencesPre : " <<  myProcessID << " " << m_normOfDistortionMatrixDifferencesPre << std::endl;// For debug
	std::cout << "PE m_normOfGainsPre : " <<  myProcessID << " " << m_normOfGainsPre << std::endl;// For debug
	std::cout << "PE m_normOfRotationsPre : " <<  myProcessID << " " << m_normOfRotationsPre << std::endl;// For debug
	std::cout << "PE m_stepLengthDampingFactorCur : " <<  myProcessID << " " << m_stepLengthDampingFactorCur << std::endl;// For debug
	std::cout << "PE m_numConsecutiveIterFunctionalDecreasing : " <<  myProcessID << " " << m_numConsecutiveIterFunctionalDecreasing << std::endl;// For debug
#endif

	if( iynConverged == 1 ){
		return AnalysisControl::INVERSIN_CONVERGED;
	}
	else{
		if( iynGoNextIteration == 1 ){
			return AnalysisControl::GO_TO_NEXT_ITERATION;
		}
		else{
			return AnalysisControl::DURING_RETRIALS;
		}
	}

}

bool AnalysisControl::checkConvergence( const double objectFunctionalCur ){

	const double criterion = m_decreaseRatioForConvegence * 0.01;

	if( m_objectFunctionalPre - objectFunctionalCur > 0.0 &&
		fabs( m_objectFunctionalPre - objectFunctionalCur ) < m_objectFunctionalPre * criterion ){
		return true;
	}

	return false;
}

// Perform convergence test 
bool AnalysisControl::checkConvergence( const double objectFunctionalCur, const double dataMisft, const double modelRoughness,
										const double normDist1, const double normDist2 ){

	const double criterion = m_decreaseRatioForConvegence * 0.01;

	if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		if( ( m_objectFunctionalPre - objectFunctionalCur ) > 0.0 &&
			fabs( m_objectFunctionalPre - objectFunctionalCur ) < m_objectFunctionalPre * criterion &&
			fabs( m_dataMisfitPre - dataMisft ) < m_dataMisfitPre * criterion &&
			fabs( m_modelRoughnessPre - modelRoughness ) < m_modelRoughnessPre * criterion &&
			fabs( m_normOfDistortionMatrixDifferencesPre - normDist1 ) < m_normOfDistortionMatrixDifferencesPre * criterion ){
			return true;
		}
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		if( ( m_objectFunctionalPre - objectFunctionalCur ) > 0.0 &&
			fabs( m_objectFunctionalPre - objectFunctionalCur ) < m_objectFunctionalPre * criterion &&
			fabs( m_dataMisfitPre - dataMisft ) < m_dataMisfitPre * criterion &&
			fabs( m_modelRoughnessPre - modelRoughness ) < m_modelRoughnessPre * criterion &&
			fabs( m_normOfGainsPre - normDist1 ) < m_normOfGainsPre * criterion &&
			fabs( m_normOfRotationsPre - normDist2 ) < m_normOfRotationsPre * criterion ){
			return true;
		}
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		if( ( m_objectFunctionalPre - objectFunctionalCur ) > 0.0 &&
			fabs( m_objectFunctionalPre - objectFunctionalCur ) < m_objectFunctionalPre * criterion &&
			fabs( m_dataMisfitPre - dataMisft ) < m_dataMisfitPre * criterion &&
			fabs( m_modelRoughnessPre - modelRoughness ) < m_modelRoughnessPre * criterion &&
			fabs( m_normOfGainsPre - normDist1 ) < m_normOfGainsPre * criterion ){
			return true;
		}
	}
	else{
		if( ( m_objectFunctionalPre - objectFunctionalCur ) > 0.0 &&
			fabs( m_objectFunctionalPre - objectFunctionalCur ) < m_objectFunctionalPre * criterion &&
			fabs( m_dataMisfitPre - dataMisft ) < m_dataMisfitPre * criterion &&
			fabs( m_modelRoughnessPre - modelRoughness ) < m_modelRoughnessPre * criterion ){
			return true;
		}
	}

	return false;
}

// Return flag specifing whether sensitivity is calculated or not
bool AnalysisControl::doesCalculateSensitivity( const int iter ) const{

	return ( m_iterationNumMax > iter ) ? true : false;

}

// Get pointer to the object of class Forward3D
Forward3D* AnalysisControl::getPointerOfForward3D() const{

	if( m_typeOfMesh == MeshData::HEXA ){

		if( m_ptrForward3DBrickElement0thOrder != NULL ){
			return static_cast<Forward3D*>( m_ptrForward3DBrickElement0thOrder );
		}else{
			OutputFiles::m_logFile << "Error : m_ptrForward3DBrickElement0thOrderv is NULL." << std::endl;
			exit(1);
		}

	}else if( m_typeOfMesh == MeshData::TETRA ){

		if( m_ptrForward3DTetraElement0thOrder != NULL ){
			return static_cast<Forward3D*>( m_ptrForward3DTetraElement0thOrder );
		}else{
			OutputFiles::m_logFile << "Error : m_ptrForward3DTetraElement0thOrder is NULL." << std::endl;
			exit(1);
		}

	}else if( m_typeOfMesh == MeshData::NONCONFORMING_HEXA ){

		if( m_ptrForward3DNonConformingHexaElement0thOrder != NULL ){
			return static_cast<Forward3D*>( m_ptrForward3DNonConformingHexaElement0thOrder );
		}else{
			OutputFiles::m_logFile << "Error : m_ptrForward3DNonConformingHexaElement0thOrder is NULL." << std::endl;
			exit(1);
		}

	}else{
		OutputFiles::m_logFile << "Error : Type of mesh is wrong !! : " << m_typeOfMesh << "." << std::endl;
		exit(1);
	}

	return NULL;

}
