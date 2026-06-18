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
#ifndef DBLDEF_ANALYSIS_CONTROL
#define DBLDEF_ANALYSIS_CONTROL

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <set>
#include "Inversion.h"
#include "Forward3D.h"
#include "Forward3DBrickElement0thOrder.h"
#include "Forward3DTetraElement0thOrder.h"
#include "Forward3DNonConformingHexaElement0thOrder.h"
#include "MeshData.h"
#include "ResistivityBlockIsotropic.h"
#include "ResistivityBlockAnisotropic.h"

// Class of Analysis Control
class AnalysisControl{

	public:
		// Type of boundary condition at the bottom of the model
		static const int BOUNDARY_BOTTOM_ONE_DIMENSIONAL       = 0; // Boundary condition on which one dimensional relation holds for the EM field at the bottom
		static const int BOUNDARY_BOTTOM_PERFECT_CONDUCTOR     = 1; // Boundary condition on which electric field is null at the bottom

		// IDs of output parameter which can be outputed to VTK file
		enum outputParameterIDsForVTK{
			OUTPUT_RESISTIVITY_VALUES_TO_VTK = 0,
			OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK,
			OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK,
			OUTPUT_CURRENT_DENSITY,
			OUTPUT_SENSITIVITY,
			OUTPUT_SENSITIVITY_DENSITY,
		};

		// The way of numbering
		enum numbering{
			NOT_ASSIGNED = -1,
			XYZ=0,
			YZX=1,
			ZXY=2,
		};

		// Flags specifing which backward or forward element is used for calculating EM field
		enum BackwardOrForwardElement{
			BACKWARD_ELEMENT = 0,
			FORWARD_ELEMENT,
		};

		// Flags specifing the way of creating roughning matrix
		enum TypeOfRoughningMatrix{
			USER_DEFINED_ROUGHNING = -1,
			USE_ELEMENTS_SHARE_FACES = 0,
			USE_RESISTIVITY_BLOCKS_SHARE_FACES,
			USE_ELEMENTS_SHARE_FACES_AREA_VOL_RATIO,
			EndOfTypeOfRoughningMatrix// This must be written at the end
		};

		// Type of galvanic distortion
		enum TypeOfDistortion{
			DISTORTION_TYPE_UNDEFINED = -1,
			NO_DISTORTION = 0,
			ESTIMATE_DISTORTION_MATRIX_DIFFERENCE,
			ESTIMATE_GAINS_AND_ROTATIONS,
			ESTIMATE_GAINS_ONLY,
		};

		// Flag of convergence behaviors
		enum ConvergenceBehaviors{
			DURING_RETRIALS = 0,
			GO_TO_NEXT_ITERATION,
			INVERSIN_CONVERGED,
		};
		
		struct UseBackwardOrForwardElement{
			enum BackwardOrForwardElement directionX;
			enum BackwardOrForwardElement directionY;
		};

		// Type of electric field used to calculate response functions
		enum TypeOfElectricField{
			USE_HORIZONTAL_ELECTRIC_FIELD = 0,
			USE_TANGENTIAL_ELECTRIC_FIELD = 1,
		};
		
		// Type of owner element
		enum TypeOfOwnerElement{
			USE_LOWER_ELEMENT = 0,
			USE_UPPER_ELEMENT = 1,
		};

		// Option about treatment of apparent resistivity & phase
		enum AppResPhaseTreatmentOption{
			NO_SPECIAL_TREATMENT_APP_AND_PHASE = 0,
			USE_Z_IF_SIGN_OF_RE_Z_DIFFER = 1,
		};

		// Type of data space algorithm
		enum TypeOfDataSpaceAlgorithm{
			NEW_DATA_SPACE_ALGORITHM = 1,
			NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX = 2,
		};

		// Options of GCV calculation
		enum OptionOfGCVCalculation {
			GCV_WITHOUT_USING_DAMPING_FACTOR = 0,
			GCV_CONSIDERING_DAMPING_FACTOR,
			END_OF_OPTION_OF_GCV_CALCULATION
		};

		// Return the the instance of the class
		static AnalysisControl* getInstance();

		// Run analysis
		void run();

		// Calculate and output elapsed time
		std::string outputElapsedTime() const;

		// Get type of boundary condition at the bottom of the model
		int getBoundaryConditionBottom() const;

		// Get order of finite element
		int getOrderOfFiniteElement() const;

		// Get process ID
		int getMyPE() const;

		// Get total number of processes
		int getTotalPE() const;

		// Get total number of threads
		int getNumThreads() const;

		// Get flag specifing either incore or out-of-core version of PARDISO is used
		int getModeOfPARDISO() const;

		// Get flag specifing the way of numbering of edges or nodess
		int getNumberingMethod() const;

		// Get flag specifing whether the results of 2D forward calculations are outputed
		bool getIsOutput2DResult() const;

		// Get initial iteration number
		int getIterationNumInit() const;
		
		// Get current iteration number
		int getIterationNumCurrent() const;

		// Get maximum iteration number
		int getIterationNumMax() const;

		// Get member variable specifing which backward or forward element is used for calculating EM field
		const UseBackwardOrForwardElement getUseBackwardOrForwardElement() const;

		// Get whether the specified parameter is outputed to VTK file
		bool doesOutputToVTK( const int paramID ) const;

		// Get trade-off parameter for the spatial variation of resistivity value
		double getTradeOffParameterForResistivityValue() const;

		// Get trade-off parameter for the degrees of anistropy
		double getTradeOffParameterForDegreeOfAnisotropy() const;

		// Get trade-off parameter for the spatial variation of the anistropy direction
		double getTradeOffParameterForAnisotropyDirection() const;

		// Get trade-off parameter for distortion matrix complexity
		double getTradeOffParameterForDistortionMatrixComplexity() const;

		// Get trade-off parameter for gains of distortion matrix
		double getTradeOffParameterForGainsOfDistortionMatrix() const;

		// Get trade-off parameter for rotations of distortion matrix
		double getTradeOffParameterForRotationsOfDistortionMatrix() const;

		// Get current factor of step length damping
		double getStepLengthDampingFactorCur() const;

		// Get maximum number of cutbacks.
		int getNumCutbackMax() const;

		// Get flag whether memory of solver is held after forward calculation
		bool holdMemoryForwardSolver() const;

		// Get type of mesh
		int getTypeOfMesh() const;

		// Get flag specifing whether distortion matrix is estimated or not
		bool estimateDistortionMatrix() const;

		// Get type of galvanic distortion
		int getTypeOfDistortion() const;

		// Get flag specifing the way of creating roughning matrix
		int geTypeOfRoughningMatrix() const;

		// Get type of the electric field used to calculate response functions
		int getTypeOfElectricField() const;

		// Flag specifing whether type of the electric field of each site is specified indivisually
		bool isTypeOfElectricFieldSetIndivisually() const;

		// Tyep of owner element of observation sites
		int getTypeOfOwnerElement() const;

		// Flag specifing whether the type of owner element of each site is specified indivisually
		bool isTypeOfOwnerElementSetIndivisually() const;

		// Get division number of right-hand sides at solve phase in forward calculation
		int getDivisionNumberOfMultipleRHSInForward() const;

		// Get division number of right-hand sides at solve phase in inversion
		int getDivisionNumberOfMultipleRHSInInversion() const;

		// Get weighting factor of alpha
		double getAlphaWeight(const int iDir) const;

		// Get flag specifing whether the coefficient matrix of the normal equation is positive definite or not
		bool getPositiveDefiniteNormalEqMatrix() const; 

		// Get flag specifing whether output file for paraview is binary or ascii
		bool writeBinaryFormat() const;

		// Get inversion method
		int getInversionMethod() const;

		// Get flag specifing whether observation point is moved to the horizontal center of the element including it
		int getIsObsLocMovedToCenter() const;

		// Get option about treatment of apparent resistivity & phase
		int getApparentResistivityAndPhaseTreatmentOption() const;

		// Get flag specifing whether roughening matrix is outputted
		bool getIsRougheningMatrixOutputted() const;

		// Get type of data space algorithm
		int getTypeOfDataSpaceAlgorithm() const;

		// Get flag specifing whether Lp optimization with difference filter is used
		bool useDifferenceFilter() const;

		// Get degree of Lp optimization
		int getDegreeOfLpOptimization() const;

		// Get lower limit of the difference of log10(rho) for Lp optimization
		double getLowerLimitOfDifflog10RhoForLpOptimization() const;

		// Get upper limit of the difference of log10(rho) for Lp optimization
		double getUpperLimitOfDifflog10RhoForLpOptimization() const;

		// Get maximum iteration number of IRWLS for Lp optimization
		int getMaxIterationIRWLSForLpOptimization() const;

		// Get threshold value for deciding convergence about IRWLS for Lp optimization
		double getThresholdIRWLSForLpOptimization() const;

		// Get directory of out-of-core files for the sensitivitry matrix
		std::string getDirectoryOfOutOfCoreFilesForSensitivityMatrix() const;

		// Get flag specifing whether anisotropic conductivity is considered
		bool isAnisotropyConsidered() const;

		// Get flag specifing whether bottom resistivity is included in roughning
		bool isBottomResistivityIncluded() const;

		// Get resistivity of the bottom of the model
		double getBottomResistivity() const;

		// Get roughning factor at the bottom of the model
		double getRoughningFactorAtBottom() const;

		// Get flag specifing whether small value is added to diagonals of roughening matrix
		bool isSmallValueToRougheningMatrixDiagonals() const;

		// Set small value added to diagonals
		double getSmallValueAddedToDiagonals() const;

		// Get minimum distance between current resistivity and resistivity bounds in common logarithm scale
		double getMinDistanceToBounds() const;

		// Get type of bound constraints
		int getTypeBoundConstraints() const;

		// Get positive real factor of inverse distance weighting
		double getInverseDistanceWeightingFactor() const;

		// Get upper limit of the absolute value of the angle updates for anisotropy
		double getUpperLimitOfAbsAngleUpdatesFoAnisotropy() const;

		// Is the GCV calculation needed?
		bool needsGCVCalculation() const;

		// Get option of GCV calculation
		int getOptionOfGCVCalculation() const;

		// Get pointer to the object of class MeshData
		const MeshData* getPointerOfMeshData() const;

		// Get pointer to the object of class MeshDataBrickElement
		const MeshDataBrickElement* getPointerOfMeshDataBrickElement() const;

		// Get pointer to the object of class MeshDataTetraElement
		const MeshDataTetraElement* getPointerOfMeshDataTetraElement() const;

		// Get pointer to the object of class MeshDataNonConformingHexaElement
		const MeshDataNonConformingHexaElement* getPointerOfMeshDataNonConformingHexaElement() const;

		// Get pointer to the object of class ResistivityBlock
		ResistivityBlock* getPointerOfResistivityBlock() const;

		// Get pointer to the object of class ResistivityBlockIsotropic
		ResistivityBlockIsotropic* getPointerOfResistivityBlockIsotropic() const;

		// Get pointer to the object of class ResistivityBlockAnisotropic
		ResistivityBlockAnisotropic* getPointerOfResistivityBlockAnisotropic() const;

private:
		// Constructer
		AnalysisControl();

		// Destructer
		~AnalysisControl();

		// Copy constructer
		AnalysisControl(const AnalysisControl& rhs){
			std::cerr << "Error : Copy constructer of the class AnalysisControl is not implemented." << std::endl;
			exit(1);
		};

		// Copy assignment operator
		AnalysisControl& operator=(const AnalysisControl& rhs){
			std::cerr << "Error : Assignment operator of the class AnalysisControl is not implemented." << std::endl;
			exit(1);
		}

		// Read analysis control data from "control.dat"
		void inputControlData();

		// Flag specifing whether each parameter has already read from control.dat
		enum controlParameterID{
			BOUNDARY_CONDITION_BOTTOM = 0,
			NUM_THREADS,
			FWD_SOLVER,
			MEM_LIMIT,
			OUTPUT_PARAM,
			NUMBERING_METHOD,
			OUTPUT_OPTION,
			OUTPUT_2D_RESULTS,
			TRADE_OFF_PARAM,
			ITERATION,
			DECREASE_THRESHOLD,
			CONVERGE,
			RETRIAL,
			STEP_LENGTH,
			MESH_TYPE,
			DISTORTION,
			ROUGH_MATRIX,
			ELEC_FIELD,
			DIV_NUM_RHS_FWD,
			DIV_NUM_RHS_INV,
			RESISTIVITY_BOUNDS,
			OFILE_TYPE,
			HOLD_FWD_MEM,
			ALPHA_WEIGHT,
			INV_MAT_POSITIVE_DEFINITE,
			BOTTOM_RESISTIVITY,
			BOTTOM_ROUGHNING_FACTOR,
			INV_METHOD,
			BOUNDS_DIST_THLD,
			IDW,
			SMALL_VALUE,
			MOVE_OBS_LOC,
			OWNER_ELEMENT,
			APP_PHS_OPTION,
			OUTPUT_ROUGH_MATRIX,
			DATA_SPACE_METHOD,
			ANISOTROPY,
			MAX_ANGLE_UPDATE,
			OUTPUT_GCV,
			GCV_OPTION,
			EndOfControlParameterID// This must be written at the end of controlParameterID
		};

		struct ParamsForResistivityBlock {
			bool includeBottomResistivity;// Flag specifing whether bottom resistivity is included in roughning
			double bottomResistivity;// Resistivity of the bottom of the model
			double roughningFactorAtBottom;// Roughning factor at the bottom of the model
			bool addSmallValueToDiagonals;// Flag specifing whether small value is added to diagonals
			double smallValueAddedToDiagonals;// Small value added to bottom
			int typeBoundConstraints;	// Type of bound constraints
			double minDistanceToBounds;// Minimum distance between current resistivity and resistivity bounds in common logarithm scale
			double inverseDistanceWeightingFactor;// Positive real factor of inverse distance weighting
		};

		// Total number of the parameters written in control.dat
		static const int numParamWrittenInControlFile = EndOfControlParameterID;

		// Process ID
		int m_myPE;

		// Total number of processes
		int m_totalPE;

		// Total number of threads
		int m_numThreads;

		// Type of boundary condition at the bottom of the model
		int m_boundaryConditionBottom;

		// Order of finite element
		int m_orderOfFiniteElement;

		// Flag specifing either incore or out-of-core version of PARDISO is used
		int m_modeOfPARDISO;

		// Flag specifing the way of numbering of edges or nodess
		int m_numberingMethod;

		// Parameters to be outputed to the file for visualization
		std::set<int> m_outputParametersForVis;

		// Flag specifing whether the results of 2D forward calculations are outputed
		bool m_isOutput2DResult;

		// Trade-off parameter for the spatial variation of resistivity value
		double m_tradeOffParameterForResistivityValue;

		// Trade-off parameter for the degrees of anistropy
		double m_tradeOffParameterForDegreeOfAnisotropy;

		// Trade-off parameter for the spatial variation of the anistropy direction
		double m_tradeOffParameterForAnisotropyDirection;

		// Trade-off parameter for distortion matrix complexity
		double m_tradeOffParameterForDistortionMatrixComplexity;

		// Trade-off parameter for gains of distortion matrix
		double m_tradeOffParameterForDistortionGain;

		// Trade-off parameter for rotation of distortion matrix
		double m_tradeOffParameterForDistortionRotation;

		// The time the class instanced
		time_t m_startTime;

		// Member variable specifing which backward or forward element is used for calculating EM field
		UseBackwardOrForwardElement m_useBackwardOrForwardElement;

		// Initial iteration number
		int m_iterationNumInit;

		// Maximun iteration number
		int m_iterationNumMax;

		// Current iteration number
		int m_iterationNumCurrent;

		// Threshold value for decreasing
		double m_thresholdValueForDecreasing;

		// Criterion of decrease ratio [%].
		// Iteration If the objective functions and its components decrease by more than criterion,
		// the inversion is considered to be converged. 
		double m_decreaseRatioForConvegence;

		// Current factor of step length damping
		double m_stepLengthDampingFactorCur;

		// Minimum factor of step length damping
		double m_stepLengthDampingFactorMin;

		// Maximum factor of step length damping
		double m_stepLengthDampingFactorMax;

		// If value of objective function decrease specified times in a row, factor of step length damping is increases.
		int m_numOfIterIncreaseStepLength;

		// Factor decreasing step length damping factor
		double m_factorDecreasingStepLength;

		// Factor increasing step length damping factor
		double m_factorIncreasingStepLength;

		// Maximum number of cutbacks.
		int m_numCutbackMax;

		// Hold memory of solver after forward calculation
		bool m_holdMemoryForwardSolver;

		// Pointer to the object of class Forward3DBrickElement0thOrder
		Forward3DBrickElement0thOrder* m_ptrForward3DBrickElement0thOrder;	

		// Pointer to the object of class Forward3DTetraElement0thOrder
		Forward3DTetraElement0thOrder* m_ptrForward3DTetraElement0thOrder;

		// Pointer to the object of class Forward3DNonConformingHexaElement0thOrder
		Forward3DNonConformingHexaElement0thOrder* m_ptrForward3DNonConformingHexaElement0thOrder;

		// Pointer to the object of class Inversion  
		Inversion* m_ptrInversion;

		// Pointer to the object of class ResistivityBlock  
		ResistivityBlock* m_ptrResistivityBlock;

		// Value of objective functional of previous iteration
		double m_objectFunctionalPre;

		// Number of consecutive iteration of which the value of objective functional decrase from the one of previous iteration
		int m_numConsecutiveIterFunctionalDecreasing;

		// Flag specifing whether increment iteration without cutback
		bool m_continueWithoutCutback; 

		// Maximum value of the memory used by out-of-core mode of PARDISO 
		int m_maxMemoryPARDISO;

		// Type of mesh
		int m_typeOfMesh;

		// Flag specifing the way of creating roughning matrix
		int m_typeOfRoughningMatrix;

		// Type of the electric field used to calculated response functions
		int m_typeOfElectricField;

		// Flag specifing whether type of the electric field of each site is specified indivisually
		bool m_isTypeOfElectricFieldSetIndivisually;

		// Type of owner element of observation sites
		int m_typeOfOwnerElement;

		// Flag specifing whether the owner element (upper or lower) of each site is specified indivisually
		bool m_isTypeOfOwnerElementSetIndivisually;

		// Division number of right-hand sides at solve phase in forward calculation
		int m_divisionNumberOfMultipleRHSInForward;

		// Division number of right-hand sides at solve phase in inversion
		int m_divisionNumberOfMultipleRHSInInversion;

		// Flag specifing whether output file for paraview is binary or ascii
		bool m_binaryOutput;

		// Type of galvanic distortion
		int m_typeOfDistortion;

		// Weighting factor of alpha
		double m_alphaWeight[3];

		// Flag specifing whether the coefficient matrix of the normal equation is positive definite or not
		bool m_positiveDefiniteNormalEqMatrix; 

		// Inversion method
		int m_inversionMethod;

		// Flag specifing whether the observation point is moved to the center of the element
		bool m_isObsLocMovedToCenter;

		// Option about treatment of apparent resistivity & phase
		int m_apparentResistivityAndPhaseTreatmentOption;

		// Flag specifing whether roughening matrix is outputed
		bool m_isRougheningMatrixOutputted;

		// Flag specifing whether anisotropic conductivity is considered
		bool m_isAnisotropyConsidered;

		// Flag specifing whether bottom resistivity is included in roughning
		bool m_isBottomResistivityincluded;

		// Resistivity of the bottom of the model
		double m_bottomResistivity;

		// Roughning factor at the bottom of the model
		double m_roughningFactorAtBottom;

		// Flag specifing whether small value is added to diagonals of roughening matrix
		bool m_isSmallValueAddedToRougheningMatrixDiagonals;

		// Small value added to diagonals of roughening matrix
		double m_smallValueAddedToRougheningMatrixRougheningMatrixDiagonals;

		// Type of bound constraints
		int m_typeBoundConstraints;

		// Minimum distance between current resistivity and resistivity bounds in common logarithm scale
		double m_minDistanceToBounds;

		// Positive real factor of inverse distance weighting
		double m_inverseDistanceWeightingFactor;

		// Upper limit of the absolute value of the angle updates for anisotropy
		double m_upperLimitOfAbsAngleUpdatesFoAnisotropy;

		// Is the GCV calculation needed?
		bool m_needsGCVCalculation;

		// Option of GCV calculation
		int m_optionOfGCVCalculation;

		// Type of data space algorithm
		int m_typeOfDataSpaceAlgorithm;

		// Flag specifing whether Lp optimization with difference filter is used
		bool m_useDifferenceFilter;

		// Degree of Lp optimization
		int m_degreeOfLpOptimization;

		// Lower limit of the difference of log10(rho) for Lp optimization
		double m_lowerLimitOfDifflog10RhoForLpOptimization;

		// Upper limit of the difference of log10(rho) for Lp optimization
		double m_upperLimitOfDifflog10RhoForLpOptimization;

		// Maximum iteration number of IRWLS for Lp optimization
		int m_maxIterationIRWLSForLpOptimization;

		// Threshold value for deciding convergence about IRWLS for Lp optimization
		double m_thresholdIRWLSForLpOptimization;

		// Directory of out-of-core files for the sensitivitry matrix
		std::string m_directoryOfOutOfCoreFilesForSensitivityMatrix;

		// Calculate forward computation
		void calcForwardComputation( const int iter );

		// Adjust factor of step length damping and output convergence data to cnv file
		AnalysisControl::ConvergenceBehaviors adjustStepLengthDampingFactor( const int iterCur, const int iCutbackCur );

		// Check convergece of the inversion
		bool checkConvergence( const double objectFunctionalCur );

		// Calculate regularization terms and write them to the cnv file for isotropic resistivity
		// @note: Objective function is returned
		const double calcRegularizationTermsAndWriteThemToCnvFileIsotropic(const int numDataTotal, const double dataMisfit) const;

		// Calculate regularization terms and write them to the cnv file for anisotropic resistivity
		// @note: Objective function is returned
		const double calcRegularizationTermsAndWriteThemToCnvFileAnisotropic(const int numDataTotal, const double dataMisfit, const double dataMisfitPredicted) const;

		// Return flag specifing whether sensitivity is calculated or not
		bool doesCalculateSensitivity( const int iter ) const;

		// Get pointer to the object of class Forward3D
		Forward3D* getPointerOfForward3D() const;

};

#endif
