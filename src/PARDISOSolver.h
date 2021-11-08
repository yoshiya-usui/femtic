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
#ifndef DBLDEF_PARDISO_SOLVER
#define DBLDEF_PARDISO_SOLVER

#include <string>

// Class of PARDISO solver
class PARDISOSolver{

public:
	enum solutionStage{
		MEMORY_RELEASED = 0,
		INITIALIZED,
		ANALYZED,
		FACTORIZED,
		SOLVED,
	};

	const static int INCORE_MODE = 0;
	const static int SELECT_MODE_AUTOMATICALLY = 1;
	const static int OUT_OF_CORE_MODE = 2;
		
	const static int REAL_AND_STRUCTURALLY_SYMMETRIC = 1;
	const static int REAL_AND_SYMMETRIC_POSITIVE_DEFINITE = 2;
	const static int REAL_AND_SYMMETRIC_INDEFINITE = -2;
	const static int COMPLEX_AND_STRUCTURALLY_SYMMETRIC =3;
	const static int COMPLEX_AND_HERMITIAN_POSITIVE_DEFINITE = 4;
	const static int COMPLEX_AND_HERMITIAN_INDEFINITE = -4;
	const static int COMPLEX_AND_SYMMETRIC_MATRIX = 6;
	const static int REAL_AND_UNSYMMETRIC_MATRIX = 11;
	const static int COMPLEX_AND_UNSYMMETRIC_MATRIX = 13;

	// Default constructer
	explicit PARDISOSolver();

	// Constructer
	explicit PARDISOSolver(const long long int matrixType);

	// Destructer
	virtual ~PARDISOSolver();
	
	// Initialize PARDISOr solver
	void initialize( const std::string& oocHeaderName, const long long int imode, const long long int type );

	// Analysis phase of PARDISO solver
	void analysis( long long int nEq, long long int* rowIndex, long long int* columns );

	// Release memory of PARDISO solver
	void releaseMemory();

	// Get memory required by PARDISO solver
	void writeMemoryRequired() const;

	// Get stage of PARDISO solver
	int getSolutionStage() const;

	// Set stage of PARDISO solver
	void setSolutionStage( const int stage );

protected:
	// Internal solver memory pointer of PARDISO Solver
	void* m_pt[64];

	// Maximum number of factors with identical nonzero sparsity structure
	long long int m_maxfct;

	// The actual matrix to be factorized at the solution phase
	long long int m_mnum;

	// Type of matrix
	long long int m_mtype;

	// Number of equations
	long long int m_numEquations;

	// Parameters which control PARDISO
	long long int m_iparm[64];

	// Message level of PARDISO
	long long int m_msglvl;

	// Total peak memory [KByte] that needs during the analysis and symbolic factorization phase
	double m_peakMemorySymbolicFactorization;

	// The permanent memory [KByte] that needs from the analysis and symbolic factorization phase to the factorization and solve phases.
	double m_permanetMemorySymbolicFactorization;

	// The total memory consumed by in-core PARDISO for internal float point arrays
	double m_memoryForNumericalFactorizationIncore;

	// The total memory consumed by out-core PARDISO for internal float point arrays
	double m_memoryForNumericalFactorizationOutcore;

	// Solution stage of PARDISO solver
	int m_solutionStage;

	// Output error messages
	void outputErrorMessages( const int ier ) const;

private:
	// Copy constructer
	PARDISOSolver(const PARDISOSolver& rhs);

	// Copy assignment operator
	PARDISOSolver& operator=(const PARDISOSolver& rhs);

};

#endif
