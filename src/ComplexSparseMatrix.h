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
#ifndef DBLDEF_COMPLEX_SPARSE_MATRIX
#define DBLDEF_COMPLEX_SPARSE_MATRIX

#include <complex>
//#include <set>
#include <map>

class ComplexSparseMatrix{

public:

	//Default Constructer
	explicit ComplexSparseMatrix();

	//Constructer
	explicit ComplexSparseMatrix( const int nrows, const int ncols, const int nrhs = 1 );

	//Destructer
	virtual ~ComplexSparseMatrix();

	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Set matrix structure ( locations of non-zero components ) by triplet format
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const std::complex<double>& val );

	//Convert matrix from triplet format to CRS format
	void convertToCRSFormat();

	//Add non-zero value to matrix
	virtual void addNonZeroValues( const int row, const int col, const std::complex<double>& val );

	//Check input data and get element number of the array containing non-zero values
	virtual int checkAndGetLocationNonZeroValue( const int row, const int col );

	//Add non-zero value to matrix by specifing element number of the array directly 
	virtual void addNonZeroValuesWithoutSearchingLocation( const int loc, const std::complex<double>& val );

	//Zero clear non-zero values of matrix stored by CSR format
	void zeroClearNonZeroValues();

	//Add non-zero value to the right hand side vector
	void addRightHandSideVector( const int row, const std::complex<double>& val, const int irhs = 0 );

	//Zero clear non-zero values of the right hand side vector
	void zeroClearRightHandSideVector();

	//Initialize matrix and right-hand side vectors
	virtual void initializeMatrixAndRhsVectors( const int nrows, const int ncols, const int nrhs );

	// Get number of rows
	int getNumRows() const;

	// Get number of columns
	int getNumColumns() const;

	//Get total number of right hand side vector 
	int getNumRightHandSideVectors() const;

	//Return whether matrix has already been converted to CRS format
	bool hasConvertedToCRSFormat() const;

	//Reallocate memory for right hand side vectors
	void reallocateMemoryForRightHandSideVectors( const int nrhs );

	////Set number of right hand side vectors
	//void setNumRightHandSideVectors( const int nrhs );

	//Release memory
	virtual void releaseMemory();

	////Multiply matrix by inputed vector and add calculated vector to another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndAddResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const;

	////Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndSubtractResult( const std::complex<double>* vecIn, const int* convertArray, std::complex<double>* vecOut ) const;

	//Copy right-hand side vector to another vector
	void copyRhsVector( std::complex<double>* vecOut ) const;

	//Copy specified components of right-hand side vector to another vector
	void copyRhsVector( const int numCompsCopied, const int* const compsCopied, std::complex<double>* vecOut ) const;

	//Debug write matrix componets
	void debugWriteMatrix() const;

	//Debug write componets of right hand side vector 
	void debugWriteRightHandSide() const;

	//Debug write non-zero componets of right hand side vector 
	void debugWriteNonZeroRightHandSide() const;

protected:
	//Delete the matrix of triplet ( Coordinate ) format
	void deleteTripletMatrix();

	//Total number of rows
	int m_numRows;

	//Total number of columns
	int m_numColumns;

	//Total number of non-zero elements in the matrix
	int m_numNonZeros;

	//Total number of right hand side vectors
	int m_numRightHandSideVectors;

	//Flag indicating whether matrix has already been converted to CRS format
	bool m_hasConvertedToCRSFormat;

	//Matrix stored by triplet ( Coordinate ) format
	//std::set<int>* m_matrixTripletFormat;
	std::map< int, std::complex<double> >* m_matrixTripletFormat;

	//Row indices of the compressed row storage format
	long long int* m_rowIndex;

	//Columns in which non-zero compnents exist
	long long int* m_columns;

	//Values of non-zero compnents
	std::complex<double>* m_values;

	//Right hand side vector
	std::complex<double>* m_rightHandSideVector;
	
private:
	//Copy constructer
	ComplexSparseMatrix(const ComplexSparseMatrix &matrix );

	// Assignment operator
	ComplexSparseMatrix& operator=(const ComplexSparseMatrix& rhs);

};

#endif
