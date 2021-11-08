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
#ifndef DBLDEF_UTIL
#define DBLDEF_UTIL

#include <iostream>
#include <string>
#include <algorithm>
#include <complex>

#include "CommonParameters.h"

// Flag specifing comparative operators
enum comparativeOperators{
	EQUAL=0,
	LARGE_LEFT_HAND_SIDE,
	SMALL_LEFT_HAND_SIDE,
};

// Sort elements by its key value with quick sort
void quickSort( const int numOfIDs, int* ids, const double* values );

// Sort elements by values of its three keys with quick sort
void quickSortThreeKeys( const int numOfIDs, int* ids,
	const double* firstKeyValues, const double* secondKeyValues, const double* thirdKeyValues );

// Compare values by its three keys
int compareValueByThreeKeys( const int lhsID, const int rhsID, const double* firstKeyValues, const double* secondKeyValues, const double* thirdKeyValues );

// Calculate matrix product for 2 x 2 double matrix
void calcProductFor2x2DoubleMatrix( const CommonParameters::DoubleMatrix2x2& matInA, const CommonParameters::DoubleMatrix2x2& matInB, CommonParameters::DoubleMatrix2x2& matOut );

// Calculate 3D vectors
CommonParameters::Vector3D calcVector3D( const CommonParameters::locationXYZ& startCoords, const CommonParameters::locationXYZ& endCoords );

// Calculate outer product of 3D vectors
CommonParameters::Vector3D calcOuterProduct( const CommonParameters::Vector3D& vec1, const CommonParameters::Vector3D& vec2 );

// Calculate inner product of 3D vectors
double calcInnerProduct( const CommonParameters::Vector3D& vec1, const CommonParameters::Vector3D& vec2 );

// Calculate impedance tensor component from apparent resisitivity and phase
void calcImpedanceTensorComponentFromApparentResistivityAndPhase( const double freq, const double appRes, const double appResError,
																 const double phase, const double phaseError,
																 std::complex<double>& Z, 	CommonParameters::DoubleComplexValues& ZError );

#ifdef _ANISOTOROPY
void rotateTensor( double tensorRotated[3][3], const double rotationTensor[3][3] );
#endif

#endif
