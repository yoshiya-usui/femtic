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
#ifndef DBLDEF_COMMON_PARAMETERS
#define DBLDEF_COMMON_PARAMETERS

#include <stdio.h>
#include <math.h>

namespace CommonParameters{

static const int EX_POLARIZATION = 0;
static const int EY_POLARIZATION = 1;

static const int TM_MODE = 0;
static const int TE_MODE = 1;

struct locationXY{
	double X;
	double Y;
};

struct locationYZ{
	double Y;
	double Z;
};

struct locationZX{
	double Z;
	double X;
};

struct locationXYZ{
	double X;
	double Y;
	double Z;
};

struct locationDipole{
	CommonParameters::locationXY startPoint;
	CommonParameters::locationXY endPoint;
};

struct DoubleComplexValues{
	double realPart;
	double imagPart;
};

struct InitComplexValues{
	int realPart;
	int imagPart;
};

struct DoubleMatrix2x2{
	double comp11;
	double comp12;
	double comp21;
	double comp22;
};

struct DoubleMatrix3x3{
	double comp11;
	double comp12;
	double comp13;
	double comp21;
	double comp22;
	double comp23;
	double comp31;
	double comp32;
	double comp33;
};

struct Vector3D{
	double X;
	double Y;
	double Z;
};

struct AreaCoords{
	double coord0;
	double coord1;
	double coord2;
};

struct VolumeCoords{
	double coord0;
	double coord1;
	double coord2;
	double coord3;
};

struct CoordPair{
	double first;
	double second;
};

// Circular constant
static const double PI = 3.14159265358979;

// Factor converting values from radians to degrees
static const double rad2deg = 180.0 / 3.14159265358979;

// Factor converting values from degrees to radians
static const double deg2rad = 3.14159265358979 / 180.0;

// Magnetic permeability
static const double mu = 12.566370614e-07;

// Value of source electric field
static const double sourceValueElectric = 1000.0;

// Electric permittivity
//static const double epsilon = 8.854187817e-12;
static const double epsilon = 0.0;

// Very small value
static const double EPS = 1e-12;

// Abscissas of one point Gauss quadrature
static const double abscissas1Point[1] = { 0.0 };

// Abscissas of two point Gauss quadrature
static const double abscissas2Point[2] = { -1.0/sqrt(3.0), 1.0/sqrt(3.0) };

// Abscissas of three point Gauss quadrature
static const double abscissas3Point[3] = { -0.774596669241483, 0.0, 0.774596669241483 };

// Weights of one point Gauss quadrature
static const double weights1Point[1] = { 2.0 };

// Weights of two point Gauss quadrature
static const double weights2Point[2] = { 1.0, 1.0 };

// Weights of three point Gauss quadrature
static const double weights3Point[3] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };

// Factor converting value from kilo-meter to meter
static const double convKilometerToMeter = 1000.0;

static char programName[]="femtic";

// [MajorVersion#].[MinorVersion#].[Revision#]
// x.x.xa -> alpha version
// x.x.xb -> beta version
static char versionID[]="4.2.5";

}

#endif
