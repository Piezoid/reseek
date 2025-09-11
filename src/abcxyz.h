#pragma once

#include "myutils.h"
#include "arrays.h"

double dot(const double a[3], const double b[3]);
void transform(const double t[3], const double u[3][3],
  const double x[3], double x_transformed[3]);

static const double PI = 3.1415926535;
static const int X = 0;
static const int Y = 1;
static const int Z = 2;

static const int A = 0;
static const int B = 1;
static const int C = 2;

const uint AL = 12;
const uint BL = 14;
const uint CL = 8;
const uint PALMCOREM = 150;

extern const uint MotifLVec[3];

void LogVec(const string &Msg, const vector<double> &v);
void LogMx(const string &Msg, const SquareMatrix<double> &Mx);
double GetMxDeterminant(const SquareMatrix<double> &Mx);
void InvertMx(const SquareMatrix<double> &Mx,
  SquareMatrix<double> &InvMx);
void MulMxVec(const SquareMatrix<double> &Mx,
  const vector<double> &Vec,
  vector<double> &Result);
void MulMx(
  const SquareMatrix<double> &A,
  const SquareMatrix<double> &B,
  SquareMatrix<double> &Prod);
void GetBasisR(const SquareMatrix<double> &Basis,
  SquareMatrix<double> &R);
void LogMx(const string &Msg, const SquareMatrix<double> &Mx);
void RotateMx(const SquareMatrix<double> &Mx,
  uint Axis, double Theta, SquareMatrix<double> &RotatedMx);
void CrossProduct(
  const vector<double> &a,
  const vector<double> &b,
  vector<double> &Prod);
void GetTriForm(
  SquareMatrix<double> &MotifCoords,
  vector<double> &t,
  SquareMatrix<double> &R);

void GetTriangleCentroid(const SquareMatrix<double> &MotifCoords,
  vector<double> &CentroidCoords);
void GetTriangleBasis(const SquareMatrix<double> &MotifCoords,
  vector<double> &CentroidCoords, SquareMatrix<double> &Basis);

static inline uint GetOtherAxis_i(uint Axis) { return (Axis+1)%3;  }
static inline uint GetOtherAxis_j(uint Axis) { return (Axis+2)%3;  }

static inline void Resize3(vector<double> &v) { v.resize(3); }
static inline void Resize3(vector<float> &v) { v.resize(3); }

static inline SquareMatrix<double> Allocate3x3()
	{
	return SquareMatrix<double>::Allocate(3);
	}

static inline void MulVecScalar(const vector<double> &v,
  double s, vector<double> &r)
	{
	asserta(SIZE(v) == 3);
	r.resize(3);
	r[0] = v[0]*s;
	r[1] = v[1]*s;
	r[2] = v[2]*s;
	}

static inline double GetDist2(
  const vector<double> &Pt1,
  const vector<double> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	double dx = Pt1[X] - Pt2[X];
	double dy = Pt1[Y] - Pt2[Y];
	double dz = Pt1[Z] - Pt2[Z];
	double d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

static inline double GetDist(
  const vector<double> &Pt1,
  const vector<double> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	double dx = Pt1[X] - Pt2[X];
	double dy = Pt1[Y] - Pt2[Y];
	double dz = Pt1[Z] - Pt2[Z];
	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

static inline double GetDist3D(
  double x1, double y1, double z1,
  double x2, double y2, double z2)
	{
	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;
	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

static inline float GetDist3D(
  float x1, float y1, float z1,
  float x2, float y2, float z2)
	{
	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;
	float d2 = dx*dx + dy*dy + dz*dz;
	float d = sqrt(d2);
	return d;
	}

static inline double GetDist2_Mxij(const SquareMatrix<double> &Mx,
  uint i, uint j)
	{
	double xi = Mx[i][X];
	double yi = Mx[i][Y];
	double zi = Mx[i][Z];

	double xj = Mx[j][X];
	double yj = Mx[j][Y];
	double zj = Mx[j][Z];

	double dx = xi - xj;
	double dy = yi - yj;
	double dz = zi - zj;

	double d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

static inline double GetDist_Mxij(const SquareMatrix<double> &Mx,
  uint i, uint j)
	{
	double xi = Mx[i][X];
	double yi = Mx[i][Y];
	double zi = Mx[i][Z];

	double xj = Mx[j][X];
	double yj = Mx[j][Y];
	double zj = Mx[j][Z];

	double dx = xi - xj;
	double dy = yi - yj;
	double dz = zi - zj;

	double d2 = dx*dx + dy*dy + dz*dz;
	double d = sqrt(d2);
	return d;
	}

static inline double GetMod_xyz(double x, double y, double z)
	{
	double d2 = x*x + y*y + z*z;
	double Mod = sqrt(d2);
	return Mod;
	}

static inline double GetMod_Vec(const vector<double> &v)
	{
	double x = v[X];
	double y = v[Y];
	double z = v[Z];
	double Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

static inline void NormalizeVec(vector<double> &v)
	{
	double Mod = GetMod_Vec(v);
	assert(Mod > 0);
	v[X] /= Mod;
	v[Y] /= Mod;
	v[Z] /= Mod;
	assert(feq(GetMod_Vec(v), 1));
	}

static inline double GetMod_Mxi(const SquareMatrix<double> &Mx,
  uint i)
	{
	double x = Mx[i][X];
	double y = Mx[i][Y];
	double z = Mx[i][Z];
	double Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

static inline double GetTheta3D(
  double xi, double yi, double zi,
  double xj, double yj, double zj)
	{
	double DotProd = xi*xj + yi*yj + zi*zj;
	double Modi = GetMod_xyz(xi, yi, zi);
	double Modj = GetMod_xyz(xj, yj, zj);
	if (fabs(Modi*Modj) < 1e-6)
		return 0;
	double cos_theta = DotProd/(Modi*Modj);
	asserta(cos_theta >= -1.02 && cos_theta <= 1.02);
	if (cos_theta < -1)
		cos_theta = -1;
	else if (cos_theta > 1)
		cos_theta = 1;
	double theta = acos(cos_theta);
	return theta;
	}

static inline double GetTheta3D(
  const vector<double> &vi,
  const vector<double> &vj)
	{
	asserta(SIZE(vi) == 3);
	asserta(SIZE(vj) == 3);
	return GetTheta3D(vi[0], vi[1], vi[2], vj[0], vj[1], vj[2]);
	}

static inline void Sub_Vecs(const vector<double> &vi,
  const vector<double> &vj, vector<double> &Diff)
	{
	Resize3(Diff);
	Diff[X] = vi[X] - vj[X];
	Diff[Y] = vi[Y] - vj[Y];
	Diff[Z] = vi[Z] - vj[Z];
	}

static inline void Add_Vecs(const vector<double> &vi,
  const vector<double> &vj, vector<double> &Sum)
	{
	Resize3(Sum);
	Sum[X] = vi[X] + vj[X];
	Sum[Y] = vi[Y] + vj[Y];
	Sum[Z] = vi[Z] + vj[Z];
	}

static inline double GetTheta_Vecs(const vector<double> &vi,
  const vector<double> &vj)
	{
	double xi = vi[X];
	double yi = vi[Y];
	double zi = vi[Z];

	double xj = vj[X];
	double yj = vj[Y];
	double zj = vj[Z];

	double theta = GetTheta3D(xi, yi, zi, xj, yj, zj);
	return theta;
	}

static inline double GetTheta_Mxij(const SquareMatrix<double> &Mx,
  uint i, uint j)
	{
	vector<double> vi(3), vj(3);
	for (uint k = 0; k < 3; ++k)
		{
		vi[k] = Mx[i][k];
		vj[k] = Mx[j][k];
		}
	double theta = GetTheta_Vecs(vi, vj);
	return theta;
	}

static inline double degrees(double Radians) { return Radians*180.0/PI; }

static inline double degrees_0_to_360(double Radians)
	{
	double Deg = Radians*180.0/PI;
	Deg = fmod(Deg, 360);
	if (Deg < 0)
		Deg += 360;
	assert(Deg >= 0 && Deg < 360);
	return Deg;
	}

void GetIdentityMx(SquareMatrix<double> &Mx);

void XFormPt(
  const vector<double> &Pt,
  const vector<double> &t,
  const SquareMatrix<double> &R,
  vector<double> &XPt);

void XFormMx(
  const SquareMatrix<double> &Mx,
  const vector<double> &t,
  const SquareMatrix<double> &R,
  SquareMatrix<double> &XMx);

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const SquareMatrix<float> &R,
  vector<float> &XPt);

void XFormMx(
  const SquareMatrix<float> &Mx,
  const vector<float> &t,
  const SquareMatrix<float> &R,
  SquareMatrix<float> &XMx);

#if DEBUG
static void AssertCanonicalUnitBasis(const SquareMatrix<double> &Basis)
	{
	assert(Basis.Rows() == 3);
	assert(Basis.Cols() == 3);

	assert(feq(Basis[0][X], 1));
	assert(feq(Basis[0][Y], 0));
	assert(feq(Basis[0][Z], 0));

	assert(feq(Basis[1][X], 0));
	assert(feq(Basis[1][Y], 1));
	assert(feq(Basis[1][Z], 0));

	assert(feq(Basis[2][X], 0));
	assert(feq(Basis[2][Y], 0));
	assert(feq(Basis[2][Z], 1));
	}
#else
#define AssertCanonicalUnitBasis(x)	0
#endif // DEBUG

static void AssertUnitBasisA(const SquareMatrix<double> &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		double Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		double theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}

#if DEBUG
static void AssertUnitBasis(const SquareMatrix<double> &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		double Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		double theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}
#else
#define AssertUnitBasis(x)	0
#endif // DEBUG

#if DEBUG
static void AssertSameLengths(const SquareMatrix<double> &Mx1,
  const SquareMatrix<double> &Mx2)
	{
	for (uint k = 0; k < 3; ++k)
		{
		double Modk1 = GetMod_Mxi(Mx1, k);
		double Modk2 = GetMod_Mxi(Mx2, k);
		assert(feq(Modk1, Modk2));
		}
	}
#else
#define AssertSameLengths(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertSameAngles(const SquareMatrix<double> &Mx1,
  const SquareMatrix<double> &Mx2)
	{
	for (uint i = 0; i < 3; ++i)
		{
		for (uint j = i + 1; j < 3; ++j)
			{
			double Theta1 = GetTheta_Mxij(Mx1, i, j);
			double Theta2 = GetTheta_Mxij(Mx2, i, j);
			asserta(feq(Theta1, Theta2));
			}
		}
	}
#else
#define AssertSameAngles(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertMx3D(const SquareMatrix<double> &Mx)
	{
	asserta(Mx.Rows() == 3);
	asserta(Mx.Cols() == 3);
	}
#else
#define AssertMx3D(x)	0
#endif
