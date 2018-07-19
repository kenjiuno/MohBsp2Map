
#include "StdAfx.h"
#include <algorithm>
#include <math.h>
#include "Texc.h"

using namespace std;

namespace
{
	// 
	inline double calcTexc_Length(
		double x, double y
		)
	{
		double z = sqrt(x * x + y * y);
		return z;
	}
	// 
	inline double calcDEGFrom0To360(double a) { while ((int)floor(a + 0.5) < 0) a += 360; while ((int)floor(a + 0.5) >= 360) a -= 360; return a; }
	// 
	inline double calcDEGFromRAD(double a) { return a / 3.14159265358979 * 180.0; }
	// 
	double calcTexc_AngleXY(
		double x, double y
		)
	{
		if (fabs(x) < 1E-5) {
			if (y < 0)
				//270‹
				return -3.1415926535897932384626433832795 * 0.5;
			else
				// 90‹
				return +3.1415926535897932384626433832795 * 0.5;
		} else {
			if (x > 0)
				return atan(y / x);
			else
				return atan(y / x) + 3.1415926535897932384626433832795;
		}
	}
	// 
	inline double calcTexc_AngleXYd(
		double x, double y
		)
	{
		double a = calcDEGFrom0To360(calcDEGFromRAD(calcTexc_AngleXY(x, y)));
		return a;
	}
	// 
	inline double calcRADFromDEG(double a) { return a / 180.0 * 3.14159265358979; }
	// 
	inline void calcTexc_rotateXYd(double fX, double fY, double &fXr, double &fYr, double f)
	{
		f = calcRADFromDEG(f);
		double eM11 = +cos(f);
		double eM12 = +sin(f);
		double eM21 = -eM12;
		double eM22 = +eM11;
		fXr = fX * eM11 + fY * eM21;
		fYr = fX * eM12 + fY * eM22;
	}
	// 
	inline void calcTexc_rotateTriPos(
		double fAx0, double fAy0,
		double fAx1, double fAy1,
		double fAx2, double fAy2,
		UINT cxTex, UINT cyTex,
		double fTexsX, double fTexsY, double fTexAngle,
		double &fEx, double &fEy
	)
	{
		if (fabs(fTexsX) < 1E-5 || fabs(fTexsY) < 1E-5) return;

		double fAxMin = __min(fAx0, fAx1); fAxMin = __min(fAxMin, fAx2);
		double fAyMin = __min(fAy0, fAy1); fAyMin = __min(fAyMin, fAy2);
		double fAxOff = (int)(fAxMin / cxTex); fAxOff = 0;
		double fAyOff = (int)(fAyMin / cyTex); fAyOff = 0;

		double fSx0 = (double)cxTex / (double)cyTex;
		double fSx1 = fabs(fTexsX) * 2;
		double fSy1 = fabs(fTexsY) * 2;
		double fTx0 = ((fAx0 - fAxOff) / cxTex);
		double fTy0 =-((fAy0 - fAyOff) / cyTex);
		double fTx1 = ((fAx1 - fAxOff) / cxTex);
		double fTy1 =-((fAy1 - fAyOff) / cyTex);
		double fTx2 = ((fAx2 - fAxOff) / cxTex);
		double fTy2 =-((fAy2 - fAyOff) / cyTex);
		double fTxMin = __min(fTx0, fTx1); fTxMin = __min(fTxMin, fTx2);
		double fTyMin = __min(fTy0, fTy1); fTyMin = __min(fTyMin, fTy2);
		double fTxMax = __max(fTx0, fTx1); fTxMax = __max(fTxMax, fTx2);
		double fTyMax = __max(fTy0, fTy1); fTyMax = __max(fTyMax, fTy2);
		double fTyLen = ceil(fTyMax) - floor(fTyMin);
		while (fTyMin < -fabs(fTexsY)) { fTyMin += fSy1, fTy0 += fSy1, fTy1 += fSy1, fTy2 += fSy1; }
		while (fTxMin < -fabs(fTexsX)) { fTxMin += fSx1, fTx0 += fSx1, fTx1 += fSx1, fTx2 += fSx1; }
		while (fTyMin > +fabs(fTexsY)) { fTyMin -= fSy1, fTy0 -= fSy1, fTy1 -= fSy1, fTy2 -= fSy1; }
		while (fTxMin > +fabs(fTexsX)) { fTxMin -= fSx1, fTx0 -= fSx1, fTx1 -= fSx1, fTx2 -= fSx1; }
		fEx = fTx0;
		fEy = fTy0 /*- fTyLen*/;
		fEx *= fSx0;
		calcTexc_rotateXYd(fEx, fEy, fEx, fEy, fTexAngle);
		fEx /= fSx0;
		fEx /= fTexsX;
		fEy /= fTexsY;
		fEx += 0;
	//	fEy += fTyLen;
	}
};

bool MB2MM::Texc::Decode()
{
	ASSERT(cxTex != 0 && cyTex != 0);

	if (fMirror) {
		swap(vv[1][0], vv[2][0]);
		swap(tv[1][0], tv[2][0]);
	}
	double fTs[][3][2] = {
		tv[0][0], tv[0][1], tv[1][0], tv[1][1], tv[2][0], tv[2][1],
		tv[0][0], tv[0][1], tv[2][0], tv[2][1], tv[1][0], tv[1][1],
		tv[1][0], tv[1][1], tv[2][0], tv[2][1], tv[0][0], tv[0][1],
		tv[1][0], tv[1][1], tv[0][0], tv[0][1], tv[2][0], tv[2][1],
		tv[2][0], tv[2][1], tv[0][0], tv[0][1], tv[1][0], tv[1][1],
		tv[2][0], tv[2][1], tv[1][0], tv[1][1], tv[0][0], tv[0][1],
	};
	double fAs[][3][2] = {
		vv[0][0], vv[0][1], vv[1][0], vv[1][1], vv[2][0], vv[2][1],
		vv[0][0], vv[0][1], vv[2][0], vv[2][1], vv[1][0], vv[1][1],
		vv[1][0], vv[1][1], vv[2][0], vv[2][1], vv[0][0], vv[0][1],
		vv[1][0], vv[1][1], vv[0][0], vv[0][1], vv[2][0], vv[2][1],
		vv[2][0], vv[2][1], vv[0][0], vv[0][1], vv[1][0], vv[1][1],
		vv[2][0], vv[2][1], vv[1][0], vv[1][1], vv[0][0], vv[0][1],
	};
	for (int iXc = 0; iXc < 6; iXc++) {
		tv[0][0] = fTs[iXc][0][0], tv[0][1] = fTs[iXc][0][1];
		tv[1][0] = fTs[iXc][1][0], tv[1][1] = fTs[iXc][1][1];
		tv[2][0] = fTs[iXc][2][0], tv[2][1] = fTs[iXc][2][1];
		vv[0][0] = fAs[iXc][0][0], vv[0][1] = fAs[iXc][0][1];
		vv[1][0] = fAs[iXc][1][0], vv[1][1] = fAs[iXc][1][1];
		vv[2][0] = fAs[iXc][2][0], vv[2][1] = fAs[iXc][2][1];

		if (!(tv[0][1] <= tv[1][1] && tv[1][1] <= tv[2][1]))
			continue;
		tv[1][0] -= tv[0][0], tv[1][1] -= tv[0][1], tv[2][0] -= tv[0][0], tv[2][1] -= tv[0][1], tv[0][0] = tv[0][1] = 0;
		vv[1][0] -= vv[0][0], vv[1][1] -= vv[0][1], vv[2][0] -= vv[0][0], vv[2][1] -= vv[0][1], vv[0][0] = vv[0][1] = 0;
		double fTxPoint = tv[1][1] / tv[2][1];
		double fTxDeltaTx = tv[2][0] * fTxPoint - tv[1][0];
		double fTxDeltaTy = 0;
		double fTxDeltaAx = (fTxDeltaTx < 0) ? (vv[1][0] - vv[2][0] * fTxPoint) : (vv[2][0] * fTxPoint - vv[1][0]);
		double fTxDeltaAy = (fTxDeltaTx < 0) ? (vv[1][1] - vv[2][1] * fTxPoint) : (vv[2][1] * fTxPoint - vv[1][1]);
		double fTxLengthA = calcTexc_Length(fTxDeltaAx, fTxDeltaAy);
		double fTxLengthT = calcTexc_Length(fTxDeltaTx * cxTex, 0);
		tsx = fTxLengthA / fTxLengthT;
		BOOL bMirrorX = (iXc & 1) ? TRUE : FALSE;
		double fTxAngleA = calcTexc_AngleXYd(fTxDeltaAx, fTxDeltaAy);

		for (int iYc = 0; iYc < 6; iYc++) {
			tv[0][0] = fTs[iYc][0][0], tv[0][1] = fTs[iYc][0][1];
			tv[1][0] = fTs[iYc][1][0], tv[1][1] = fTs[iYc][1][1];
			tv[2][0] = fTs[iYc][2][0], tv[2][1] = fTs[iYc][2][1];
			vv[0][0] = fAs[iYc][0][0], vv[0][1] = fAs[iYc][0][1];
			vv[1][0] = fAs[iYc][1][0], vv[1][1] = fAs[iYc][1][1];
			vv[2][0] = fAs[iYc][2][0], vv[2][1] = fAs[iYc][2][1];

			if (!(tv[0][0] <= tv[1][0] && tv[1][0] <= tv[2][0]))
				continue;
			tv[1][0] -= tv[0][0], tv[1][1] -= tv[0][1], tv[2][0] -= tv[0][0], tv[2][1] -= tv[0][1], tv[0][0] = tv[0][1] = 0;
			vv[1][0] -= vv[0][0], vv[1][1] -= vv[0][1], vv[2][0] -= vv[0][0], vv[2][1] -= vv[0][1], vv[0][0] = vv[0][1] = 0;
			double fTyPoint = tv[1][0] / tv[2][0];
			double fTyDeltaTx = 0;
			double fTyDeltaTy = tv[2][1] * fTyPoint - tv[1][1];
			double fTyDeltaAx = (fTyDeltaTy < 0) ? (vv[2][0] * fTyPoint - vv[1][0]) : (vv[1][0] - vv[2][0] * fTyPoint);
			double fTyDeltaAy = (fTyDeltaTy < 0) ? (vv[2][1] * fTyPoint - vv[1][1]) : (vv[1][1] - vv[2][1] * fTyPoint);
			double fTyLengthA = calcTexc_Length(fTyDeltaAx, fTyDeltaAy);
			double fTyLengthT = calcTexc_Length(0, fTyDeltaTy * cyTex);
			tsy = fTyLengthA / fTyLengthT;
			BOOL bMirrorY = (iYc & 1) ? TRUE : FALSE;
			double fTyAngleA = calcTexc_AngleXYd(fTyDeltaAx, fTyDeltaAy);

			double fTxTyAngleA = calcDEGFrom0To360(fTxAngleA - fTyAngleA);

			angle = fTxAngleA;

			if (89 < fTxTyAngleA && fTxTyAngleA < 91) {
				// -sx +sy
				tsx = -tsx;
				angle += 180;
			} else if (269 < fTxTyAngleA && fTxTyAngleA < 271) {
				// +sx +sy
			} else {
				// ???
				return false;
			}
			angle = calcDEGFrom0To360(angle);

			double fEx, fEy;
			calcTexc_rotateTriPos(
				fAs[0][0][0], fAs[0][0][1],
				fAs[0][1][0], fAs[0][1][1],
				fAs[0][2][0], fAs[0][2][1],
				cxTex, cyTex, tsx, tsy, angle, fEx, fEy
				);
			fEx -= fTs[0][0][0];
			fEy -= fTs[0][0][1];
			if (fabs(fEx) < 1E-5) fEx = 0;
			if (fabs(fEy) < 1E-5) fEy = 0;
			tx = (int)floor((-1.0 * fEx * cxTex) + 0.5);
			ty = (int)floor((-1.0 * fEy * cyTex) + 0.5);
			tx -= (int)(tx / cxTex) * cxTex;
			ty -= (int)(ty / cyTex) * cyTex;
			while ((int)floor(tx + 0.5) < 0) tx += cxTex;
			while ((int)floor(ty + 0.5) < 0) ty += cyTex;
			while ((int)floor(tx + 0.5) >= (int)cxTex) tx -= cxTex;
			while ((int)floor(ty + 0.5) >= (int)cyTex) ty -= cyTex;

			return true;
		}
	}
	return false;
}
