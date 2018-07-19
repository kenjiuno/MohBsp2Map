
#pragma once

#pragma warning(disable: 4786)

#include <afxwin.h>
#include <afxtempl.h>
#include <set>
#include <map>
#include <list>
#include <math.h>
#include <algorithm>
#include "SizeVuff.h"
#include "Observer.h"

using namespace std;

//#define _DEBUG_ONLY_WORLD_SPAWN
//#define _DEBUG_TRACE_SURFACE_TRIS
//#define _DEBUG_TRACE_BRUSH_STAGE_POLYS
//#define _DEBUG_TRACE_JOIN_TRI
//#define _DEBUG_TRACE_JOIN_POLY_INPUT
#define _DEBUG_DENY_COMBINE_ENCLOSED_POLY
//#define _DEBUG_TRACE_CONVEX_TEST_AFTER

#define USE_NEW_JOIN_TRI_TEST

#define EPS 1E-5

/* /////////////////////////////////////////////////////////

///////////////////////////////////////////////////////// */

#pragma pack(push, 1)

namespace MB2MM
{
	// 
	enum n2fType {
		n2fPXZ,
		n2fPYZ,
		n2fPXY,
		n2fNXZ,
		n2fNYZ,
		n2fNXY,
	};

	// 
	typedef struct model3_t
	{	float mins[3];
		float maxs[3];
		DWORD face, n_faces;
		DWORD brush, n_brushes;
	}	model3_t;
	// 
	typedef struct surface3_t
	{	DWORD texture;
		DWORD effect;
		DWORD type;
		DWORD vertex;
		DWORD n_vertexes;
		DWORD meshvert;
		DWORD n_meshverts;
		DWORD lm_index;
		DWORD lm_start[2];
		DWORD lm_size[2];
		float lm_origin[3];
		float lm_vecs[2][3];
		float normal[3];
		DWORD size[2];
	}	surface3_t;
	// 
	typedef struct vertex3_t
	{	float position[3];
		float texcoord[2][2];
		float normal[3];
		BYTE color[4];
	}	vertex3_t;
	// 
	typedef struct tex3_t
	{	BYTE name[64];
		DWORD surf, contents;
	}	tex3_t;
	// 
	typedef struct terrain_t
	{	BYTE flags; // 0x00
		BYTE dummy1[0x23]; // 0x01
		__int8 x; // 0x24
		__int8 y; // 0x25
		__int16 z; // 0x26
		WORD texture; // 0x28
		BYTE dummy2[0x106]; // 0x2A
		BYTE m[9][9]; // 0x130
		// 0x181
	}	terrain_t;
};

#pragma pack(pop)

namespace MB2MM
{
	// 
	class myVtx_t
	{
	public:
		// 
		double x, y, z;

		// 
		myVtx_t(): x(0), y(0), z(0) { }
		// 
		myVtx_t(double x_, double y_, double z_): x(x_), y(y_), z(z_) { }
		// 
		myVtx_t operator *(double d) const
		{
			return myVtx_t(x * d, y * d, z * d);
		}
		// 
		myVtx_t operator /(double d) const
		{
			return myVtx_t(x / d, y / d, z / d);
		}
		// 
		myVtx_t operator +(const myVtx_t &src) const
		{
			return myVtx_t(x + src.x, y + src.y, z + src.z);
		}
		// 
		myVtx_t operator +=(const myVtx_t &src)
		{
			return *this = myVtx_t(x + src.x, y + src.y, z + src.z);
		}
		// 
		myVtx_t operator -(const myVtx_t &src) const
		{
			return myVtx_t(x - src.x, y - src.y, z - src.z);
		}
		// 
		bool operator ==(const myVtx_t &src) const
		{
			return (x == src.x && y == src.y && z == src.z) ? true : false;
		}
		// 
		bool operator !=(const myVtx_t &src) const
		{
			return !(this->operator ==(src));
		}
		// 
		bool IsZeroVector() const
		{
			return (x == 0 && y == 0 && z == 0) ? true : false;
		};
		// 
		double CalcLength() const
		{
			double f = sqrt(x * x + y * y);
			double g = sqrt(f * f + z * z);
			return g;
		}
		// 
		myVtx_t CalcNorm() const
		{
			double l = CalcLength();
			if (l == 0)
				return myVtx_t();
			return myVtx_t(x / l, y / l, z / l);
		}
		// 
		bool operator <(const myVtx_t &src) const
		{
			return CompareTo(src) < 0;
		}
		// 
		bool operator >(const myVtx_t &src) const
		{
			return CompareTo(src) > 0;
		}
		// 
		int CompareTo(const myVtx_t &src) const
		{
			if (x < src.x) return -1;
			if (x > src.x) return 1;
			if (y < src.y) return -1;
			if (y > src.y) return 1;
			if (z < src.z) return -1;
			if (z > src.z) return 1;
			return 0;
		}
		// 
		double Vec() const
		{
			return sqrt(x * x + y * y + z * z);
		}
		// 
		void Normalize()
		{
			double t = 1.0 / Vec();
			x *= t;
			y *= t;
			z *= t;
		}
		// 
		void Clear()
		{
			x = y = z = 0;
		}
		// 
		double DotProduct(const myVtx_t &v2) const
		{
			return DotProduct(*this, v2);
		}
		// 
		myVtx_t operator -() const
		{
			return myVtx_t(-x, -y, -z);
		}

		// 
		static double DotProduct(const myVtx_t &v1, const myVtx_t &v2)
		{
			return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
		}
		// 
		static myVtx_t CrossProduct(const myVtx_t &v1, const myVtx_t &v2)
		{
			myVtx_t t;

			t.x = v1.y * v2.z - v1.z * v2.y;
			t.y = v1.z * v2.x - v1.x * v2.z;
			t.z = v1.x * v2.y - v1.y * v2.x;

			return t;
		}
		// 
		static n2fType Normal2Facing(const myVtx_t &v);
	};
	// 
	class myVtx2_t : public myVtx_t
	{
	public:
		// 
		float tu, tv;

		// 
		myVtx2_t(): myVtx_t() { }
		// 
		myVtx2_t(double x, double y, double z): myVtx_t(x, y, z) { }
		// 
		myVtx2_t(double x, double y, double z, float tu, float tv): myVtx_t(x, y, z), tu(tu), tv(tv) { }
		// 
		explicit myVtx2_t(const myVtx_t &s): myVtx_t(s)
		{

		}
	};
//	// 
//	extern void calcN1Vec(myVtx_t &v);
	// 
	class myTri_t
	{
	public:
		// 
		myVtx2_t v[3];
		// 
		myTri_t()
		{

		}
		// 
		myTri_t(const myVtx2_t &v0, const myVtx2_t &v1, const myVtx2_t &v2)
		{
			v[0] = v0, v[1] = v1, v[2] = v2;
		}
		// 
		bool IsValidTri() const
		{
			if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0])
				return false;
			double fXmin = __min(v[0].x, __min(v[1].x, v[2].x));
			double fXmax = __max(v[0].x, __max(v[1].x, v[2].x));
			double fYmin = __min(v[0].y, __min(v[1].y, v[2].y));
			double fYmax = __max(v[0].y, __max(v[1].y, v[2].y));
			double fZmin = __min(v[0].z, __min(v[1].z, v[2].z));
			double fZmax = __max(v[0].z, __max(v[1].z, v[2].z));
			int iIdentical = 0;
			if (fabs(fXmin - fXmax) < EPS) iIdentical++;
			if (fabs(fYmin - fYmax) < EPS) iIdentical++;
			if (fabs(fZmin - fZmax) < EPS) iIdentical++;
			if (iIdentical >= 2)
				return false;
#if 0
			myVtx_t v10 = v[1] - v[0];
			myVtx_t v21 = v[2] - v[1];
			myVtx_t v02 = v[0] - v[2];
			if (fabs(fabs(myVtx_t::DotProduct(v10, v21)) - 1.0) < EPS)
				return false;
			if (fabs(fabs(myVtx_t::DotProduct(v21, v02)) - 1.0) < EPS)
				return false;
			if (fabs(fabs(myVtx_t::DotProduct(v02, v10)) - 1.0) < EPS)
				return false;
#endif
			return true;
		}
		// 
		void Norm(myVtx_t &t)
		{
			myVtx_t(v[0].x-v[1].x,v[0].y-v[1].y,v[0].z-v[1].z).CrossProduct(
				myVtx_t(v[2].x-v[0].x,v[2].y-v[0].y,v[2].z-v[0].z),
				t
				);
			t.Normalize();
		}
	};
	// 
	class myMtx_t
	{
	public:
		// 
		union {
			// 
			double m[4][4];
			// 
			double z[16];
		};
		// 
		void Identity()
		{
			for (char y = 0; y < 4; y++)
				for (char x = 0; x < 4; x++)
					m[y][x] = (x == y) ? 1.0 : 0;
		}
		// 
		friend myVtx_t operator *(const myVtx_t &v, const myMtx_t &r)
		{
			double xyzw[4] = {v.x, v.y, v.z, 1.0};
			myVtx_t t;
			for (char y = 0; y < 3; y++) {
				double a = 0;
				for (char x = 0; x < 4; x++) {
					a += xyzw[x] * r.m[x][y];
				}
				switch (y) {
				case 0: t.x = a; break;
				case 1: t.y = a; break;
				case 2: t.z = a; break;
				}
			}
			return t;
		}
		// 
		const myMtx_t &operator *=(const myMtx_t &v)
		{
			for (char y = 0; y < 4; y++) {
				for (char x = 0; x < 4; x++) {
					double a = 0;
					for (char c = 0; c < 4; c++) a += m[y][c] * v.m[c][x];
					m[y][x] = a;
				}
			}
			return *this;
		}
		// 
		myMtx_t operator *(const myMtx_t &v) const
		{
			myMtx_t t;
			for (char y = 0; y < 4; y++) {
				for (char x = 0; x < 4; x++) {
					double a = 0;
					for (char c = 0; c < 4; c++) a += m[y][c] * v.m[c][x];
					t.m[y][x] = a;
				}
			}
			return t;
		}
		// 
		void rotateX(double rad)
		{
			myMtx_t t;
			t.Identity();
			t.m[1][1] = cos(rad);
			t.m[1][2] = sin(rad);
			t.m[2][1] = -t.m[1][2];
			t.m[2][2] =  t.m[1][1];
			*this *= t;
		}
	};
	// 
	class myPoly_t
	{
	public:
		// 
		CArray<myVtx2_t, myVtx2_t> arrVerts;
		// 
		myVtx_t vtxNorm;
		// 
		CArray<myTri_t, myTri_t> arrTris;

		// 
		myPoly_t()
		{
			arrVerts.SetSize(0, 64);
			arrTris.SetSize(0, 15);
		}
		// 
		myPoly_t(const myPoly_t &s)
		{
			*this = s;
		}

		// 
		const myPoly_t &operator =(const myPoly_t &s)
		{
			vtxNorm = s.vtxNorm;
			arrVerts.Copy(s.arrVerts);
			arrTris.Copy(s.arrTris);
			return *this;
		}
		// 
		bool operator ==(const myPoly_t &s) const
		{
			UINT n = arrVerts.GetSize();
			if (s.arrVerts.GetSize() != n) return false;

			UINT i;
			for (i = 0; i < n; i++) {
				if (arrVerts[i] != s.arrVerts[i])
					return false;
			}
			return true;
		}

		// 
		BOOL IsEmpty() const
		{
			return (arrVerts.GetSize() == 0) ? TRUE : FALSE;
		}
		// 
		BOOL JoinTri(const myTri_t &s)
		{
			int iVert, nVerts = arrVerts.GetSize();
			ASSERT(nVerts == 0 || nVerts >= 3);
			if (nVerts == 0) {
				arrVerts.SetSize(3);
				for (iVert = 0; iVert < 3; iVert++) arrVerts[iVert] = s.v[iVert];
				OnDevAdd(3);
				return TRUE;
			}
			for (iVert = 0; iVert < nVerts; iVert++) {
				const UINT iVert0 = (iVert    );
				const UINT iVert1 = (iVert + 1) % nVerts;
				const myVtx_t &rV0 = arrVerts.ElementAt(iVert0);
				const myVtx_t &rV1 = arrVerts.ElementAt(iVert1);
				ASSERT(rV0 != rV1);
				if (rV0 == s.v[1] && rV1 == s.v[0]) {
					ASSERT(rV0 != s.v[2] && rV1 != s.v[2]);
					if (!JoinOk(iVert + 1, s.v[2])) break;
					arrVerts.InsertAt(iVert + 1, s.v[2]);
					OnDevAdd(iVert+1);
					RemoveRift(); // RemoveJoint(iVert+1);
					OnJoinDone(s);
					return TRUE;
				}
				if (rV0 == s.v[2] && rV1 == s.v[1]) {
					ASSERT(rV0 != s.v[0] && rV1 != s.v[0]);
					if (!JoinOk(iVert + 1, s.v[0])) break;
					arrVerts.InsertAt(iVert + 1, s.v[0]);
					OnDevAdd(iVert+1);
					RemoveRift(); // RemoveJoint(iVert+1);
					OnJoinDone(s);
					return TRUE;
				}
				if (rV0 == s.v[0] && rV1 == s.v[2]) {
					ASSERT(rV0 != s.v[1] && rV1 != s.v[1]);
					if (!JoinOk(iVert + 1, s.v[1])) break;
					arrVerts.InsertAt(iVert + 1, s.v[1]);
					OnDevAdd(iVert+1);
					RemoveRift(); // RemoveJoint(iVert+1);
					OnJoinDone(s);
					return TRUE;
				}
			}
			return FALSE;
		}
		// 
		bool ReduceNeighbor()
		{
			UINT n = arrVerts.GetSize();
			UINT m = n;
			ASSERT(3 <= n);
			if (n < 3) return false;
			UINT i;
			for (i = 0; i < n; ) {
				UINT v0 = (i + 0) % n;
				UINT v1 = (i + 1) % n;
				UINT v2 = (i + 2) % n;
				myVtx_t v1_0 = (arrVerts[v1] - arrVerts[v0]);
				myVtx_t v2_1 = (arrVerts[v2] - arrVerts[v1]);
				v1_0.Normalize();
				v2_1.Normalize();
				if (fabs(v1_0.DotProduct(v2_1) - 1.0) < EPS) {
					//!ASSERT(3 < n);
					if (n <= 3) return false;

					UINT fm;
					for (fm = 0; fm < i; fm++) arrVerts.Add(arrVerts[fm]);
					arrVerts.RemoveAt(0, i);

					//myVtx_t t0 = arrVerts[v0];
					//myVtx_t t1 = arrVerts[v1];
					//myVtx_t t2 = arrVerts[v2];

					arrVerts.RemoveAt(1);
					i = 0;
					n = arrVerts.GetSize();
				} else {
					i++;
				}
			}
			return true;
		}
		// 
		void OnDevFinalReview(UINT x)
		{
#ifdef _DEBUG_TRACE_BRUSH_STAGE_POLYS
			{
				TRACE0("=-=-=-=-=-=\n");
				TRACE1("#%d\n", x);
				UINT i, n;
				i, n = arrVerts.GetSize();
				for (i = 0; i < n; i++) {
					TRACE3("   %5d,%5d,%5d\n"
						, (int)arrVerts[i].x
						, (int)arrVerts[i].y
						, (int)arrVerts[i].z
						);
					Sleep(1);
				}
			}
#endif
		}
		// 
		void Reduce()
		{
#if 1
			ASSERT(FALSE);
#else
//			int nRepair;
//			do {
//				nRepair = 0;
//				UINT i, n = arrVerts.GetSize();
//				for (i = 0; i < n && n >= 2; i++) {
//					const myVtx_t &v0 = arrVerts.ElementAt((i + 0));
//					const myVtx_t &v1 = arrVerts.ElementAt((i + 1) % n);
//					if (v0 == v1) {
//						arrVerts.RemoveAt((i + 1) % n);
//						i--;
//						n--;
//						nRepair++;
//					}
//				}
//				for (i = 0; i < n && n >= 3; i++) {
//					const myVtx_t &v0 = arrVerts.ElementAt((i + 0));
//					const myVtx_t &v1 = arrVerts.ElementAt((i + 1) % n);
//					const myVtx_t &v2 = arrVerts.ElementAt((i + 2) % n);
//					if (v0 == v2) {
//						arrVerts.RemoveAt((i + 1) % n);
//						if (!(i < (i + 1) % n))
//							i--;
//						n--;
//						nRepair++;
//						arrVerts.RemoveAt((i + 0));
//						i--;
//						n--;
//						nRepair++;
//					}
//				}
//				for (i = 0; i < n && n >= 3; i++) {
//					const myVtx_t &v0 = arrVerts.ElementAt((i + 0));
//					const myVtx_t &v1 = arrVerts.ElementAt((i + 1) % n);
//					const myVtx_t &v2 = arrVerts.ElementAt((i + 2) % n);
//					myVtx_t c0 = v1 - v0;
//					myVtx_t c1 = v2 - v1;
//					calcN1Vec(c0);
//					calcN1Vec(c1);
//					double rv0 = acos(c0.x) / 3.1415926535897932384626433832795 * 180.0;
//					double rw0 = acos(c0.y) / 3.1415926535897932384626433832795 * 180.0;
//					if (c0.z < 0)
//						rv0 = 360 - rv0,
//						rw0 = 360 - rw0;
//					double rv1 = acos(c1.x) / 3.1415926535897932384626433832795 * 180.0;
//					double rw1 = acos(c1.y) / 3.1415926535897932384626433832795 * 180.0;
//					if (c1.z < 0)
//						rv1 = 360 - rv1,
//						rw1 = 360 - rw1;
//					if (fabs(rv0 - rv1) < 1.0 && fabs(rw0 - rw1) < 1.0) {
//						arrVerts.RemoveAt((i + 1) % n);
//						i--;
//						n--;
//						nRepair++;
//					}
//				}
//			} while (nRepair != 0);
#endif
		}
#if 0
//		// 
//		BOOL JoinPoly(myPoly_t &s)
//		{
//			const nPolys1 = arrVerts.GetSize();
//			const nPolys2 = s.arrVerts.GetSize();
//			for (int iPoly1 = 0; iPoly1 < nPolys1; iPoly1++) {
//				const myVtx_t &s0v0 = arrVerts.ElementAt( iPoly1               );
//				const myVtx_t &s0v1 = arrVerts.ElementAt((iPoly1 + 1) % nPolys1);
//				double fXd = s0v1.x - s0v0.x;
//				double fYd = s0v1.y - s0v0.y;
//				double fZd = s0v1.z - s0v0.z;
//				for (int iPoly2 = 0; iPoly2 < nPolys2; iPoly2++) {
//					const myVtx_t &s1v0 = s.arrVerts.ElementAt( iPoly2               );
//					const myVtx_t &s1v1 = s.arrVerts.ElementAt((iPoly2 + 1) % nPolys2);
//					if (s0v0 == s1v1 && s0v1 == s1v0) {
//						for (int iJoin1 = 0; iJoin1 < nPolys2 - 2; iJoin1++) {
//							const iInsertAt = (iPoly1 + 1 + iJoin1) % nPolys1;
//							arrVerts.InsertAt(iInsertAt, s.arrVerts.ElementAt((iPoly2 + 2 + iJoin1) % nPolys2));
//						}
//						return TRUE;
//					}
//					double fXf0 = clampVal((s1v0.x - s0v0.x) / fXd, 3);
//					double fYf0 = clampVal((s1v0.y - s0v0.y) / fYd, 3);
//					double fZf0 = clampVal((s1v0.z - s0v0.z) / fZd, 3);
//					double fXf1 = clampVal((s1v1.x - s0v0.x) / fXd, 3);
//					double fYf1 = clampVal((s1v1.y - s0v0.y) / fYd, 3);
//					double fZf1 = clampVal((s1v1.z - s0v0.z) / fZd, 3);
//					BOOL bOnl0 = (fXf0 == fYf0 && fYf0 == fZf0) ? TRUE : FALSE;
//					BOOL bOnl1 = (fXf1 == fYf1 && fYf1 == fZf1) ? TRUE : FALSE;
//					if (!bOnl0 && !bOnl1)
//						continue;
//					printf("");
//				}
//			}
//			return FALSE;
//		}
#endif
		// 
		void OnDevAdd(UINT x)
		{
#ifdef _DEBUG_TRACE_JOIN_TRI
			TRACE0("---\n");
			TRACE1("?:%p\n", this);
			UINT i, n = arrVerts.GetSize();
			for (i = 0; i < n; i++) {
				::AfxTrace("%c  %5d,%5d,%5d\n"
					, (i == x) ? '+' : ' '
					, (int)arrVerts[i].x
					, (int)arrVerts[i].y
					, (int)arrVerts[i].z
					);
				Sleep(1);
			}
			TRACE0("===\n");
#endif
		}
#if 0
//		// 
//		void RemoveJoint(UINT i)
//		{
//			if (arrVerts.GetSize() <= i + 2) return;
//			if (arrVerts[i] != arrVerts[i + 2]) return;
//
//			arrVerts.RemoveAt(i, 2);
//
//			OnDevAdd(-1);
//		}
#endif
		// 
		bool JoinPoly(const myPoly_t &s)
		{
			const UINT n1 = arrVerts.GetSize();
			const UINT n2 = s.arrVerts.GetSize();
			const UINT nn = __min(n1, n2);

#if _DEBUG_TRACE_JOIN_POLY_INPUT
			{
				TRACE0("==========\n");
				UINT i, n;
				TRACE1("   %p\n", &arrVerts);
				i, n = arrVerts.GetSize();
				for (i = 0; i < n; i++) {
					TRACE3("1: %5d,%5d,%5d\n"
						, (int)arrVerts[i].x
						, (int)arrVerts[i].y
						, (int)arrVerts[i].z
						);
					Sleep(1);
				}
				TRACE0("---\n");
				TRACE1("   %p\n", &s.arrVerts);
				i, n = s.arrVerts.GetSize();
				for (i = 0; i < n; i++) {
					TRACE3("2: %5d,%5d,%5d\n"
						, (int)s.arrVerts[i].x
						, (int)s.arrVerts[i].y
						, (int)s.arrVerts[i].z
						);
					Sleep(1);
				}
				TRACE0("---\n");
			}
#endif

			bool fReady = false;
			UINT iStart1, iStart2, nLen;

			// try to join two polys with joint sides.
			UINT f1;
			for (f1 = 0; f1 < n1; f1++) {
				UINT f2a, f2s, f2t;
				for (f2a = 0; f2a < n2; f2a++) {
					if (arrVerts[f1] != s.arrVerts[f2a]) continue;

					for (f2s = 0; f2s < nn; f2s++) {
						if (arrVerts[CalcMod(n1 + f1 - f2s      , n1)] == s.arrVerts[CalcMod(     f2a + f2s      , n2)]) continue;

						break;
					}

					ASSERT(f2s != 0);
					f2s--;

					for (f2t = 0; f2t < nn; f2t++) {
						if (arrVerts[CalcMod(n1 + f1 - f2s + f2t, n1)] == s.arrVerts[CalcMod(n2 + f2a + f2s - f2t, n2)]) continue;

						break;
					}

					if (f2t < 2)
						continue;

					if (fReady) {
						UINT t1 = f1 + n1;
						if ((iStart1 <= f1 && f1 < iStart1 + nLen) || (iStart1 <= t1 && t1 < iStart1 + nLen)) {
							UINT t2 = f2a + n2;
							if ((n2 + iStart2 - nLen + 1 <= f2a && f2a <= n2 + iStart2) || (n2 + iStart2 - nLen + 1 <= t2 && t2 <= n2 + iStart2)) {
								continue;
							} else {
								printf("");
							}
						} else {
							printf("");
						}

						// poly cannot be connected with two or more joints.
						//!ASSERT(false); // very rarely case
						return false;
					}

#ifdef _DEBUG_DENY_COMBINE_ENCLOSED_POLY
					if (f2t == nn) {
						// poly cannot be enclosed completely.
						ASSERT(false); // rarely case
						return false;
					}
#endif

					fReady = true;
					iStart1 = CalcMod(n1 + f1 - f2s, n1);
					iStart2 = CalcMod(    f2a + f2s, n2);
					nLen = f2t;

					f2a += f2t;

					//!ASSERT(nLen != n2);
				}
			}

			if (!fReady) return false; // usual case

			if (nLen == n2) return true; // enclosing complete identical region

			UINT nAdd = n2 - nLen + 2;

			UINT fm;
			for (fm = 0; fm < iStart1; fm++) arrVerts.Add(arrVerts[fm]);
			arrVerts.RemoveAt(0, iStart1);

			iStart1 = 0;

			arrVerts.RemoveAt(iStart1, nLen);
			arrVerts.InsertAt(iStart1, myVtx2_t(), nAdd);

			UINT fa;
			for (fa = 0; fa < nAdd; fa++) {
				arrVerts[iStart1 + fa] = s.arrVerts[(iStart2 + fa) % n2];
			}

			RemoveRift();
			OnJoinDone(s);
			return true;
		}
		// 
		void RemoveRift()
		{
			UINT i, n = arrVerts.GetSize();
			ASSERT(3 <= n);
			if (n < 3) return;

			for (i = 0; i < n; ) {
				if (arrVerts[i] == arrVerts[(i + 1) % n]) {
					UINT fm;
					for (fm = 0; fm < i; fm++) arrVerts.Add(arrVerts[fm]);
					arrVerts.RemoveAt(0, i);

					arrVerts.RemoveAt(0);

					i = 0;
					n--;
					if (n < 3) break;
					continue;
				}
				if (arrVerts[i] == arrVerts[(i + 2) % n]) {
					UINT fm;
					for (fm = 0; fm < i; fm++) arrVerts.Add(arrVerts[fm]);
					arrVerts.RemoveAt(0, i);

					arrVerts.RemoveAt(0, 2);

					i = 0;
					n -= 2;
					if (n < 3) break;
					continue;
				}
				i++;
			}
			ASSERT(3 <= n);
		}
#ifdef USE_NEW_JOIN_TRI_TEST
		// 
		bool JoinOk(UINT i, const myVtx_t &vi)
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			const myVtx_t &v0 = arrVerts[(n + i - 3) % n];
			const myVtx_t &v1 = arrVerts[(n + i - 2) % n];
			const myVtx_t &v2 = arrVerts[(n + i - 1) % n];
			//---//
			const myVtx_t &v3 = arrVerts[(    i    ) % n];
			const myVtx_t &v4 = arrVerts[(    i + 1) % n];
			const myVtx_t &v5 = arrVerts[(    i + 2) % n];

			myVtx_t vx1 = myVtx_t::CrossProduct(    vtxNorm, v2 - v1);
			if (myVtx_t::DotProduct(vx1, vi) < -EPS)
				return false;

			myVtx_t vx2 = myVtx_t::CrossProduct(    vtxNorm, vi - v2);
			if (myVtx_t::DotProduct(vx2, v3) < -EPS)
				return false;

			myVtx_t vx3 = myVtx_t::CrossProduct(    vtxNorm, v3 - vi);
			if (myVtx_t::DotProduct(vx3, v4) < -EPS)
				return false;

			myVtx_t vx4 = myVtx_t::CrossProduct(    vtxNorm, v4 - v3);
			if (myVtx_t::DotProduct(vx4, v5) < -EPS)
				return false;

			return true;
		}
#else 
		// 
		bool JoinOk(UINT i, const myVtx_t &v)
		{
			return true
				&& JoinOk1(i, v)
				&& JoinOk2(i, v)
				&& JoinOk3(i, v)
				&& JoinOk4(i, v)
				;
		}
		// 
		bool JoinOk1(UINT i, const myVtx_t &v3)
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			const myVtx_t &v0 = arrVerts[(n + i - 3) % n];
			const myVtx_t &v1 = arrVerts[(n + i - 2) % n];
			const myVtx_t &v2 = arrVerts[(n + i - 1) % n];

			myVtx_t vx1 = myVtx_t::CrossProduct(v2 - v1, v0 - v1);
			myVtx_t vx2 = myVtx_t::CrossProduct(    vx1, v2 - v1);

			double r = myVtx_t::DotProduct(vx2, v3);

			if (0 < r)
				return true;
			return false;
		}
		// 
		bool JoinOk2(UINT i, const myVtx_t &v2)
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			const myVtx_t &v0 = arrVerts[(n + i - 2) % n];
			const myVtx_t &v1 = arrVerts[(n + i - 1) % n];

			const myVtx_t &v3 = arrVerts[(n + i - 0) % n];

			myVtx_t vx1 = myVtx_t::CrossProduct(v2 - v1, v0 - v1);
			myVtx_t vx2 = myVtx_t::CrossProduct(    vx1, v2 - v1);

			double r = myVtx_t::DotProduct(vx2, v3);

			if (0 < r)
				return true;
			return false;
		}
		// 
		bool JoinOk3(UINT i, const myVtx_t &v1)
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			const myVtx_t &v0 = arrVerts[(n + i - 1) % n];

			const myVtx_t &v2 = arrVerts[(n + i - 0) % n];
			const myVtx_t &v3 = arrVerts[(n + i + 1) % n];

			myVtx_t vx1 = myVtx_t::CrossProduct(v2 - v1, v0 - v1);
			myVtx_t vx2 = myVtx_t::CrossProduct(    vx1, v2 - v1);

			double r = myVtx_t::DotProduct(vx2, v3);

			if (0 < r)
				return true;
			return false;
		}
		// 
		bool JoinOk4(UINT i, const myVtx_t &v0)
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			const myVtx_t &v1 = arrVerts[(n + i - 0) % n];
			const myVtx_t &v2 = arrVerts[(n + i + 1) % n];
			const myVtx_t &v3 = arrVerts[(n + i + 2) % n];

			myVtx_t vx1 = myVtx_t::CrossProduct(v2 - v1, v0 - v1);
			myVtx_t vx2 = myVtx_t::CrossProduct(    vx1, v2 - v1);

			double r = myVtx_t::DotProduct(vx2, v3);

			if (0 < r)
				return true;
			return false;
		}
#endif
		// 
		bool JoinPolySafely(const myPoly_t &s)
		{
			CArray<myVtx2_t, myVtx2_t> arr;
			arr.Copy(this->arrVerts);

			ASSERT(IsConvexShape());

			if (JoinPoly(s)) {
				if (IsConvexShape()) {
					return true;
				}
				arrVerts.Copy(arr);
			}
			return false;
		}
		// 
		bool IsConvexShape()
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);
			if (n < 3) return true;

			UINT i;
			for (i = 0; i < n; i++) {
				const myVtx_t &v0 = arrVerts[(i    )    ];
				const myVtx_t &v1 = arrVerts[(i + 1) % n];
				const myVtx_t &v2 = arrVerts[(i + 2) % n];

				myVtx_t vx2 = myVtx_t::CrossProduct(vtxNorm, v1 - v0);

				if (myVtx_t::DotProduct(vx2, v2 - v1) < -EPS) {
#ifdef _DEBUG_TRACE_CONVEX_TEST_AFTER
					{
						TRACE0("=-=-=-=-=-=\n");
						TRACE1("  %p\n", this);
						UINT t, n;
						t, n = arrVerts.GetSize();
						for (t = 0; t < n; t++) {
							TRACE3("   %5d,%5d,%5d\n"
								, (int)arrVerts[t].x
								, (int)arrVerts[t].y
								, (int)arrVerts[t].z
								);
							Sleep(1);
						}
					}
					double w = myVtx_t::DotProduct(vx2, v2 - v1);
					myVtx_t t10 = v1 - v0; t10.Normalize();
					myVtx_t t21 = v2 - v1; t21.Normalize();
#endif
					return false;
				}
			}

			return true;
		}
		// 
		bool IsNeedleBrush()
		{
			const UINT n = arrVerts.GetSize();
			ASSERT(3 <= n);

			UINT i = 0;
			{
				myVtx_t v0 = arrVerts[(i    )    ];
				myVtx_t v1 = arrVerts[(i + 1) % n];
				myVtx_t v2 = arrVerts[(i + 2) % n];

				myVtx_t v21 = v2 - v1;
				myVtx_t v10 = v1 - v0;

				myVtx_t vUpa = myVtx_t::CrossProduct(v10, v21);

				v21.x = (int)(0.5 +v21.x);
				v21.y = (int)(0.5 +v21.y);
				v21.z = (int)(0.5 +v21.z);
				v10.x = (int)(0.5 +v10.x);
				v10.y = (int)(0.5 +v10.y);
				v10.z = (int)(0.5 +v10.z);

				myVtx_t vUpb = myVtx_t::CrossProduct(v10, v21);

				vUpa.Normalize();
				vUpb.Normalize();

				double x = myVtx_t::DotProduct(vUpa, vUpb);

				if (x < 0.98) return true;
			}

			return false;
		}

		// 
		static UINT CalcMod(int x, int y)
		{
			if (x < 0) {
				return y - (x % y);
			}

			return x % y;
		}

	private:
		// 
		double clampVal(double d, int n)
		{
			int x = n;
			for (; x > 0; x--)
				d *= 10.0;
			if (d < 0)
				d = ceil(d);
			else
				d = floor(d);
			for (; n > 0; n--)
				d /= 10.0;
			return d;
		}
		// 
		double fix1P0(double d)
		{
			if (memcmp(&d, "\xFF\xFF\xFF\xFF\xFF\xFF\xEF\x3F", 8) == 0)
				memcpy(&d, "\x00\x00\x00\x00\x00\x00\xF0\x3F", 8);
			if (memcmp(&d, "\xFF\xFF\xFF\xFF\xFF\xFF\xEF\xBF", 8) == 0)
				memcpy(&d, "\x00\x00\x00\x00\x00\x00\xF0\xBF", 8);
			return d;
		}

	protected:
		// 
		inline void OnJoinDone(const myTri_t &t)
		{
			arrTris.Add(t);
		}
		// 
		inline void OnJoinDone(const myPoly_t &t)
		{
			arrTris.Append(t.arrTris);
		}
	};
	// 
	class myPolyKey_t
	{
	public:
		// 
		DWORD texture;
		// 
		myVtx_t norm;

		// 
		int CompareTo(const myPolyKey_t &s) const
		{
			if (texture < s.texture) return -1;
			if (texture > s.texture) return 1;
			if (norm < s.norm) return -1;
			if (norm > s.norm) return 1;
			return 0;
		}
		// 
		bool operator <(const myPolyKey_t &s) const
		{
			return CompareTo(s) < 0;
		}
	};
	// 
	class myPolySet_t
	{
	public:
		// 
		myPoly_t z;

		// 
		bool operator ==(const myPolySet_t &s) const
		{
			return z == s.z;
		}
	};
	// 
	class myPolySetList_t : public multimap<myPolyKey_t, myPolySet_t>
	{
	public:
		// 
		bool JoinTri(const myTri_t &t, const myVtx_t &norm, UINT texture, bool fReduce = true);
		// 
		void JoinPolys();
		// 
		void ReduceNeighbor();
		// 
		void OnDevFinalReview();

		// 
		void Close()
		{
			clear();
		}
	};
	// 
	class myTexture_t
	{
	public:
		// 
		CString strName;
		// 
		DWORD surf, contents;

		// 
		CString FormatName() const
		{
			if (_tcsncmp(strName, "textures/", 9) == 0)
				return (LPCTSTR)strName + 9;
			return strName;
		}
		// 
		CString FormatContents() const
		{
			CString strRet;
			if (contents & 0x8000000) {
				strRet += "+surfaceparm detail ";
			}
			return strRet;
		}
	};
};

namespace MB2MM
{
	// 
	bool ReadInt(LPCTSTR psz, int &x);
	// 
	struct TextureMetric {
		// 
		UINT cx;
		// 
		UINT cy;
		// 
		CString strName;
	};
	// 
	typedef set<CString> TriggerSet;
	// 
	typedef map<CString, TextureMetric> TextureMetricMap;
	// 
	struct Knowledge {
		// 
		TriggerSet triggers;
		// 
		TextureMetricMap m;

		// 
		bool LoadTexMetrFiles(CString strDir);
		// 
		bool LoadBBoxFiles(CString strDir);

		// 
		bool IsTriggerEntity(CString strKey) const
		{
			if (triggers.find(strKey) != triggers.end())
				return true;
			return false;
		}

		// 
		void Close()
		{
			triggers.clear();
			m.clear();
		}
	};
	// 
	class ReadBSP
	{
		// 
		CStdioFile f;

	public:
		// 
		~ReadBSP()
		{
			Close();
		}
		// 
		void Close()
		{
			f.Abort();
		}
		// 
		bool Open(LPCTSTR pszIn)
		{
			Close();

			if (f.Open(pszIn, 0 |CFile::modeRead |CFile::shareDenyWrite |CFile::typeBinary)) {
				return true;
			}
			return false;
		}
		// 
		bool ReadEntry(int i, SizeBuff &sb)
		{
			FILE *f = this->f.m_pStream;

			if (fseek(f, 12 +8*i, SEEK_SET) != 0)
				return false;
			DWORD x[2];
			if (fread(x, 8, 1, f) != 1)
				return false;
			if (x[1] == 0) {
				sb.Free();
				return true;
			}
			if (fseek(f, x[0], SEEK_SET) != 0)
				return false;
			if (!sb.Alloc(x[1]))
				return false;
			if (fread(sb.GetData(), x[1], 1, f) != 1)
				return false;
			return true;
		}
	};
	// 
	struct TextParsing {
		// 
		LPCSTR psz;
		// 
		int iPos, nLen;
		// 
		char c;

		// 
		TextParsing()
		{

		}
		// 
		TextParsing(LPCSTR psz)
			: psz(psz), iPos(0), nLen(strlen(psz))
		{

		}
		// 
		bool Next()
		{
			if (iPos < nLen) {
				c = psz[iPos];
				iPos++;
				return true;
			}
			c = 0;
			return false;
		}
		// 
		bool IsEOF() const
		{
			return !(iPos < nLen);
		}
		// 
		void SkipWs()
		{
			if (!IsEOF()) {
				while (isspace((BYTE)c) && Next());
			}
		}
		// 
		bool ReadOpenBracket()
		{
			SkipWs();
			if (c == '{') {
				Next();
				return true;
			}
			return false;
		}
		// 
		bool ReadCloseBracket()
		{
			SkipWs();
			if (c == '}') {
				Next();
				return true;
			}
			return false;
		}
		// 
		int ReadQuotedText(CString &str)
		{
			str.Empty();
			SkipWs();
			if (c == '\"') {
				if (!Next())
					return -1;
				while (c != '\"') {
					str += c;
					if (!Next())
						return -1;
				}
				Next();
				return +1;
			}
			return 0;
		}

		// 
		bool ParseModelNo(int &x)
		{
			if (false
				|| !Next()
				|| c != '*'
				|| !ReadInt(psz +iPos, x)
			) {
				return false;
			}
			return true;
		}
		// 
		bool ParseVertex(myVtx_t &t)
		{
			int v[3];
			if (false
				|| sscanf(psz, "%d %d %d", v+0, v+1, v+2) != 3
			) {
				return false;
			}
			t.x = v[0];
			t.y = v[1];
			t.z = v[2];
			return true;
		}
	};
	// 
	struct Entity1 : map<CString, CString> {
		// 
		CString className;
		// 
		CString origin;
		// 
		int modelNo;
		// 
		myVtx_t v_origin;
		// 
		CString modelName;

		// 
		void Finalize()
		{
			className = GetValue("classname");
			origin = GetValue("origin");
			if (TextParsing(GetValue("model")).ParseModelNo(modelNo)) { erase("model"); } else { modelNo = -1; }
			if (TextParsing(origin).ParseVertex(v_origin)) { } else { v_origin.Clear(); }

			modelName = GetValue("model");
		}
		// 
		CString GetValue(const CString &strKey)
		{
			iterator
				iterPos = find(strKey),
				iterEnd = end();
			if (iterPos != iterEnd)
				return iterPos->second;
			return _T("");
		}
		// 
		bool AddValue(const CString &strKey, const CString &strVal)
		{
			iterator
				iterPos = find(strKey),
				iterEnd = end();
			if (iterPos == iterEnd) {
				(*this)[strKey] = strVal;
				return true;
			}
			return false;
		}
	};
	// 
	struct ReadEntity {
		// 
		CArray<Entity1, Entity1> entities;

		// 
		void Close()
		{
			entities.RemoveAll();
		}

		// 
		bool Read(const SizeBuff &sb)
		{
			TextParsing z;
			z.psz = (LPCSTR)sb.GetData();
			z.iPos = 0;
			z.nLen = sb.GetSize();
			return Read(z);
		}
		// 
		bool Read(TextParsing &z)
		{
			list<Entity1> t;
			if (!z.Next())
				return false;
			while (true) {
				if (!z.ReadOpenBracket()) {
					if (z.IsEOF())
						break;
					return false;
				}
				t.push_back(Entity1());
				while (true) {
					CString strKey;
					int r = z.ReadQuotedText(strKey);
					if (r < 0)
						return false;
					if (r == 0) {
						t.back().Finalize();
						break;
					}
					CString strVal;
					if (z.ReadQuotedText(strVal) <= 0)
						return false;
					t.back().AddValue(strKey, strVal);
				}
				if (!z.ReadCloseBracket()) {
					return false;
				}
			}
			entities.SetSize(t.size());
			for (UINT i = 0; t.size() != 0; i++, t.pop_front()) entities[i] = t.front();
			return true;
		}
	};
	// 
	struct Texfc {
		// 
		int tx, ty;
		// 
		double angle;
		// 
		double tsx, tsy;

		// 
		Texfc() : tx(0), ty(0), angle(0), tsx(1.0), tsy(1.0)
		{

		}
	};
	// 
	struct BrushAttr : Texfc {
		// 
		CString strName;
		// 
		DWORD f[3];
		// 
		CString strFlags;

		// 
		BrushAttr()
		{
			ZeroMemory(f, sizeof(f));
		}
		// 
		void SetTexture(const myTexture_t &s)
		{
			strName = s.FormatName();
			strFlags = s.FormatContents();
		}
		// 
		void UndefTexture(bool fClearFlags = true)
		{
			strName = "common/black";

			if (fClearFlags) strFlags.Empty();

			tx = ty = 0;
			angle = 0;
			tsx = tsy = 1;
			f[0] = f[1] = f[2] = 0;
		}
	};
	// 
	struct BrushTri : BrushAttr {
		// 
		myTri_t tri;

		// 
		BrushTri() { }
		// 
		BrushTri(const BrushAttr &s)
		{
			BrushAttr::operator =(s);
		}
	};
	// 
	struct PatchDef2Attr {
		// 
		CString strName;
		// 
		UINT cx, cy;
		// 
		DWORD f[3];
		// 
		CString strFlags;

		// 
		void SetTexture(const myTexture_t &s)
		{
			strName = s.FormatName();
			strFlags = s.FormatContents();
		}
		// 
		void UndefTexture()
		{
			strName = "common/black";
			strFlags.Empty();
		}
	};
	// 
	struct PatchDef2Vtx {
		// 
		myVtx2_t v;
	};
	// 
	struct PatchDef2VtxArray {
		// 
		CArray<PatchDef2Vtx, PatchDef2Vtx &> z;

		// 
		PatchDef2VtxArray()
		{

		}
		// 
		PatchDef2VtxArray(const PatchDef2VtxArray &s)
		{
			*this = s;
		}
		// 
		const PatchDef2VtxArray &operator =(const PatchDef2VtxArray &s)
		{
			z.Copy(s.z);

			return *this;
		}
	};
	// 
	struct PatchDef2Set {
		// 
		PatchDef2Attr a;
		// 
		PatchDef2VtxArray z;
	};
	// 
	struct PatchDef2SetList : list<PatchDef2Set> {
		// 
		void Close()
		{
			clear();
		}
	};
	// 
	struct MatchTri {
		// 
		myVtx_t v[3];

		// 
		bool operator <(const MatchTri &s) const
		{
			return memcmp(v, s.v, sizeof(v)) < 0;
		}
	};
	// 
	struct MatchTriList : set<MatchTri> {
		// 
		bool Exist(const MatchTri &s) const
		{
			return find(s) != end();
		}
		// 
		bool Add(const MatchTri &t)
		{
			if (Exist(t)) return false;

			MatchTri x;
			char c;
			for (c = 0; c < 3; c++) {
				x.v[(c + 0) % 3] = t.v[0];
				x.v[(c + 1) % 3] = t.v[1];
				x.v[(c + 2) % 3] = t.v[2];
				VERIFY(insert(x).second);
			}
			return true;
		}
	};
	// 
	struct LODTerrainHead {
		// 
		UINT cx, cy, f;
		// 
		myVtx_t origin;

		// 
		LODTerrainHead()
		{
			cx = cy = f = 0;
		}
	};
	// 
	struct LODTerrainPart {
		// 
		DWORD f1[2];
		// 
		CString strName;
		// 
		int tx, ty;
		// 
		double angle;
		// 
		double tsx, tsy;
		// 
		DWORD f2;
		// 
		DWORD f3[3];
		// 
		CString strFlags;

		// 
		LODTerrainPart()
		{
			ZeroMemory(f1, sizeof(f1));
			tx = ty = 0;
			angle = 0;
			tsx = tsy = 1.0;
			f2 = 0;
			ZeroMemory(f3, sizeof(f3));
		}
		// 
		void SetTexture(const myTexture_t &s)
		{
			strName = s.FormatName();
			strFlags = s.FormatContents();
		}
		// 
		void UndefTexture()
		{
			strName = "common/black";
			strFlags.Empty();
		}
	};
	// 
	struct LODTerrainUnit {
		// 
		int z;
		// 
		CString str1;
		// 
		CString str2;

		// 
		LODTerrainUnit()
		{
			z = 0;
		}
	};
	// 
	class WriteMap
	{
		// 
		FILE *f;
		// 
		int iIndent;

		// 
		WriteMap(const WriteMap &);
		// 
		void operator =(const WriteMap &);

		// 
		inline static int vc2i(double f)
		{
			return (int)(f + 0.5);
		}

	public:
		// 
		WriteMap()
			: f(NULL)
		{

		}
		// 
		~WriteMap()
		{
			Close();
		}
		// 
		void Close()
		{
			if (f)
				fclose(f), f = NULL;

			iIndent = 0;
		}
		// 
		bool Open(LPCTSTR psz)
		{
			Close();

			if (f = fopen(psz, "wt")) {
				return true;
			}
			return false;
		}
		// 
		bool WriteOpenBracket()
		{
			if (false
				|| !FillIndent()
				|| EOF == fputs("{\n", f)
			) {
				return false;
			}
			iIndent++;
			return true;
		}
		// 
		bool WriteCloseBracket()
		{
			iIndent--;
			if (false
				|| !FillIndent()
				|| EOF == fputs("}\n", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteAttr(LPCSTR pszKey, LPCSTR pszVal)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f, "\"%s\" \"%s\"\n", pszKey, pszVal) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteCommentNo(LPCSTR pszWhat, int i)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f, "// %s %d\n", pszWhat, i) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteBrush(const BrushTri &t)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f
					, "( %d %d %d ) ( %d %d %d ) ( %d %d %d ) %s %d %d %lf %lf %lf %u %u %u %s\n"
					, (int)vc2i(t.tri.v[0].x)
					, (int)vc2i(t.tri.v[0].y)
					, (int)vc2i(t.tri.v[0].z)
					, (int)vc2i(t.tri.v[1].x)
					, (int)vc2i(t.tri.v[1].y)
					, (int)vc2i(t.tri.v[1].z)
					, (int)vc2i(t.tri.v[2].x)
					, (int)vc2i(t.tri.v[2].y)
					, (int)vc2i(t.tri.v[2].z)
					, (LPCSTR)t.strName
					, (int)t.tx
					, (int)t.ty
					, t.angle
					, t.tsx
					, t.tsy
					, t.f[0]
					, t.f[1]
					, t.f[2]
					, t.strFlags
					) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool FillIndent()
		{
			for (int n = iIndent; n > 0; n--)
				if (EOF == fputc(' ', f))
					return false;
			return true;
		}
		// 
		bool WritePatchDef2()
		{
			if (false
				|| !FillIndent()
				|| EOF == fputs("patchDef2\n", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteOpenMBracket()
		{
			if (false
				|| !FillIndent()
				|| EOF == fputs("(\n", f)
			) {
				return false;
			}
			iIndent++;
			return true;
		}
		// 
		bool WriteCloseMBracket()
		{
			iIndent--;
			if (false
				|| !FillIndent()
				|| EOF == fputs(")\n", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteMBracketBegin()
		{
			if (false
				|| !FillIndent()
				|| EOF == fputs("(", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteMBracketEnd()
		{
			if (false
				|| EOF == fputs(" )\n", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WritePatchDef2Attr(const PatchDef2Attr &a)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f
					, "%s\n"
					, (LPCSTR)a.strName
					) < 0
				|| !FillIndent()
				|| fprintf(f
					, "( %d %d %u %u %u %s )\n"
					, (int)a.cy
					, (int)a.cx
					, (UINT)a.f[0]
					, (UINT)a.f[1]
					, (UINT)a.f[2]
					, (LPCSTR)a.strFlags
					) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WritePatchDef2Vertex(const PatchDef2Vtx &t)
		{
			if (false
				|| fprintf(f
					, " ( %d %d %d %lf %lf )"
					, vc2i(t.v.x)
					, vc2i(t.v.y)
					, vc2i(t.v.z)
					, (double)t.v.tu
					, (double)t.v.tv
					) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteTerrainDef()
		{
			if (false
				|| !FillIndent()
				|| EOF == fputs("terrainDef\n", f)
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteTerrainDef2Head(const LODTerrainHead &a)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f
					, "%u %u %u\n"
					, a.cx
					, a.cy
					, a.f
					) < 0
				|| !FillIndent()
				|| fprintf(f
					, "%lf %lf %lf\n"
					, (double)a.origin.x
					, (double)a.origin.y
					, (double)a.origin.z
					) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteTerrainDef2Part(const LODTerrainPart &a)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f, "%u %u ( %s %d %d %5.2lf %u %f %f %u %u %u%s )\n"
				, (UINT)a.f1[0]
				, (UINT)a.f1[1]
				, (LPCSTR)a.strName
				, a.tx
				, a.ty
				, a.angle
				, (UINT)a.f2
				, a.tsx
				, a.tsy
				, (UINT)a.f3[0]
				, (UINT)a.f3[1]
				, (UINT)a.f3[2]
				, (LPCSTR)a.strFlags
				) < 0
			) {
				return false;
			}
			return true;
		}
		// 
		bool WriteTerrainDef2Unit(const LODTerrainUnit &a)
		{
			if (false
				|| !FillIndent()
				|| fprintf(f, "%f (%s ) (%s )\n"
				, (double)a.z
				, (LPCSTR)a.str1
				, (LPCSTR)a.str2
				) < 0
			) {
				return false;
			}
			return true;
		}
	};
	// 
	struct Optz {
		// 
		int nThick;
		// 
		int nReduction;
		// 
		bool fUniqueTris;
		// 
		bool fRemoveSharp;
		// 
		bool fRecoverTexfc;
		// 
		bool fSkew;
	};
	// 
	class Decompiler
	{
		// 
		SizeBuff
			sbTexture,
			sbSurface,
			sbVertex,
			sbModel,
			sbEntity,
			sbMesh,
			sbLOD;
		// 
		Knowledge &kb;
		// 
		ReadEntity entityset;
		// 
		WriteMap wm;
		// 
		CArray<myTexture_t, myTexture_t> fineTextures;
		// 
		myPolySetList_t polysets;
		// 
		Optz optz;
		// 
		DataWatcher &mo;

		// 
		PBYTE pModels;
		// 
		PBYTE pSurfaces;
		// 
		PBYTE pVertices;
		// 
		PBYTE pTextures;
		// 
		PBYTE pMeshes;
		// 
		PBYTE pLODs;

		// 
		UINT nModels;
		// 
		UINT nSurfaces;
		// 
		UINT nVertices;
		// 
		UINT nTextures;
		// 
		UINT nMeshes;
		// 
		UINT nLODs;

		// 
		model3_t *GetModel(UINT i)
		{
			if (i < nModels) return (model3_t *)(pModels +40*i);
			return NULL;
		}
		// 
		surface3_t *GetSurface(UINT i)
		{
			if (i < nSurfaces) return (surface3_t *)(pSurfaces +108*i);
			return NULL;
		}
		// 
		vertex3_t *GetVertex(UINT i)
		{
			if (i < nVertices) return (vertex3_t *)(pVertices +44*i);
			return NULL;
		}
		// 
		tex3_t *GetTex(UINT i)
		{
			if (i < nTextures) return (tex3_t *)(pTextures +140*i);
			return NULL;
		}
		// 
		DWORD *GetMesh(UINT i)
		{
			if (i < nMeshes) return (DWORD *)(pMeshes +4*i);
			return NULL;
		}
		// 
		myTexture_t *GetFineTex(UINT i)
		{
			if (i < nTextures) return &fineTextures[i];
			return NULL;
		}
		// 
		bool FindTexture(UINT i, TextureMetric &tm);
		// 
		void EnsureTexture(UINT i, TextureMetric &tm);
		// 
		terrain_t *GetLODt(UINT i)
		{
			if (i < nLODs) return (terrain_t *)(pLODs +388*i);
			return NULL;
		}

	public:
		// 
		Decompiler(Knowledge &kb, DataWatcher &mo, const Optz &optz)
			: kb(kb)
			, mo(mo)
			, brushNo(0)
			, optz(optz)
		{

		}
		// 
		bool Decompile(CString strIn, CString strOut);
		// 
		void Close()
		{
			entityset.Close();
			wm.Close();
			fineTextures.RemoveAll();
			polysets.Close();

			brushNo = 0;
		}

	protected:
		// 
		int brushNo;

		// 
		bool writeModel(UINT entityNo, int iModel);
		// 
		bool writeBrushFromBBox(int iModel, myVtx_t origin);
		// 
		bool writeBoxBrush(myVtx_t bboxMin, myVtx_t bboxMax, const BrushAttr &brushAttr);
		// 
		bool summarize();
		// 
		bool writeLevelOfDetailedTerrain();
		// 
		bool formatTexfc(const myPoly_t &poly, UINT iTex, Texfc &t);
	};
};

/* /////////////////////////////////////////////////////////

///////////////////////////////////////////////////////// */

#define MB2MML_F_MAKE_BMP	(0x00000001)
#define MB2MML_F_ENABLE_BBOX_TREATS	(0x00000002)
#define MB2MML_F_CALC_TEXC_ACTUALLY	(0x00000004)

struct MB2MMLTexMetr
{
	// 
	CString tstrName;
	// 
	UINT cx, cy;
};

struct MB2MMLOptical
{
	// 
	UINT nMask;
	// 
	SIZE sizeMakeBM;
	// 
	HBITMAP hMadeBM;
	// 
	CStringArray arrBBoxTreats;
	// 
	CArray<MB2MMLTexMetr, MB2MMLTexMetr> arrTexMetr;
	// 
	MB2MMLOptical(): nMask(0) { }

};

/* /////////////////////////////////////////////////////////
  óAèo
///////////////////////////////////////////////////////// */

extern int _cdecl _MohBSP2Map(LPCTSTR lpszIn1, LPCTSTR lpszInto1);
extern int _cdecl _MohBSP2Map(LPCTSTR lpszIn1, LPCTSTR lpszInto1, MB2MMLOptical &rOpts);
extern int g_bVerbose;

/* /////////////////////////////////////////////////////////
  Ç≤Ç…ÇÂÇ≤Ç…ÇÂ
///////////////////////////////////////////////////////// */
