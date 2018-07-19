
#pragma warning(disable: 4786)

#include <afxwin.h>
#include <afxole.h>
#include <atlconv.h>
#include <math.h>
#include <float.h>
#include <afxtempl.h>
#include <limits.h>
#include "MohBSP2MAP.h"

#include "CSVrw.h"
#include "OSP.h"
#include "Texc.h"

using namespace OSP;

/* /////////////////////////////////////////////////////////
  ::MB2MM
///////////////////////////////////////////////////////// */

bool MB2MM::ReadInt(LPCTSTR psz, int &x)
{
	LPTSTR pszEnd = NULL;
	x = _tcstol(psz, &pszEnd, 10);
	return (*pszEnd == 0);
}

/* /////////////////////////////////////////////////////////
  ::MB2MM::Knowledge
///////////////////////////////////////////////////////// */

bool MB2MM::Knowledge::LoadTexMetrFiles(CString strDir)
{
	CString strIn;
	for (int i = 0; i < 99; i++) {
		strIn.Format("texmetr%d.csv", i);
		CSVR_t csvr;
		if (!csvr.Open(OSP_JoinPath(strDir, strIn))) break;
		CString str[3];
		while (csvr.ReadNextLine()) {
			TextureMetric tm;
			if (false
				|| !csvr.ReadNextToken(str[0])
				|| !csvr.ReadNextToken(str[1])
				|| !ReadInt(str[1], (int &)tm.cx)
				|| !csvr.ReadNextToken(str[2])
				|| !ReadInt(str[2], (int &)tm.cy)
			) {
				continue;
			}
			tm.strName = str[0];
			m[str[0]] = tm;
		}
	}
	return true;
}

bool MB2MM::Knowledge::LoadBBoxFiles(CString strDir)
{
	CString strIn;
	for (int i = 0; i < 99; i++) {
		strIn.Format("bbox%d.csv", i);
		CSVR_t csvr;
		if (!csvr.Open(OSP_JoinPath(strDir, strIn))) break;
		CString str[1];
		while (csvr.ReadNextLine()) {
			if (false
				|| !csvr.ReadNextToken(str[0])
				|| str[0].IsEmpty()
			) {
				continue;
			}
			triggers.insert(str[0]);
		}
	}
	return true;
}

/* /////////////////////////////////////////////////////////
  ::MB2MM::myVtx_t
///////////////////////////////////////////////////////// */

MB2MM::n2fType MB2MM::myVtx_t::Normal2Facing(const myVtx_t &v)
{
	ASSERT(-1 <= v.x && v.x <= +1);
	ASSERT(-1 <= v.y && v.y <= +1);
	ASSERT(-1 <= v.z && v.z <= +1);
	if (v.x > +EPS && v.y > +EPS && fabs(+v.x - +v.y) < EPS) return n2fPYZ;
	if (v.x < -EPS && v.y > +EPS && fabs(-v.x - +v.y) < EPS) return n2fNXZ;
	if (v.x > +EPS && v.y < -EPS && fabs(+v.x - -v.y) < EPS) return n2fPXZ;
	if (v.x < -EPS && v.y < -EPS && fabs(-v.x - -v.y) < EPS) return n2fNYZ;
	if (fabs(v.x) <= fabs(v.y)) {
		// North/South
		if (fabs(v.y) <= fabs(v.z))
			// Upper/Down
			if (0 <= v.z)
				// Upper
				return n2fPXY;
			else
				// Down
				return n2fNXY;
		if (0 <= v.y)
			// North
			return n2fNXZ;
		else
			// South
			return n2fPXZ;
	} else {
		// West/East
		if (fabs(v.x) <= fabs(v.z))
			// Upper/Down
			if (0 <= v.z)
				// Upper
				return n2fPXY;
			else
				// Down
				return n2fNXY;
		if (0 <= v.x)
			// East
			return n2fPYZ;
		else
			// West
			return n2fNYZ;
	}
}

/* /////////////////////////////////////////////////////////
  ::MB2MM::myPolySetList_t
///////////////////////////////////////////////////////// */

bool MB2MM::myPolySetList_t::JoinTri(const myTri_t &t, const myVtx_t &norm, UINT texture, bool fReduce)
{
	myPolyKey_t k;
	k.texture = texture;
	k.norm = norm;
	iterator
		iterPos = lower_bound(k),
		iterEnd = upper_bound(k);

	if (fReduce) {
		for (; iterPos != iterEnd; iterPos++) {
			myPolySet_t &polyset = iterPos->second;
			if (polyset.z.JoinTri(t)) {
				ASSERT(polyset.z.IsConvexShape());
				return true;
			}
		}
	}
	iterPos = insert(value_type(k, myPolySet_t()));
	{
		myPolySet_t &polyset = iterPos->second;
		polyset.z.JoinTri(t);
		polyset.z.vtxNorm = -norm;
		polyset.z.arrTris.Add(t);
	}
	return true;
}

void MB2MM::myPolySetList_t::JoinPolys()
{
	set<myPolyKey_t> polykeyset;

	myPolySetList_t::iterator iterEnd = this->end();
	
	{
		myPolySetList_t::iterator iterPos = this->begin();

		for (; iterPos != iterEnd; iterPos = this->upper_bound(iterPos->first)) {
			VERIFY(static_cast<bool>(polykeyset.insert(iterPos->first).second));
		}
	}

	{
		while (true) {
			int nRepair = 0;

			set<myPolyKey_t>::iterator
				iterUPos = polykeyset.begin(),
				iterUEnd = polykeyset.end();

			int nKeys = distance(iterUPos, iterUEnd);

			for (; iterUPos != iterUEnd; ) {

				myPolySetList_t::iterator
					iter1Pos = this->lower_bound(*iterUPos),
					iter1End = this->upper_bound(*iterUPos);

				int nLen = distance(iter1Pos, iter1End);

				for (; iter1Pos != iter1End; iter1Pos++) {

					myPolySetList_t::iterator
						iter2Pos = iter1Pos,
						iter2End = iter1End;
					for (iter2Pos++; iter2Pos != iter2End; ) {
						if (iter1Pos->second.z.JoinPolySafely(iter2Pos->second.z)) {
							iter2Pos = this->erase(iter2Pos);
							nRepair++;
						} else {
							iter2Pos++;
						}
					}
				}

				iterUPos++;
			}

			if (nRepair == 0) break;
		}
	}
}

void MB2MM::myPolySetList_t::ReduceNeighbor()
{
	myPolySetList_t::iterator
		iter1Pos = begin(),
		iter1End = end();
	for (; iter1Pos != iter1End; ) {
		if (iter1Pos->second.z.ReduceNeighbor()) {
			iter1Pos++;
		} else {
			iter1Pos = erase(iter1Pos);
		}
	}
}

void MB2MM::myPolySetList_t::OnDevFinalReview()
{
#ifdef _DEBUG
	{
		UINT i = 0;
		myPolySetList_t::iterator
			iter1Pos = begin(),
			iter1End = end();
		for (; iter1Pos != iter1End; iter1Pos++, i++) {
			iter1Pos->second.z.OnDevFinalReview(i);
		}
	}
#endif
}

/* /////////////////////////////////////////////////////////
  ::MB2MM::Decompiler
///////////////////////////////////////////////////////// */

bool MB2MM::Decompiler::FindTexture(UINT i, TextureMetric &tm)
{
	myTexture_t *pTex = GetFineTex(i);
	if (pTex == NULL) return false;

	TextureMetricMap::iterator
		iterPos = kb.m.find(pTex->strName),
		iterEnd = kb.m.end();
	if (iterPos == iterEnd) return false;

	tm = iterPos->second;
	return true;
}

void MB2MM::Decompiler::EnsureTexture(UINT iTex, TextureMetric &tm)
{
	tm.strName = "common/black";
	tm.cx = 128;
	tm.cy = 128;

	if (!FindTexture(iTex, tm)) {
		myTexture_t *pTexture = GetFineTex(iTex);
		if (pTexture != NULL) {
			tm.strName = pTexture->strName;
		}
	}
}

bool MB2MM::Decompiler::Decompile(CString strIn, CString strOut)
{
	if (strIn.IsEmpty() || strOut.IsEmpty()) return false;

	ReadBSP rb;
	if (false
		|| !rb.Open(strIn)
		|| !rb.ReadEntry( 0, sbTexture)
		|| !rb.ReadEntry( 3, sbSurface)
		|| !rb.ReadEntry( 4, sbVertex)
		|| !rb.ReadEntry( 5, sbMesh)
		|| !rb.ReadEntry(13, sbModel)
		|| !rb.ReadEntry(14, sbEntity)
		|| !rb.ReadEntry(22, sbLOD)
		|| !entityset.Read(sbEntity)
	) {
		return false;
	}
	rb.Close();

	pModels = sbModel.GetData(), nModels = sbModel.GetSize() / 40;
	pSurfaces = sbSurface.GetData(), nSurfaces = sbSurface.GetSize() / 108;
	pVertices = sbVertex.GetData(), nVertices = sbVertex.GetSize() / 44;
	pTextures = sbTexture.GetData(), nTextures = sbTexture.GetSize() / 140;
	pMeshes = sbMesh.GetData(), nMeshes = sbMesh.GetSize() / 4;
	pLODs = sbLOD.GetData(), nLODs = sbLOD.GetSize() / 388;

	if (!summarize())
		return false;

	if (!wm.Open(strOut))
		return false;

	UINT i;
	for (i = 0; i < entityset.entities.GetSize(); i++) {
		Entity1 &r = entityset.entities[i];
		mo.fProcessingWorld = (i == 0);
		mo.AddEntityClass(r.className);
		mo.AddModelName(r.modelName);
		if (!wm.WriteCommentNo("entity", i)) return false;
		if (!wm.WriteOpenBracket()) return false;
		//
		Entity1::iterator
			iterPos = r.begin(),
			iterEnd = r.end();
		for (; iterPos != iterEnd; iterPos++) {
			if (!wm.WriteAttr(iterPos->first, iterPos->second)) return false;
		}
		//
		if (!writeModel(i, r.modelNo)) return false;
		//
		if (!wm.WriteCloseBracket()) return false;
	}
	wm.Close();
	return true;
}

// Unofficial Quake 3 BSP Format
// http://www.gametutorials.com/Tutorials/OpenGL/Quake3Format.htm

// Unofficial Quake 3 Map Specs
// http://graphics.stanford.edu/~kekoa/q3/

// number of brushes, number of entities, file size, number of duplicate planes, etc... 

bool MB2MM::Decompiler::writeModel(UINT entityNo, int modelNo)
{
	if (entityNo == 0)
		modelNo = 0;
	if (modelNo < 0)
		return true;
	const model3_t *pModel = GetModel(modelNo);
	ASSERT(pModel != NULL);
	if (pModel == NULL) return false;

	const Entity1 &entity = entityset.entities[entityNo];

	brushNo = 0;

	if (kb.IsTriggerEntity(entity.className)) {
		if (!writeBrushFromBBox(modelNo, entity.v_origin)) return false;
		return true;
	}

#ifdef _DEBUG_ONLY_WORLD_SPAWN
	if (entityNo != 0)
		return true;
#endif

	PatchDef2SetList patchDef2sets;

	const bool fUniqueTris = optz.fUniqueTris;
	const bool fRemoveSharp = optz.fRemoveSharp;
	const bool fRecoverTexfc = optz.fRecoverTexfc;
	const bool fSkew = optz.fSkew;

	// Surfaces to poly set list
	{
		MatchTriList match3List;

		UINT iSurface;
		const UINT iSurfaceOffset = pModel->face;
		const UINT nSurfaces = pModel->n_faces;
		for (iSurface = 0; iSurface < nSurfaces; iSurface++) {
			surface3_t *pSurface = GetSurface(iSurfaceOffset + iSurface);
			ASSERT(pSurface != NULL);
			if (pSurface == NULL) return false;

			if (pSurface->type == 1) {
				const UINT iVertOffset = pSurface->vertex;
				const UINT nVerts = pSurface->n_vertexes;
				const UINT iMeshVertOffset = pSurface->meshvert;
				const UINT nMeshVerts = pSurface->n_meshverts;
				const UINT texture = (nTextures <= pSurface->texture) ? -1 : pSurface->texture;
				ASSERT((nMeshVerts % 3) == 0);
				myVtx_t norm(pSurface->normal[0], pSurface->normal[1], pSurface->normal[2]);
				UINT iTri;
				const UINT nTris = nMeshVerts / 3;
				for (iTri = 0; iTri < nTris; iTri++) {
					MatchTri mt;
					myTri_t tri;
					BYTE iVert;
					for (iVert = 0; iVert < 3; iVert++) {
						DWORD *pMeshVert = GetMesh(iMeshVertOffset + 3 * iTri + iVert);
						if (pMeshVert == NULL) return false;
						vertex3_t *pVertex = GetVertex(iVertOffset + *pMeshVert);
						if (pVertex == NULL) return false;
						tri.v[iVert].x = pVertex->position[0];
						tri.v[iVert].y = pVertex->position[1];
						tri.v[iVert].z = pVertex->position[2];
						tri.v[iVert].tu = pVertex->texcoord[0][0];
						tri.v[iVert].tv = pVertex->texcoord[0][1];
						tri.v[iVert] += entity.v_origin;
					}

#ifdef _DEBUG_TRACE_SURFACE_TRIS
					{
						TRACE0("==============================\n");
						for (char i = 0; i < 3; i++) {
							TRACE3("3: %5d,%5d,%5d\n"
								, (int)tri.v[i].x
								, (int)tri.v[i].y
								, (int)tri.v[i].z
								);
							Sleep(1);
						}
						TRACE0("---\n");
						TRACE3("n: %5.2f,%5.2f,%5.2f\n"
							, norm.x
							, norm.y
							, norm.z
							);
						Sleep(1);
						TRACE0("---\n");
					}
#endif
					mt.v[0] = tri.v[0];
					mt.v[1] = tri.v[1];
					mt.v[2] = tri.v[2];

					if (!tri.IsValidTri()) {
						mo.AddBadTri();
						continue;
					}
					if (fUniqueTris && !match3List.Add(mt)) {
						mo.AddJunkTri();
						continue;
					}

					//tri.Norm(norm);
					VERIFY(polysets.JoinTri(tri, norm, texture, 1 <= optz.nReduction));

					mo.AddTri();
				}
			} else if (pSurface->type == 2) {
				const UINT cx = pSurface->size[1];
				const UINT cy = pSurface->size[0];

				const UINT texture = (nTextures <= pSurface->texture) ? -1 : pSurface->texture;

				if (cx == 0 || cy == 0)
					return false;

				patchDef2sets.push_back(PatchDef2Set());

				PatchDef2Set &patchdef2s = patchDef2sets.back();

				patchdef2s.a.cx = cx;
				patchdef2s.a.cy = cy;

				if (texture < nTextures) {
					patchdef2s.a.SetTexture(fineTextures[texture]);
				} else {
					patchdef2s.a.UndefTexture();
				}

				const UINT iVertOffset = pSurface->vertex;

				UINT y;
				for (y = 0; y < cy; y++) {
					UINT x;
					for (x = 0; x < cx; x++) {
						vertex3_t *pVertex = GetVertex(iVertOffset +y +cy*x);
						if (pVertex == NULL) return false;

						PatchDef2Vtx v;
						v.v.x = pVertex->position[0];
						v.v.y = pVertex->position[1];
						v.v.z = pVertex->position[2];
						v.v.tu = pVertex->texcoord[0][0];
						v.v.tv = pVertex->texcoord[0][1];

						patchdef2s.z.z.Add(v);
					}
				}
			}
		}
		printf(""); // match3List
	}

	if (entityNo == 0) {
		if (!writeLevelOfDetailedTerrain()) return false;
	}

	polysets.OnDevFinalReview();

	if (2 <= optz.nReduction) polysets.JoinPolys();

	polysets.ReduceNeighbor();

	{
		PatchDef2SetList::iterator
			iterPos = patchDef2sets.begin(),
			iterEnd = patchDef2sets.end();
		for (; iterPos != iterEnd; iterPos = patchDef2sets.erase(iterPos)) {
			const PatchDef2Set &patchDef2s = *iterPos;

			if (!wm.WriteCommentNo("brush", brushNo)) return false;
			if (!wm.WriteOpenBracket()) return false;
			if (!wm.WritePatchDef2()) return false;
			if (!wm.WriteOpenBracket()) return false;
			if (!wm.WritePatchDef2Attr(patchDef2s.a)) return false;
			if (!wm.WriteOpenMBracket()) return false;
			//
			const UINT cx = patchDef2s.a.cx;
			const UINT cy = patchDef2s.a.cy;
			UINT y;
			for (y = 0; y < cy; y++) {
				if (!wm.WriteMBracketBegin()) return false;
				//
				UINT x;
				for (x = 0; x < cx; x++) {
					if (!wm.WritePatchDef2Vertex(patchDef2s.z.z[x + cx*y])) return false;
				}
				//
				if (!wm.WriteMBracketEnd()) return false;
			}
			//
			if (!wm.WriteCloseMBracket()) return false;
			if (!wm.WriteCloseBracket()) return false;
			if (!wm.WriteCloseBracket()) return false;
			brushNo++;
			mo.AddPMesh();
		}
	}

	{
		myPolySetList_t::iterator
			iterPos = polysets.begin(),
			iterEnd = polysets.end();
		for (; iterPos != iterEnd; iterPos = polysets.erase(iterPos)) {
			const myPolyKey_t &polykey = iterPos->first;
			myPolySet_t &polyset = iterPos->second;

			//if (!(brushNo < 1)) continue;

			//polyset.z.Reduce();

			if (polyset.z.arrVerts.GetSize() < 3)
				continue;
			if (fRemoveSharp && polyset.z.IsNeedleBrush())
				continue;

			CArray<myVtx2_t, myVtx2_t> &verts = polyset.z.arrVerts;

			const UINT nVerts = verts.GetSize();
			const UINT iTex = polykey.texture;

			myVtx_t backwardNorm = polykey.norm * -optz.nThick;

			BrushTri t;

			if (iTex < nTextures) {
				t.SetTexture(fineTextures[iTex]);
			} else {
				t.UndefTexture();
			}

			if (!wm.WriteCommentNo("brush", brushNo)) return false;
			if (!wm.WriteOpenBracket()) return false;

			CArray<myVtx_t, myVtx_t> vertsBackwardVec;

			{
				UINT i;
				const UINT n = verts.GetSize();
				vertsBackwardVec.InsertAt(0, backwardNorm, n);

				if (fSkew) {
					for (i = 0; i < n; i++) {
						const myVtx_t &v0 = verts[(n + i - 1) % n];
						const myVtx_t &v1 = verts[(    i    ) % n];
						const myVtx_t &v2 = verts[(    i + 1) % n];

						myVtx_t vv1 = myVtx_t::CrossProduct(polykey.norm, v2 - v1);
						myVtx_t vv2 = myVtx_t::CrossProduct(polykey.norm, v1 - v0);

						myVtx_t vvx = vv1 + vv2;
						vvx.Normalize();
						vvx = vvx * -optz.nThick;

						vertsBackwardVec[i] = backwardNorm + vvx;
					}
				}
			}

			// Front face
			{
				if (fRecoverTexfc) {
					formatTexfc(polyset.z, iTex, t);
				}

				BYTE v;
				for (v = 0; v < 3; v++) {
					t.tri.v[v] = verts[v];
				}
				if (!wm.WriteBrush(t)) return false;
			}

			t.UndefTexture(false);

			// Back face
			{
				BYTE v;
				for (v = 0; v < 3; v++) {
					t.tri.v[v] = myVtx2_t(verts[2 - v] + vertsBackwardVec[2 - v]);
				}
				if (!wm.WriteBrush(t)) return false;
			}
			// Side faces
			{
				UINT i;
				for (i = 0; i < nVerts; i++) {
					t.tri.v[0] = verts[(i+0)       ];
					t.tri.v[2] = verts[(i+1)%nVerts];
					t.tri.v[1] = myVtx2_t(t.tri.v[2] + vertsBackwardVec[(i+1)%nVerts]);

					if (!wm.WriteBrush(t)) return false;
				}
			}
			if (!wm.WriteCloseBracket()) return false;
			brushNo++;
			mo.AddBrush();
		}
	}
	//
	return true;
}

bool MB2MM::Decompiler::writeBrushFromBBox(int modelNo, myVtx_t origin)
{
	const model3_t *pModel = GetModel(modelNo);
	ASSERT(pModel != NULL);

	myVtx_t bboxMin(pModel->mins[0], pModel->mins[1], pModel->mins[2]);
	myVtx_t bboxMax(pModel->maxs[0], pModel->maxs[1], pModel->maxs[2]);

	BrushAttr b;
	b.strName = "common/trigger";

	return writeBoxBrush(bboxMin + origin, bboxMax + origin, b);
}

bool MB2MM::Decompiler::writeBoxBrush(myVtx_t bboxMin, myVtx_t bboxMax, const BrushAttr &brushAttr)
{
	static const char iBoxTbl[6][3] = {
		{4,6,5}, // Upper
		{1,3,0}, // Lower
		{1,5,3}, // East
		{2,6,0}, // West
		{0,4,1}, // South
		{3,7,2}, // North
	};

	BrushTri brushNew = brushAttr;

	if (!wm.WriteCommentNo("brush", brushNo)) return false;
	if (!wm.WriteOpenBracket()) return false;
	//
	for (int tri = 0; tri < 6; tri++) {
		for (int v = 0; v < 3; v++) {
			int f = iBoxTbl[tri][v];
			brushNew.tri.v[v].x = (f & 1) ? bboxMin.x : bboxMax.x;
			brushNew.tri.v[v].y = (f & 2) ? bboxMax.y : bboxMin.y;
			brushNew.tri.v[v].z = (f & 4) ? bboxMin.z : bboxMax.z;
		}
		if (!wm.WriteBrush(brushNew)) return false;
	}
	//
	if (!wm.WriteCloseBracket()) return false;
	return true;
}

bool MB2MM::Decompiler::summarize()
{
	UINT i;

	fineTextures.SetSize(nTextures);

	for (i = 0; i < nTextures; i++) {
		const tex3_t *pTex = GetTex(i);
		ASSERT(pTex != NULL);
		myTexture_t &tex = fineTextures[i];

		tex.strName = pTex->name;
		tex.surf = pTex->surf;
		tex.contents = pTex->contents;
	}

	return true;
}

bool MB2MM::Decompiler::writeLevelOfDetailedTerrain()
{
	UINT i;
	for (i = 0; i < nLODs; i++) {
		terrain_t *pGEOMetric = GetLODt(i);
		ASSERT(pGEOMetric != NULL);

		UINT texture = pGEOMetric->texture;

		LODTerrainHead header;
		header.cx = 9;
		header.cy = 9;
		header.origin.x = 64 * pGEOMetric->x;
		header.origin.y = 64 * pGEOMetric->y;
		header.origin.z = pGEOMetric->z;
		header.f = (pGEOMetric->flags & 0x40) ? 1 : 0;

		LODTerrainPart partition;

		myTexture_t *pTexture = GetFineTex(texture);

		if (pTexture != NULL) {
			partition.SetTexture(*pTexture);
		} else {
			partition.UndefTexture();
		}

		if (false
			|| !wm.WriteCommentNo("brush", brushNo)
			|| !wm.WriteOpenBracket()
			|| !wm.WriteTerrainDef()
			||  !wm.WriteOpenBracket()
			||  !wm.WriteTerrainDef2Head(header)
			||   !wm.WriteOpenBracket()
			||    !wm.WriteTerrainDef2Part(partition)
			||    !wm.WriteTerrainDef2Part(partition)
			||    !wm.WriteTerrainDef2Part(partition)
			||    !wm.WriteTerrainDef2Part(partition)
			||   !wm.WriteCloseBracket()
			||   !wm.WriteOpenBracket()
		) {
			return false;
		}
		//
		LODTerrainUnit t;
		char y;
		for (y = 0; y < 9; y++) {
			char x;
			for (x = 0; x < 9; x++) {
				t.z = 2 * pGEOMetric->m[y][x];

				if (!wm.WriteTerrainDef2Unit(t)) return false;
			}
		}
		//
		if (false
			||   !wm.WriteCloseBracket()
			||  !wm.WriteCloseBracket()
			|| !wm.WriteCloseBracket()
		) {
			return false;
		}

		mo.AddLODt();
	}

	return true;
}

bool MB2MM::Decompiler::formatTexfc(const myPoly_t &poly, UINT iTex, Texfc &t)
{
	ASSERT(3 <= poly.arrVerts.GetSize());

	TextureMetric tm;
	EnsureTexture(iTex, tm);

	if (tm.cx == 0 || tm.cy == 0) return false;

	Texc texc;

	texc.cxTex = tm.cx;
	texc.cyTex = tm.cy;

	n2fType n2fPoly = myVtx_t::Normal2Facing(poly.vtxNorm);

	const UINT nTrys = poly.arrTris.GetSize();

	UINT iTry;
	for (iTry = 0; iTry < nTrys; iTry++) {
		myTri_t s = poly.arrTris[iTry];
		switch (n2fPoly) {
		case n2fPXY:
		case n2fNXY:
			{
				for (char i = 0; i < 3; i++) texc.vv[i][0] = s.v[i].x, texc.vv[i][1] = s.v[i].y;
				break;
			}
		case n2fPYZ:
		case n2fNYZ:
			{
				for (char i = 0; i < 3; i++) texc.vv[i][0] = s.v[i].y, texc.vv[i][1] = s.v[i].z;
				break;
			}
		case n2fPXZ:
		case n2fNXZ:
			{
				for (char i = 0; i < 3; i++) texc.vv[i][0] = s.v[i].x, texc.vv[i][1] = s.v[i].z;
				break;
			}
		}

		for (char i = 0; i < 3; i++) texc.tv[i][0] = s.v[i].tu, texc.tv[i][1] = s.v[i].tv;

		switch (n2fPoly) {
		case n2fNXY:
		case n2fNYZ:
		case n2fNXZ:
			texc.fMirror = true;
			break;
		case n2fPXY:
		case n2fPYZ:
		case n2fPXZ:
			texc.fMirror = false;
			break;
		}

		if (texc.Decode()) {
			t.tx = texc.tx;
			t.ty = texc.ty;
			t.angle = texc.angle;
			t.tsx = texc.tsx;
			t.tsy = texc.tsy;
			return true;
		}
	}
	return false;
}

/* /////////////////////////////////////////////////////////
  OBSOLETE
///////////////////////////////////////////////////////// */

#if 0 // ->->->->->->->->->->->->->->->->->->->->

using namespace MB2MM;

typedef MB2MM::surface3_t face3_t;

void calcOPVec(myVtx_t &vOP, const myVtx_t &v1, const myVtx_t &v2)
{
	vOP.x = v1.y * v2.z - v1.z * v2.y;
	vOP.y = v1.z * v2.x - v1.x * v2.z;
	vOP.z = v1.x * v2.y - v1.y * v2.x;
}

void MB2MM::calcN1Vec(myVtx_t &v)
{
	double dXY = sqrt(v.x * v.x + v.y * v.y);
	double dNZ = sqrt(dXY * dXY + v.z * v.z);
	if (dNZ == 0)
		return;
	double dRev = 1.0 / dNZ;
	v.x *= dRev;
	v.y *= dRev;
	v.z *= dRev;
}

void calcNormVec(myVtx_t &vNorm, const myVtx_t &vA, const myVtx_t &vB, const myVtx_t &vC)
{
	myVtx_t vT1(
		vB.x - vA.x,
		vB.y - vA.y,
		vB.z - vA.z
		);
	myVtx_t vT2(
		vC.x - vA.x,
		vC.y - vA.y,
		vC.z - vA.z
		);
	calcOPVec(vNorm, vT2, vT1);
	calcN1Vec(vNorm);
}

#if 0
//	int calcNVoS(const myVtx_t &v)
//	{
//		ASSERT(-1 <= v.x && v.x <= 1);
//		double fRz = acos(v.x) / 3.1415926535897932384626433832795 * 180.0;
//		if (v.y < 0)
//			fRz = 360 - fRz;
//		double fRx = asin(v.z) / 3.1415926535897932384626433832795 * 180.0;
//		int iZ = -1, iX = -1;
//		double fPx =(fRx < 0) ? -fRx : fRx;
//		if      ( 45.0 < fRz && fRz <= 135.0)
//			iZ = 1;
//		else if (135.0 < fRz && fRz <= 225.0)
//			iZ = 2;
//		else if (225.0 < fRz && fRz <= 315.0)
//			iZ = 3;
//		else if (315.0 < fRz || fRz <=  45.0)
//			iZ = 0;
//		if      ( 45.0 < fPx && fPx <= 135.0)
//			iX = 1;
//		else if (135.0 < fPx && fPx <= 225.0)
//			iX = 2;
//		else if (225.0 < fPx && fPx <= 315.0)
//			iX = 3;
//		else if (315.0 < fPx || fPx <=  45.0)
//			iX = 0;
//		ASSERT(iX != -1 && iZ != -1);
//		ASSERT(TRUE);
//		static const int x = -1;
//		static const int iNVoSMtx[4][4] = {
//			{1,3,4,0},//{1,0,1,0},
//			{2,2,2,2},//{2,2,2,2},
//			{x,x,x,x},//{x,x,x,x},
//			{5,5,5,5},//{2,2,2,2},
//		};
//		int iNVoS = iNVoSMtx[iX][iZ];
//		ASSERT(iNVoS == 0 || iNVoS == 1 || iNVoS == 2 || iNVoS == 3 || iNVoS == 4 || iNVoS == 5);
//		return iNVoS;
//	}
#endif

/* /////////////////////////////////////////////////////////
  ごにょごにょ
///////////////////////////////////////////////////////// */

namespace
{

class CFindTexMetr
{
public:
	// 
	virtual ~CFindTexMetr() { }
	// 
	virtual BOOL Find(LPCTSTR lpszName, UINT &cx, UINT &cy) = 0;

};

class CFindTexMetrExp : public CFindTexMetr
{
public:
	// 
	typedef CArray<MB2MMLTexMetr, MB2MMLTexMetr> MB2MMLTexMetrArray;
	// 
	CFindTexMetrExp(const MB2MMLTexMetrArray &r): m_rArr(r) { }
	// 
	BOOL Find(LPCTSTR lpszName, UINT &cx, UINT &cy)
	{
		for (int i = 0; i < m_rArr.GetSize(); i++) {
			if (m_rArr[i].tstrName != lpszName)
				continue;
			cx = m_rArr[i].cx;
			cy = m_rArr[i].cy;
			return TRUE;
		}
		return FALSE;
	}

private:
	// 
	const MB2MMLTexMetrArray &m_rArr;

};

class myTerrain_t
{
public:
	// 
	UINT nBit[9][9];
	// 
	myVtx_t origin;
	// 
	UINT iTex;

};

#define NVOS_PXZ 0
#define NVOS_PYZ 1
#define NVOS_PXY 2
#define NVOS_NXZ 3
#define NVOS_NYZ 4
#define NVOS_NXY 5

inline double calcRADFromDEG(double a) { return a / 180.0 * 3.14159265358979; }
inline double calcDEGFromRAD(double a) { return a / 3.14159265358979 * 180.0; }

inline double calcDEGFrom0To360(double a) { while ((int)floor(a + 0.5) < 0) a += 360; while ((int)floor(a + 0.5) >= 360) a -= 360; return a; }

int calcNVoS(const myVtx_t &v)
{
	ASSERT(-1 <= v.x && v.x <= +1);
	ASSERT(-1 <= v.y && v.y <= +1);
	ASSERT(-1 <= v.z && v.z <= +1);
	if (v.x > +1E-5 && v.y > +1E-5 && fabs(+v.x - +v.y) < 1E-5) return NVOS_PYZ;
	if (v.x < -1E-5 && v.y > +1E-5 && fabs(-v.x - +v.y) < 1E-5) return NVOS_NXZ;
	if (v.x > +1E-5 && v.y < -1E-5 && fabs(+v.x - -v.y) < 1E-5) return NVOS_PXZ;
	if (v.x < -1E-5 && v.y < -1E-5 && fabs(-v.x - -v.y) < 1E-5) return NVOS_NYZ;
	if (fabs(v.x) <= fabs(v.y)) {
		// 北向き・南向き
		if (fabs(v.y) <= fabs(v.z))
			// 上向き・下向き
			if (0 <= v.z)
				// 上向き
				return NVOS_PXY;
			else
				// 下向き
				return NVOS_NXY;
		if (0 <= v.y)
			// 北向き
			return NVOS_NXZ;
		else
			// 南向き
			return NVOS_PXZ;
	} else {
		// 西向き・東向き
		if (fabs(v.x) <= fabs(v.z))
			// 上向き・下向き
			if (0 <= v.z)
				// 上向き
				return NVOS_PXY;
			else
				// 下向き
				return NVOS_NXY;
		if (0 <= v.x)
			// 東向き
			return NVOS_PYZ;
		else
			// 西向き
			return NVOS_NYZ;
	}
}

double clampVal(double d, int n)
{
	int x = n;
	for (; x > 0; x--)
		d *= 10.0;
	if (d < 0)
		d = ceil(d - 0.49);
	else
		d = floor(d + 0.49);
	for (; n > 0; n--)
		d /= 10.0;
	return d;
}

double fix1P0(double d)
{
	if (memcmp(&d, "\xFF\xFF\xFF\xFF\xFF\xFF\xEF\x3F", 8) == 0)
		memcpy(&d, "\x00\x00\x00\x00\x00\x00\xF0\x3F", 8);
	if (memcmp(&d, "\xFF\xFF\xFF\xFF\xFF\xFF\xEF\xBF", 8) == 0)
		memcpy(&d, "\x00\x00\x00\x00\x00\x00\xF0\xBF", 8);
	return d;
}

//	void calcRevAtMed(double &f0, double &f1, double &f2)
//	{
//		double fL = __min(f0, __min(f1, f2));
//		double fR = __max(f0, __min(f1, f2));
//		double fM = (fL + fR) / 2.0;
//		f0 = fM * 2 - f0;
//		f1 = fM * 2 - f1;
//		f2 = fM * 2 - f2;
//	}

//	void rotateXYc(double fX, double fY, double bX, double bY, double &fXr, double &fYr, double r)
//	{
//		r = r / 180.0 * 3.1415926535897932384626433832795;
//		double m[2][2] = {
//			{ cos(r), sin(r)},
//			{-sin(r), cos(r)}
//		};
//		fX -= bX;
//		fY -= bY;
//		fXr = fX * m[0][0] + fY * m[1][0];
//		fYr = fX * m[0][1] + fY * m[1][1];
//		fXr += bX;
//		fYr += bY;
//	}

//	double calc360RevAngle(double f)
//	{
//		while (f < 0)
//			f += 360;
//		while (f >=360.0)
//			f -= 360;
//		if (  0  < f && f < 180)
//			return 360 - f;
//		if (180 <= f && f < 360)
//			return 360 - f;
//		return f;
//	}

double calcTexc_AngleXY(
	double x, double y
	)
{
	if (fabs(x) < 1E-5) {
		if (y < 0)
			//270°
			return -3.1415926535897932384626433832795 * 0.5;
		else
			// 90°
			return +3.1415926535897932384626433832795 * 0.5;
	} else {
		if (x > 0)
			return atan(y / x);
		else
			return atan(y / x) + 3.1415926535897932384626433832795;
	}
}

double calcTexc_AngleXYd(
	double x, double y
	)
{
	double a = calcDEGFrom0To360(calcDEGFromRAD(calcTexc_AngleXY(x, y)));
	return a;
}

void calcTexc_rotateXYd(double fX, double fY, double &fXr, double &fYr, double f)
{
	f = calcRADFromDEG(f);
	double eM11 = +cos(f);
	double eM12 = +sin(f);
	double eM21 = -eM12;
	double eM22 = +eM11;
	fXr = fX * eM11 + fY * eM21;
	fYr = fX * eM12 + fY * eM22;
}

double calcTexc_Length(
	double x, double y
	)
{
	double z = sqrt(x * x + y * y);
	return z;
}

char calcTexc_IdentifiedSign(double x)
{
	if (x > +1E-5) return +1;
	if (x < -1E-5) return -1;
	return 0;
}

void calcTexc_rotateTriPos(
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

BOOL calcTexc_TouchToXY(
	double fTx0, double fTy0,
	double fTx1, double fTy1,
	double fTx2, double fTy2,
	double fAx0, double fAy0,
	double fAx1, double fAy1,
	double fAx2, double fAy2,
	UINT cxTex, UINT cyTex, UINT iMirrorTex,
	double &fTexsX, double &fTexsY, double &fTexAngle, double &fTexcX, double &fTexcY
	)
{
	if (iMirrorTex) {
		double x;
		x = fAx1; fAx1 = fAx2; fAx2 = x;
		x = fTx1; fTx1 = fTx2; fTx2 = x;
	}
	double fTs[][3][2] = {
		fTx0, fTy0, fTx1, fTy1, fTx2, fTy2,
		fTx0, fTy0, fTx2, fTy2, fTx1, fTy1,
		fTx1, fTy1, fTx2, fTy2, fTx0, fTy0,
		fTx1, fTy1, fTx0, fTy0, fTx2, fTy2,
		fTx2, fTy2, fTx0, fTy0, fTx1, fTy1,
		fTx2, fTy2, fTx1, fTy1, fTx0, fTy0,
	};
	double fAs[][3][2] = {
		fAx0, fAy0, fAx1, fAy1, fAx2, fAy2,
		fAx0, fAy0, fAx2, fAy2, fAx1, fAy1,
		fAx1, fAy1, fAx2, fAy2, fAx0, fAy0,
		fAx1, fAy1, fAx0, fAy0, fAx2, fAy2,
		fAx2, fAy2, fAx0, fAy0, fAx1, fAy1,
		fAx2, fAy2, fAx1, fAy1, fAx0, fAy0,
	};
	for (int iXc = 0; iXc < 6; iXc++) {
		fTx0 = fTs[iXc][0][0], fTy0 = fTs[iXc][0][1];
		fTx1 = fTs[iXc][1][0], fTy1 = fTs[iXc][1][1];
		fTx2 = fTs[iXc][2][0], fTy2 = fTs[iXc][2][1];
		fAx0 = fAs[iXc][0][0], fAy0 = fAs[iXc][0][1];
		fAx1 = fAs[iXc][1][0], fAy1 = fAs[iXc][1][1];
		fAx2 = fAs[iXc][2][0], fAy2 = fAs[iXc][2][1];

		if (!(fTy0 <= fTy1 && fTy1 <= fTy2))
			continue;
		fTx1 -= fTx0, fTy1 -= fTy0, fTx2 -= fTx0, fTy2 -= fTy0, fTx0 = fTy0 = 0;
		fAx1 -= fAx0, fAy1 -= fAy0, fAx2 -= fAx0, fAy2 -= fAy0, fAx0 = fAy0 = 0;
		double fTxPoint = fTy1 / fTy2;
		double fTxDeltaTx = fTx2 * fTxPoint - fTx1;
		double fTxDeltaTy = 0;
		double fTxDeltaAx = (fTxDeltaTx < 0) ? (fAx1 - fAx2 * fTxPoint) : (fAx2 * fTxPoint - fAx1);
		double fTxDeltaAy = (fTxDeltaTx < 0) ? (fAy1 - fAy2 * fTxPoint) : (fAy2 * fTxPoint - fAy1);
		double fTxLengthA = calcTexc_Length(fTxDeltaAx, fTxDeltaAy);
		double fTxLengthT = calcTexc_Length(fTxDeltaTx * cxTex, 0);
		fTexsX = fTxLengthA / fTxLengthT;
		BOOL bMirrorX = (iXc & 1) ? TRUE : FALSE;
		double fTxAngleA = calcTexc_AngleXYd(fTxDeltaAx, fTxDeltaAy);

		for (int iYc = 0; iYc < 6; iYc++) {
			fTx0 = fTs[iYc][0][0], fTy0 = fTs[iYc][0][1];
			fTx1 = fTs[iYc][1][0], fTy1 = fTs[iYc][1][1];
			fTx2 = fTs[iYc][2][0], fTy2 = fTs[iYc][2][1];
			fAx0 = fAs[iYc][0][0], fAy0 = fAs[iYc][0][1];
			fAx1 = fAs[iYc][1][0], fAy1 = fAs[iYc][1][1];
			fAx2 = fAs[iYc][2][0], fAy2 = fAs[iYc][2][1];

			if (!(fTx0 <= fTx1 && fTx1 <= fTx2))
				continue;
			fTx1 -= fTx0, fTy1 -= fTy0, fTx2 -= fTx0, fTy2 -= fTy0, fTx0 = fTy0 = 0;
			fAx1 -= fAx0, fAy1 -= fAy0, fAx2 -= fAx0, fAy2 -= fAy0, fAx0 = fAy0 = 0;
			double fTyPoint = fTx1 / fTx2;
			double fTyDeltaTx = 0;
			double fTyDeltaTy = fTy2 * fTyPoint - fTy1;
			double fTyDeltaAx = (fTyDeltaTy < 0) ? (fAx2 * fTyPoint - fAx1) : (fAx1 - fAx2 * fTyPoint);
			double fTyDeltaAy = (fTyDeltaTy < 0) ? (fAy2 * fTyPoint - fAy1) : (fAy1 - fAy2 * fTyPoint);
			double fTyLengthA = calcTexc_Length(fTyDeltaAx, fTyDeltaAy);
			double fTyLengthT = calcTexc_Length(0, fTyDeltaTy * cyTex);
			fTexsY = fTyLengthA / fTyLengthT;
			BOOL bMirrorY = (iYc & 1) ? TRUE : FALSE;
			double fTyAngleA = calcTexc_AngleXYd(fTyDeltaAx, fTyDeltaAy);

			double fTxTyAngleA = calcDEGFrom0To360(fTxAngleA - fTyAngleA);

			fTexAngle = fTxAngleA;

			if (89 < fTxTyAngleA && fTxTyAngleA < 91) {
				// -sx +sy
				fTexsX = -fTexsX;
				fTexAngle += 180;
			} else if (269 < fTxTyAngleA && fTxTyAngleA < 271) {
				// +sx +sy
			} else {
				// ???
				return FALSE;
			}
			fTexAngle = calcDEGFrom0To360(fTexAngle);

			double fEx, fEy;
			calcTexc_rotateTriPos(
				fAs[0][0][0], fAs[0][0][1],
				fAs[0][1][0], fAs[0][1][1],
				fAs[0][2][0], fAs[0][2][1],
				cxTex, cyTex, fTexsX, fTexsY, fTexAngle, fEx, fEy
				);
			fEx -= fTs[0][0][0];
			fEy -= fTs[0][0][1];
			if (fabs(fEx) < 1E-5) fEx = 0;
			if (fabs(fEy) < 1E-5) fEy = 0;
			fTexcX = (int)floor((-1.0 * fEx * cxTex) + 0.5);
			fTexcY = (int)floor((-1.0 * fEy * cyTex) + 0.5);
			while ((int)floor(fTexcX + 0.5) < 0) fTexcX += cxTex;
			while ((int)floor(fTexcY + 0.5) < 0) fTexcY += cyTex;
			while ((int)floor(fTexcX + 0.5) >= (int)cxTex) fTexcX -= cxTex;
			while ((int)floor(fTexcY + 0.5) >= (int)cyTex) fTexcY -= cyTex;

			ASSERT(TRUE);
			return TRUE;
		}
	}
	return FALSE;
}

//	BOOL calcTexc_TouchToX(
//		double fTx0, double fTy0,
//		double fTx1, double fTy1,
//		double fTx2, double fTy2,
//		double fAx0, double fAy0,
//		double fAx1, double fAy1,
//		double fAx2, double fAy2,
//		int nMirror,
//		UINT cxTex, UINT cyTex,
//		double &fTexsX
//		)
//	{
//		if (!(fTy0 <= fTy1 && fTy1 <= fTy2))
//			return FALSE;
//		fTx1 -= fTx0, fTy1 -= fTy0, fTx2 -= fTx0, fTy2 -= fTy0, fTx0 = fTy0 = 0;
//		fAx1 -= fAx0, fAy1 -= fAy0, fAx2 -= fAx0, fAy2 -= fAy0, fAx0 = fAy0 = 0;
//		double fTouch = fTy1 / fTy2;
//		double fDeltaTx = fTx2 * fTouch - fTx1;
//		double fDeltaTy = 0;
//		double fDeltaAx0 = fAx2 * fTouch - fAx1;
//		double fDeltaAy0 = fAy2 * fTouch - fAy1;
//		double fDeltaAx1 = (fDeltaTx < 0) ? (-fDeltaAx0) : (+fDeltaAx0);
//		double fDeltaAy1 = (fDeltaTx < 0) ? (-fDeltaAy0) : (+fDeltaAy0);
//		double fDelta1Angle = calcTexc_AngleXY(fDeltaAx1, fDeltaAy1) / 3.1415926535897932384626433832795 * 180.0;
//		double fLineA = calcTexc_Length(fDeltaAx0, fDeltaAy0);
//		double fLineT = calcTexc_Length(fDeltaTx * cxTex, 0);
//		fTexsX = fLineA / fLineT;
//		while (fDelta1Angle < 0) fDelta1Angle += 360;
//		return TRUE;
//	}

//	BOOL calcTexc_TouchToY(
//		double fTx0, double fTy0,
//		double fTx1, double fTy1,
//		double fTx2, double fTy2,
//		double fAx0, double fAy0,
//		double fAx1, double fAy1,
//		double fAx2, double fAy2,
//		int nMirror,
//		double &fTexsY
//		)
//	{
//		return TRUE;
//	}

//	BOOL calcTexc_TouchAngle(
//		double fTx0, double fTy0,
//		double fTx1, double fTy1,
//		double fTx2, double fTy2,
//		double fAx0, double fAy0,
//		double fAx1, double fAy1,
//		double fAx2, double fAy2,
//		BOOL bMirror,
//		double &fRA
//		)
//	{
//		// テクスチャ座標と頂点座標を用いて、描写されるテクスチャがベ－スラインに
//		// 対して何度傾いているかを計算する。
//		if (!(fTy0 <= fTy1 && fTy1 <= fTy2))
//			return FALSE;
//		fTx1 -= fTx0, fTy1 -= fTy0,
//		fTx2 -= fTx0, fTy2 -= fTy0,
//		fTx0 = fTy0 = 0;
//		fAx1 -= fAx0, fAy1 -= fAy0,
//		fAx2 -= fAx0, fAy2 -= fAy0,
//		fAx0 = fAy0 = 0;
//		if (!(0 <= fTy1 && 0 <= fAy1) && !(0 >= fTy1 && 0 >= fAy1))
//			return FALSE;
//		BOOL bFlipH = !(0 <= fTx2 && 0 <= fAx2) && !(0 >= fTx2 && 0 >= fAx2);
//		BOOL bFlipV = !(0 <= fTy2 && 0 > fAy2) && !(0 >= fTy2 && 0 < fAy2);
//	
//		ASSERT(fabs(fTy2) > 1E-5);
//		double fTf = (fTy1 / fTy2);
//		double fDx0 = fTf * fAx2;
//		double fDy0 = fTf * fAy2;
//		int v1 = (fTx0 < fTx1); // 0=右タッチ,1=左タッチ
//		int v2 = v1 ? (0 <= fAx1 - fDx0) : (0 <= fDx0 - fAx1);
//		int v3 = v1 ? (0 <= fAy1 - fDy0) : (0 <= fDy0 - fAy1);
//		int v4 = fTx2 < fTx1;
//		int v5 = fTy2 < fTy1;
//		double fDx0t = fAx1 - fDx0;
//		double fDy0t = fAy1 - fDy0;
//		//if (bFlipH) fDx0t = -fDx0t;
//		//if (bFlipV) fDy0t = -fDy0t;
//		if (v1 != 0) {
//			fRA = calcTexc_AngleXY(+fDx0t, +fDy0t);
//		} else if (v1 == 0) {
//			fRA = calcTexc_AngleXY(-fDx0t, -fDy0t);
//		} else {
//			ASSERT(FALSE);
//		}
//		ASSERT(!(bMirror && v1));
//		return TRUE;
//	}

BOOL calcTexcAndReverse(
	UINT cxTex, UINT cyTex, 
	double fAx0, double fAy0, double fTx0, double fTy0, 
	double fAx1, double fAy1, double fTx1, double fTy1,
	double fAx2, double fAy2, double fTx2, double fTy2,
	UINT iMirror,
	double &fTexcX,
	double &fTexcY,
	double &fTexAngle,
	double &fTexsX,
	double &fTexsY
	)
{
#if 0
//	double z1[][2] = {
//		fAx0, fAy0,
//		fAx1, fAy1,
//		fAx2, fAy2,
//	};
//	double z2[][2] = {
//		fTx0, fTy0,
//		fTx1, fTy1,
//		fTx2, fTy2,
//	};
//	CString s, t;
//	for (int a = 0; a < 2; a++) {
//		double (*zz)[2] = (a == 0) ? z1 : z2;
//		for (int b = 0; b < 2; b++) {
//			for (int c = 0; c < 3; c++) {
//				double x0 = zz[((c - 1) + 3) % 3][0];
//				double y0 = zz[((c - 1) + 3) % 3][1];
//				double x1 = zz[((c    ) + 3) % 3][0];
//				double y1 = zz[((c    ) + 3) % 3][1];
//				double x2 = zz[((c + 1) + 3) % 3][0];
//				double y2 = zz[((c + 1) + 3) % 3][1];
//				double l01 = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
//				double l12 = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
//
//				double r21 = acos((x2 - x1) / l12);
//				double g21 = r21 / 3.1415926535897932384626433832795 * 180.0;
//				if (y2 - y1 < 0)
//					g21 = 360 - g21;
//				double r01 = acos((x0 - x1) / l01);
//				double g01 = r01 / 3.1415926535897932384626433832795 * 180.0;
//				if (y0 - y1 < 0)
//					g01 = 360 - g01;
//				if (b == 0) {
//					t.Format("%lf,%lf,%d,%lf", x1, y1, (int)g21, (double)l12);
//				} else {
//					double rr = (g21 - g01);
//					while (rr < 0) rr += 360;
//					while (rr > 360) rr -= 360;
//					if (rr > 180)
//						rr = 360 - rr;
//					t.Format(",,%lf", (rr));
//				}
//				s += t;
//				s += "\n";
//				t.Empty();
//			}
//		}
//	}
//	TRACE1("---\n%s", s);
#endif
	{
		if (!calcTexc_TouchToXY(fTx0, fTy0, fTx1, fTy1, fTx2, fTy2, fAx0, fAy0, fAx1, fAy1, fAx2, fAy2, cxTex, cyTex, iMirror, fTexsX, fTexsY, fTexAngle, fTexcX, fTexcY)) {
			//ASSERT(FALSE);
			return FALSE;
		}
	}
	return TRUE;
#if 0
//	{
//		double fRA = 0;
//		if (!calcTexc_TouchAngle(fTx0, fTy0, fTx1, fTy1, fTx2, fTy2, fAx0, fAy0, fAx1, fAy1, fAx2, fAy2, 0, fRA)) {
//			if (!calcTexc_TouchAngle(fTx1, fTy1, fTx2, fTy2, fTx0, fTy0, fAx1, fAy1, fAx2, fAy2, fAx0, fAy0, 0, fRA)) {
//				if (!calcTexc_TouchAngle(fTx2, fTy2, fTx0, fTy0, fTx1, fTy1, fAx2, fAy2, fAx0, fAy0, fAx1, fAy1, 0, fRA)) {
//					if (!calcTexc_TouchAngle(fTx0, fTy0, fTx2, fTy2, fTx1, fTy1, fAx0, fAy0, fAx2, fAy2, fAx1, fAy1, 1, fRA)) {
//						if (!calcTexc_TouchAngle(fTx1, fTy1, fTx0, fTy0, fTx2, fTy2, fAx1, fAy1, fAx0, fAy0, fAx2, fAy2, 1, fRA)) {
//							if (!calcTexc_TouchAngle(fTx2, fTy2, fTx1, fTy1, fTx0, fTy0, fAx2, fAy2, fAx1, fAy1, fAx0, fAy0, 1, fRA)) {
//								ASSERT(FALSE);
//								return FALSE;
//							}
//						}
//					}
//				}
//			}
//		}
//		fTexAngle = fRA / 3.1415926535897932384626433832795 * 180.0;
//		if (fTexAngle < 0)
//			fTexAngle += 360;
//		if (!calcTexc_TouchToX(fTx0, fTy0, fTx1, fTy1, fTx2, fTy2, fAx0, fAy0, fAx1, fAy1, fAx2, fAy2, 0, cxTex, cyTex, fTexsX)) {
//			if (!calcTexc_TouchToX(fTx1, fTy1, fTx2, fTy2, fTx0, fTy0, fAx1, fAy1, fAx2, fAy2, fAx0, fAy0, 0, cxTex, cyTex, fTexsX)) {
//				if (!calcTexc_TouchToX(fTx2, fTy2, fTx0, fTy0, fTx1, fTy1, fAx2, fAy2, fAx0, fAy0, fAx1, fAy1, 0, cxTex, cyTex, fTexsX)) {
//					if (!calcTexc_TouchToX(fTx0, fTy0, fTx2, fTy2, fTx1, fTy1, fAx0, fAy0, fAx2, fAy2, fAx1, fAy1, 1, cxTex, cyTex, fTexsX)) {
//						if (!calcTexc_TouchToX(fTx1, fTy1, fTx0, fTy0, fTx2, fTy2, fAx1, fAy1, fAx0, fAy0, fAx2, fAy2, 1, cxTex, cyTex, fTexsX)) {
//							if (!calcTexc_TouchToX(fTx2, fTy2, fTx1, fTy1, fTx0, fTy0, fAx2, fAy2, fAx1, fAy1, fAx0, fAy0, 1, cxTex, cyTex, fTexsX)) {
//								ASSERT(FALSE);
//								return FALSE;
//							}
//						}
//					}
//				}
//			}
//		}
//		if (!calcTexc_TouchToY(fTx0, fTy0, fTx1, fTy1, fTx2, fTy2, fAx0, fAy0, fAx1, fAy1, fAx2, fAy2, 0, fTexsX)) {
//			if (!calcTexc_TouchToY(fTx1, fTy1, fTx2, fTy2, fTx0, fTy0, fAx1, fAy1, fAx2, fAy2, fAx0, fAy0, 0, fTexsX)) {
//				if (!calcTexc_TouchToY(fTx2, fTy2, fTx0, fTy0, fTx1, fTy1, fAx2, fAy2, fAx0, fAy0, fAx1, fAy1, 0, fTexsX)) {
//					if (!calcTexc_TouchToY(fTx0, fTy0, fTx2, fTy2, fTx1, fTy1, fAx0, fAy0, fAx2, fAy2, fAx1, fAy1, 1, fTexsX)) {
//						if (!calcTexc_TouchToY(fTx1, fTy1, fTx0, fTy0, fTx2, fTy2, fAx1, fAy1, fAx0, fAy0, fAx2, fAy2, 1, fTexsX)) {
//							if (!calcTexc_TouchToY(fTx2, fTy2, fTx1, fTy1, fTx0, fTy0, fAx2, fAy2, fAx1, fAy1, fAx0, fAy0, 1, fTexsX)) {
//								ASSERT(FALSE);
//								return FALSE;
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	fTexcX = fTexcY = fTexAngle = 0;
//	fTexsX = fTexsY = 1;
//	double fEx0 = fTx0, fEy0 = fTy0;
//	double fEx1 = fTx1, fEy1 = fTy1;
//	double fEx2 = fTx2, fEy2 = fTy2;
//	calcRevAtMed(fEy0, fEy1, fEy2);
//	{
//		double fXd1 = fAx1 - fAx0;
//		double fYd1 = fAy1 - fAy0;
//		double fLd1 = sqrt(fXd1 * fXd1 + fYd1 * fYd1);
//		ASSERT(fLd1 != 0);
//		double fXn1 = fXd1 / fLd1;
//		double fYn1 = fYd1 / fLd1;
//		double fXr1 = acos(fXn1) / 3.1415926535897932384626433832795 * 180.0;
//		if (fYn1 < 0)
//			fXr1 = 360 - fXr1;
//		fXr1 = clampVal(fXr1, 0);
//	
//		double fXd2 = fAx2 - fAx0;
//		double fYd2 = fAy2 - fAy0;
//		double fLd2 = sqrt(fXd2 * fXd2 + fYd2 * fYd2);
//		ASSERT(fLd2 != 0);
//		double fXn2 = fXd2 / fLd2;
//		double fYn2 = fYd2 / fLd2;
//		double fXr2 = acos(fXn2) / 3.1415926535897932384626433832795 * 180.0;
//		if (fYn2 < 0)
//			fXr2 = 360 - fXr2;
//		fXr2 = clampVal(fXr2, 0);
//	
//		double fXtd1 = clampVal(fEx1 - fEx0, 3) * cxTex;
//		double fYtd1 = clampVal(fEy1 - fEy0, 3) * cyTex;
//		double fLtd1 = sqrt(fXtd1 * fXtd1 + fYtd1 * fYtd1);
//		ASSERT(fLtd1 != 0);
//		double fXtld1 = (fEx1 - fEx0) * fTexsX;
//		double fYtld1 = (fEy1 - fEy0) * fTexsY;
//		double fLtld1 = sqrt(fXtld1 * fXtld1 + fYtld1 * fYtld1);
//		double fXtn1 = fXtld1 / fLtld1;
//		double fYtn1 = fYtld1 / fLtld1;
//		double fXtr1 = acos(fXtn1) / 3.1415926535897932384626433832795 * 180.0;
//		if (fYtn1 < 0)
//			fXtr1 = 360 - fXtr1;
//		fXtr1 = clampVal(fXtr1, 0);
//		double fXtrr1 = calc360RevAngle(fXtr1);
//	
//		double fXtd2 = clampVal(fEx2 - fEx0, 3) * cxTex;
//		double fYtd2 = clampVal(fEy2 - fEy0, 3) * cyTex;
//		double fLtd2 = sqrt(fXtd2 * fXtd2 + fYtd2 * fYtd2);
//		ASSERT(fLtd2 != 0);
//		double fXtld2 = (fEx2 - fEx0) * fTexsX;
//		double fYtld2 = (fEy2 - fEy0) * fTexsY;
//		double fLtld2 = sqrt(fXtld2 * fXtld2 + fYtld2 * fYtld2);
//		double fXtn2 = fXtld2 / fLtld2;
//		double fYtn2 = fYtld2 / fLtld2;
//		double fXtr2 = acos(fXtn2) / 3.1415926535897932384626433832795 * 180.0;
//		if (fYtn2 < 0)
//			fXtr2 = 360 - fXtr2;
//		fXtr2 = clampVal(fXtr2, 0);
//		double fXtrr2 = calc360RevAngle(fXtr2);
//	
//		fTexAngle = fXr1 - fXtrr1;
//	
//		double fXb0 = 0;
//		double fYb0 = 0;
//	
//		rotateXYc(fEx0, fEy0, fXb0, fYb0, fEx0, fEy0, -fTexAngle);
//		rotateXYc(fEx1, fEy1, fXb0, fYb0, fEx1, fEy1, -fTexAngle);
//		rotateXYc(fEx2, fEy2, fXb0, fYb0, fEx2, fEy2, -fTexAngle);
//	
//		double fXa0 = 0;
//		double fYa0 = 0;
//		rotateXYc(fAx0, fAy0, 0, 0, fXa0, fYa0, -fTexAngle);
//	
//		ASSERT(TRUE);
//		{
//			double fXtd1 = clampVal(fEx1 - fEx0, 3) * cxTex;
//			double fYtd1 = clampVal(fEy1 - fEy0, 3) * cyTex;
//			double fLtd1 = sqrt(fXtd1 * fXtd1 + fYtd1 * fYtd1);
//			ASSERT(fLtd1 != 0);
//			double fXtld1 = (fEx1 - fEx0) * fTexsX;
//			double fYtld1 = (fEy1 - fEy0) * fTexsY;
//			double fLtld1 = sqrt(fXtld1 * fXtld1 + fYtld1 * fYtld1);
//			double fXtn1 = fXtld1 / fLtld1;
//			double fYtn1 = fYtld1 / fLtld1;
//			double fXtr1 = acos(fXtn1) / 3.1415926535897932384626433832795 * 180.0;
//			if (fYtn1 < 0)
//				fXtr1 = 360 - fXtr1;
//			fXtr1 = clampVal(fXtr1, 0);
//			double fXtrr1 = calc360RevAngle(fXtr1);
//	
//			double fXtd2 = clampVal(fEx2 - fEx0, 3) * cxTex;
//			double fYtd2 = clampVal(fEy2 - fEy0, 3) * cyTex;
//			double fLtd2 = sqrt(fXtd2 * fXtd2 + fYtd2 * fYtd2);
//			ASSERT(fLtd2 != 0);
//			double fXtld2 = (fEx2 - fEx0) * fTexsX;
//			double fYtld2 = (fEy2 - fEy0) * fTexsY;
//			double fLtld2 = sqrt(fXtld2 * fXtld2 + fYtld2 * fYtld2);
//			double fXtn2 = fXtld2 / fLtld2;
//			double fYtn2 = fYtld2 / fLtld2;
//			double fXtr2 = acos(fXtn2) / 3.1415926535897932384626433832795 * 180.0;
//			if (fYtn2 < 0)
//				fXtr2 = 360 - fXtr2;
//			fXtr2 = clampVal(fXtr2, 0);
//			double fXtrr2 = calc360RevAngle(fXtr2);
//	
//			if (fXn1 != 0) {
//				ASSERT(fXtd1 != 0);
//				fTexsX = fLd1 * fXn1 / fXtd1;
//			} else if (fXn2 != 0) {
//				ASSERT(fXtd2 != 0);
//				fTexsX = fLd2 * fXn2 / fXtd2;
//			} else {
//				ASSERT(FALSE);
//			}
//			if (fYn1 != 0) {
//				ASSERT(fYtd1 != 0);
//				fTexsY = fLd1 * fYn1 / fYtd1;
//			} else if (fYn2 != 0) {
//				ASSERT(fYtd2 != 0);
//				fTexsY = fLd2 * fYn2 / fYtd2;
//			} else {
//				ASSERT(FALSE);
//			}
//			fTexsX = clampVal(fTexsX, 2);
//			fTexsY = clampVal(fTexsY, 2);
//			ASSERT(fTexsX != 0 && fTexsY != 0);
//	
//			fTexcX = fTx0 * cxTex - fAx0 / fTexsX;
//			fTexcY = fTy0 * cyTex - fAy0 / fTexsY;
//			fTexcX = (int)clampVal(fTexcX, 2) % cxTex;
//			fTexcY = (int)clampVal(fTexcY, 2) % cyTex;
//			fTexsY =-fTexsY;
//		}
//	}
#endif
}

void calcTextureName(
	const CByteArray &arrTexture,
	UINT iIdx,
	CString &tstrTexNameLong,
	CString &tstrTexNameShort,
	CString &tstrTexAdditional
	)
{
	tstrTexNameLong = _T("common/black");
	tstrTexNameShort = _T("common/black");
	tstrTexAdditional.Empty();
	if ((UINT)arrTexture.GetSize() <= iIdx * 140)
		return;
	const tex3_t *lpTex3 = (const tex3_t *)(&arrTexture.GetData()[140 * iIdx]);
	tstrTexNameShort = (tstrTexNameLong = lpTex3->name);
	if (tstrTexNameShort.Left(9).CompareNoCase(_T("textures/")) == 0)
		tstrTexNameShort.Delete(0, 9);
	if (lpTex3->contents & 0x8000000) {
		tstrTexAdditional += _T("+surfaceparm");
		tstrTexAdditional += _T(" detail");
	}
}

BOOL ParseChunk2EntityLines(const CByteArray &arr, CStringList &l)
{
	CString tstrLine;
	CByteArray arrTemp;
	const BYTE *lpbPos = arr.GetData();
	const BYTE *lpbEnd = lpbPos + arr.GetSize();
	const BYTE *lpbRem = lpbPos;
	while (lpbPos != lpbEnd) {
		if (*lpbPos == '\n') {
			UINT n = lpbPos - lpbRem;
			arrTemp.SetSize(n + 1);
			strncpy((char *)arrTemp.GetData(), (const char *)lpbRem, n);
			arrTemp.SetAt(n, 0);
			tstrLine = arrTemp.GetData();
			tstrLine.TrimLeft();
			l.AddTail(tstrLine);

			lpbRem = lpbPos = (lpbPos + 1);
		} else {
			lpbPos++;
		}
	}
	return TRUE;
}

BOOL ParseKeyNVal(
	CString tstr,
	CString &tstrKey,
	CString &tstrVal
)
{
	int p0 = -1;
	int p1 = tstr.Find(TCHAR('\"'), p0 + 1);
	int p2 = tstr.Find(TCHAR('\"'), p1 + 1);
	int p3 = tstr.Find(TCHAR('\"'), p2 + 1);
	int p4 = tstr.Find(TCHAR('\"'), p3 + 1);
	if (p1 < 0 || p2 < 0 || p3 < 0 || p4 < 0)
		return FALSE;
	tstrKey = tstr.Mid(p1 + 1, p2 - p1 - 1);
	tstrVal = tstr.Mid(p3 + 1, p4 - p3 - 1);
	return TRUE;
}

void Write2(
	CFile &fileInto,
	LPCTSTR lpsz1Line
)
{
	LPCSTR lpsz1 = T2CA(lpsz1Line);
	fileInto.Write((void *)lpsz1, strlen(lpsz1));
	fileInto.Write((void *)"\r\n", 2);
}

BOOL WriteModel2(
	CFile &fileInto,
	int iModelNo,
	myVtx_t vtxMdlPos,
	BOOL bAddBBox,
	CFindTexMetr &rTexMetr,
	const CByteArray &arrTexture,
	const CByteArray &arrFace,
	const CByteArray &arrVertex,
	const CByteArray &arrModel,
	const CByteArray &arrMeshVerts,
	const CByteArray &arrTerr
)
{
	int iIndex = 1;
	CString tstrLine, tstrTemp, tstrVertHint, tstrPerTri;

	int iMax_PolysArray = 0;
	int iMax_VertexArray = 0;

	if (iModelNo == 0) {
		// テレイン・パッチ展開
		const nTerrs = arrTerr.GetSize() / 388;
		CArray<myTerrain_t, myTerrain_t> arrTerrs;
		arrTerrs.SetSize(nTerrs, 0);
		{
			CUIntArray arrTerBox;
			arrTerBox.SetSize(32 * 32, 0);
			for (int i = 0; i < 64; i++) {
				int nRestTotal = 0;
				CByteArray arrTerUsed;
				arrTerUsed.SetSize(nTerrs);
				{
					const ix = i % 8, iy = i / 8;
					for (int iTerr = 0; iTerr < nTerrs; iTerr++) {
						arrTerUsed[iTerr] = 0;
						const BYTE *lpbEach = arrTerr.GetData() + 388 * iTerr;
						int x = (CHAR)lpbEach[0x24];
						int y = (CHAR)lpbEach[0x25];
						int bx = x % 8;
						int by = y % 8;
						if (bx < 0) bx += 8;
						if (by < 0) by += 8;
						ASSERT(0 <= bx && bx < 8);
						ASSERT(0 <= by && by < 8);
						if (bx != ix || by != iy)
							continue;
						nRestTotal++;
					}
				}
				while (nRestTotal != 0) {
					for (int t = 0; t < 32 * 32; t++)
						arrTerBox[t] = (UINT)-1;
					const ix = i % 8, iy = i / 8;
					int iTerr;
					int nRest = 0;
					int nHeight = 0;
					for (iTerr = 0; iTerr < nTerrs; iTerr++) {
						if (arrTerUsed[iTerr] != 0)
							continue;
						const BYTE *lpbEach = arrTerr.GetData() + 388 * iTerr;
						int x = (CHAR)lpbEach[0x24];
						int y = (CHAR)lpbEach[0x25];
						int z = *(SHORT *)&lpbEach[0x26];
						int bx = x % 8;
						int by = y % 8;
						if (bx < 0) bx += 8;
						if (by < 0) by += 8;
						ASSERT(0 <= bx && bx < 8);
						ASSERT(0 <= by && by < 8);
						if (bx != ix || by != iy)
							continue;
						if (nRest != 0 && nHeight != z) {
							continue;
						}
						int e = ((BYTE)x) / 8 + ((BYTE)y) / 8 * 32;
						arrTerBox.ElementAt(e) = iTerr;
						if (nRest == 0)
							nHeight = z;
						nRest++;
						nRestTotal--;
						arrTerUsed[iTerr] = 1;
					}
					while (nRest != 0) {
						//{
						//	FILE *lpf = fopen("H:\\DEV.TXT", "wt");
						//	if (lpf) {
						//		for (int y = 0; y < 256; y++) {
						//			for (int x = 0; x < 256; x++) {
						//				const px = x;
						//				const py = y;
						//				const pe = px + py * 256;
						//				DWORD z = arrTerBox[pe];
						//				fprintf(lpf, "%03d,%03d,%5d %08lX\n", x, y, pe, z);
						//			}
						//		}
						//		fclose(lpf);
						//	}
						//}
						int maxcx = 0, maxcy = 0, maxpx = 0, maxpy = 0;
						for (int sy = 0; sy < 32; sy++) {
							for (int sx = 0; sx < 32; sx++) {
								const px = sx;
								const py = sy;
								const pe = px + py * 32;
								ASSERT(0 <= px && px < 32);
								ASSERT(0 <= py && py < 32);
								if (arrTerBox[pe] == (UINT)-1)
									continue;
								for (int fy = 1; fy < 32; fy++) {
									BOOL bT2D = FALSE;
									for (int fx = 1; fx < 32; fx++) {
										BOOL bBlocked = FALSE;
										for (int vy = 0; !bBlocked && vy < fy; vy++) {
											for (int vx = 0; !bBlocked && vx < fx; vx++) {
												const px = (sx + vx) % 32;
												const py = (sy + vy) % 32;
												ASSERT(0 <= px && px < 32);
												ASSERT(0 <= py && py < 32);
												const pe = px + py * 32;
												if (arrTerBox[pe] != (UINT)-1) {
													continue;
												}
												bBlocked = TRUE;
											}
										}
										if (!bBlocked) {
											if (maxcx * maxcy < (fx) * (fy)) {
												maxcx = fx;
												maxcy = fy;
												maxpx = sx;
												maxpy = sy;
											}
										} else {
											if (fx == 1)
												bT2D = TRUE;
											break;
										}
									}
									if (bT2D) {
										break;
									}
								}
							}
						}
						ASSERT(maxcx != 0 && maxcy != 0);

						tstrLine.Format(
							_T("// Brush #%d"), iIndex
							);
						Write2(fileInto, tstrLine);
						Write2(fileInto, _T("{"));
						Write2(fileInto, _T("terrainDef"));
						Write2(fileInto, _T("{"));

						int iDirZ;
						double ox, oy, oz;
						{
							const px = maxpx;
							const py = maxpy;
							const pe = px + py * 32;
							ASSERT(0 <= px && px < 32);
							ASSERT(0 <= py && py < 32);
							UINT iIdx = arrTerBox[pe];
							ASSERT(iIdx < (UINT)nTerrs);
							const BYTE *lpbEach = arrTerr.GetData() + 388 * iIdx;
							ox = 64 * (int)((CHAR)lpbEach[0x24]);
							oy = 64 * (int)((CHAR)lpbEach[0x25]);
							oz = *(SHORT *)&lpbEach[0x26];
							iDirZ = (lpbEach[0x00] & 0x40) ? 1 : 0;
						}

						tstrTemp.Format(
							_T("%d %d %d"), 9 + 8 * (maxcx - 1), 9 + 8 * (maxcy - 1), iDirZ
							);
						Write2(fileInto, tstrTemp);
						tstrTemp.Format(
							_T("%lf %lf %lf"), (double)ox, (double)oy, (double)oz
							);
						Write2(fileInto, tstrTemp);

						Write2(fileInto, _T("{"));
						for (int ty = 0; ty < maxcy + 1; ty++) {
							for (int tx = 0; tx < maxcx + 1; tx++) {
								const px = (tx + maxpx) % 32;
								const py = (ty + maxpy) % 32;
								const pe = px + py * 32;
								CString tstrTexName;
								CString tstrTexNamec;
								if (tx != maxcx && ty != maxcy) {
									UINT iIdx = arrTerBox[pe];
									ASSERT(iIdx < (UINT)nTerrs);
									ASSERT(iIdx != (UINT)-1);
									const BYTE *lpbEach = arrTerr.GetData() + 388 * iIdx;
									UINT iTex = *(WORD *)&lpbEach[0x28];
									tex3_t *lpTex3 = (tex3_t *)(((const BYTE *)arrTexture.GetData()) + 140 * iTex);
									tstrTexName = lpTex3->name;
									tstrTexNamec = tstrTexName;
									if (tstrTexNamec.Left(9) == _T("textures/"))
										tstrTexNamec.Delete(0, 9);
								} else {
									tstrTexName = _T("notexture");
									tstrTexNamec = _T("notexture");
								}

								tstrTemp.Format(
									_T("%ld %ld ( %s %.0lf %.0lf %.2lf %.0lf %lf %lf %ld %ld %ld )"),
									(DWORD)0, (DWORD)0,
									tstrTexNamec,
									(double)0, (double)0, (double)0, (double)0, (double)1, (double)1,
									(DWORD)0, (DWORD)0, (DWORD)0
								);
								Write2(fileInto, tstrTemp);
							}
						}
						Write2(fileInto, _T("}"));

						Write2(fileInto, _T("{"));
						const maxtx = 9 + 8 * (maxcx - 1);
						const maxty = 9 + 8 * (maxcy - 1);
						for (int gy = 0; gy < maxty; gy++) {
							for (int gx = 0; gx < maxtx; gx++) {
								const px = ((gx - 1) / 8 + maxpx) % 32;
								const py = ((gy - 1) / 8 + maxpy) % 32;
								const pe = px + py * 32;
								UINT iIdx = arrTerBox[pe];
								ASSERT(iIdx < (UINT)nTerrs);
								const BYTE *lpbEach = arrTerr.GetData() + 388 * iIdx;

								const mx = gx - ((gx - 1) / 8) * 8;
								const my = gy - ((gy - 1) / 8) * 8;
								int h = ((int)lpbEach[0x130 + mx + my * 9]) * 2;

								tstrTemp.Format(
									_T("%lf ( %s ) ( %s )"),
									(double)h,
									_T(""),
									_T("")
									);
								Write2(fileInto, tstrTemp);
							}
						}
						for (int ly = 0; ly < maxcy; ly++) {
							for (int lx = 0; lx < maxcx; lx++) {
								const px = (lx + maxpx) % 32;
								const py = (ly + maxpy) % 32;
								const pe = px + py * 32;

								arrTerBox[pe] = (UINT)-1;
								nRest--;
							}
						}
						Write2(fileInto, _T("}"));

						Write2(fileInto, _T("}"));
						Write2(fileInto, _T("}"));
						iIndex++;
					}
				}
			}
		}

#if 0
		//	for (iTerr = 0; iTerr < nTerrs; iTerr++) {
		//		const BYTE *lpbEach = arrTerr.GetData() + 388 * iTerr;
		//		myTerrain_t &r = arrTerrs.ElementAt(iTerr);
		//		r.origin.x = 64 * (int)((CHAR)lpbEach[0x24]);
		//		r.origin.y = 64 * (int)((CHAR)lpbEach[0x25]);
		//		r.origin.z = *(__int16 *)&lpbEach[0x26];
		//		r.iTex = *(unsigned __int16 *)&lpbEach[0x28];
		//		int o = 0;
		//		for (int y = 0; y < 9; y++) {
		//			for (int x = 0; x < 9; x++, o++) {
		//				r.nBit[y][x] = ((int)lpbEach[0x130 + o]) * 2;
		//			}
		//		}
		//	}
		//	for (iTerr = 0; iTerr < nTerrs; iTerr++) {
		//		myTerrain_t &r = arrTerrs.ElementAt(iTerr);
		//	
		//		tex3_t *lpTex3 = (tex3_t *)(((const BYTE *)arrTexture.GetData()) + 140 * r.iTex);
		//		CString tstrTexName = lpTex3->name;
		//		CString tstrTexNamec = tstrTexName;
		//		if (tstrTexNamec.Left(9) == _T("textures/"))
		//			tstrTexNamec.Delete(0, 9);
		//	
		//		tstrLine.Format(
		//			_T("// Brush #%d"), iIndex
		//			);
		//		Write2(fileInto, tstrLine);
		//		Write2(fileInto, _T("{"));
		//		Write2(fileInto, _T(" terrainDef"));
		//		Write2(fileInto, _T(" {"));
		//		//
		//		tstrTemp.Format(
		//			_T("  %d %d %d"), 9, 9, 0
		//			);
		//		Write2(fileInto, tstrTemp);
		//		tstrTemp.Format(
		//			_T("  %lf %lf %lf"), (double)r.origin.x, (double)r.origin.y, (double)r.origin.z
		//			);
		//		Write2(fileInto, tstrTemp);
		//		Write2(fileInto, _T("  {"));
		//		tstrTemp.Format(
		//			_T("   %ld %ld ( %s %.0lf %.0lf %.2lf %.0lf %lf %lf %ld %ld %ld )"),
		//			tstrTexNamec,
		//			(DWORD)0, (DWORD)0,
		//			(double)0, (double)0, (double)0, (double)0, (double)1, (double)1,
		//			(DWORD)0, (DWORD)0, (DWORD)0
		//		);
		//		Write2(fileInto, tstrTemp);
		//		Write2(fileInto, tstrTemp);
		//		Write2(fileInto, tstrTemp);
		//		Write2(fileInto, tstrTemp);
		//		Write2(fileInto, _T("  }"));
		//	
		//		Write2(fileInto, _T("  {"));
		//		for (int y = 0; y < 9; y++) {
		//			for (int x = 0; x < 9; x++) {
		//				tstrTemp.Format(
		//					_T("   %lf ( %s ) ( %s )"),
		//					(double)r.nBit[y][x],
		//					_T(""),
		//					_T("")
		//					);
		//				Write2(fileInto, tstrTemp);
		//			}
		//		}
		//		Write2(fileInto, _T("  }"));
		//		//
		//		Write2(fileInto, _T(" }"));
		//		Write2(fileInto, _T("}"));
		//		iIndex++;
		//		continue;
		//	}
#endif
	}

	model3_t *lpMdl3 = (model3_t *)(((const BYTE *)arrModel.GetData()) + 40 * iModelNo);
	const nFaces = lpMdl3->face + lpMdl3->n_faces;
	int iFace = lpMdl3->face;
	int nBadNorm = 0, nBadTri = 0;
	const double dThin = 8.0;
	for (; iFace < nFaces; iFace++) {
		face3_t *lpFace3 = (face3_t *)(((const BYTE *)arrFace.GetData()) + 108 * (iFace));

		CString tstrTexName;
		CString tstrTexNamec;
		CString tstrTexAdd;
		calcTextureName(arrTexture, lpFace3->texture, tstrTexName, tstrTexNamec, tstrTexAdd);

		if (tstrTexAdd)
			tstrTexAdd.Insert(0, TCHAR(' '));

		if (lpFace3->type == 2) {
			const iVertBase = lpFace3->vertex;
			const nVerts = lpFace3->n_vertexes;
			const cyMtx = lpFace3->size[0];
			const cxMtx = lpFace3->size[1];

			tstrLine.Format(
				_T("// Brush #%d"), iIndex
				);
			Write2(fileInto, tstrLine);
			Write2(fileInto, _T("{"));
			Write2(fileInto, _T("patchDef2"));
			Write2(fileInto, _T("{"));
			tstrLine.Format(
				_T("%s"), tstrTexNamec
				);
			Write2(fileInto, tstrLine);
			tstrLine.Format(
				_T("( %d %d %d %d %d%s )"), (int)cyMtx, (int)cxMtx, (int)0, (int)0, (int)0, tstrTexAdd
				);
			Write2(fileInto, tstrLine);
			Write2(fileInto, _T("("));

			for (int yMtx = 0; yMtx < cyMtx; yMtx++) {
				tstrLine = _T("(");
				for (int xMtx = 0; xMtx < cxMtx; xMtx++) {
					const iVert = cyMtx * xMtx + yMtx;
					ASSERT(0 <= iVert && iVert < nVerts);
					vertex3_t *lpVtx3 = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertBase + iVert));
					tstrTemp.Format(
						_T(" ( %lf %lf %lf %lf %lf )"),
						lpVtx3->position[0], lpVtx3->position[1], lpVtx3->position[2],
						lpVtx3->texcoord[0][0], lpVtx3->texcoord[0][1]
						);
					tstrLine += tstrTemp;
				}
				tstrLine += _T(" )");
				Write2(fileInto, tstrLine);
			}

			Write2(fileInto, _T(")"));
			Write2(fileInto, _T("}"));
			Write2(fileInto, _T("}"));
			iIndex++;
			continue;
		}
		if (lpFace3->type != 1)
			continue;
		DWORD *lpMeshVerts3 = (DWORD *)(((const BYTE *)arrMeshVerts.GetData()) + 4 * (lpFace3->meshvert));
		const iVertBase = lpFace3->vertex;
		const nVerts = lpFace3->n_vertexes;
		const nTris = lpFace3->n_meshverts / 3;
		ASSERT((lpFace3->n_meshverts % 3) == 0);
		CArray<myPoly_t, myPoly_t> arrPolys;
		CArray<myVtx_t, myVtx_t> arrPolyNorms;
		arrPolys.SetSize(0, 16);
		arrPolyNorms.SetSize(0, 16);
		myVtx_t aNorm(
			lpFace3->normal[0], lpFace3->normal[1], lpFace3->normal[2]
			);
		if (aNorm.x == 0 && aNorm.y == 0 && aNorm.z == 0) {
			nBadNorm++;
			//continue;
		}
		int iNVoS = calcNVoS(aNorm);
		BOOL bTexCoordReady = FALSE;
		double fTexcX = 0, fTexcY = 0, fTexAngle = 0, fTexsX = 1.0, fTexsY = 1.0;
		for (int iTri = 0; iTri < nTris; iTri++, lpMeshVerts3 += 3) {
			const iVertIndices[3] = {
				iVertBase + lpMeshVerts3[0],
				iVertBase + lpMeshVerts3[1],
				iVertBase + lpMeshVerts3[2]
			};
			ASSERT(lpMeshVerts3[0] < (DWORD)nVerts);
			ASSERT(lpMeshVerts3[1] < (DWORD)nVerts);
			ASSERT(lpMeshVerts3[2] < (DWORD)nVerts);
			vertex3_t *lpVtx3A = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertIndices[0]));
			vertex3_t *lpVtx3B = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertIndices[1]));
			vertex3_t *lpVtx3C = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertIndices[2]));

			myTri_t aTri(
				myVtx_t(lpVtx3A->position[0], lpVtx3A->position[1], lpVtx3A->position[2]),
				myVtx_t(lpVtx3B->position[0], lpVtx3B->position[1], lpVtx3B->position[2]),
				myVtx_t(lpVtx3C->position[0], lpVtx3C->position[1], lpVtx3C->position[2])
				);
			if (!aTri.IsValidTri()) {
				nBadTri++;
				continue;
			}
			myVtx_t aTriNorm;
			calcNormVec(
				aTriNorm,
				aTri.v[0],
				aTri.v[1],
				aTri.v[2]
				);
			aTriNorm.x = clampVal(fix1P0(aTriNorm.x), 3);
			aTriNorm.y = clampVal(fix1P0(aTriNorm.y), 3);
			aTriNorm.z = clampVal(fix1P0(aTriNorm.z), 3);
			int iPoly = 0, nPolys = arrPolys.GetSize();
			for (; iPoly < nPolys; iPoly++) {
				myPoly_t &r = arrPolys.ElementAt(iPoly);
				ASSERT(!r.IsEmpty());
				myVtx_t &rNorm = arrPolyNorms.ElementAt(iPoly);
				if (aTriNorm != rNorm)
					continue;
				if (!r.JoinTri(aTri))
					continue;
				break;
			}
			if (iPoly == nPolys) {
				myPoly_t aPoly;
				aPoly.JoinTri(aTri);
				arrPolys.Add(aPoly);
				arrPolyNorms.Add(aTriNorm);
			}
			if (!bTexCoordReady) {
				UINT cxTex, cyTex;
				if (rTexMetr.Find(tstrTexName, cxTex, cyTex) && cxTex != 0 && cyTex != 0) {
					int i = 0;
					double fX[3], fY[3];
					if (iNVoS == NVOS_PXY || iNVoS == NVOS_NXY)
						for (; i < 3; i++)
							fX[i] = aTri.v[i].x + vtxMdlPos.x, fY[i] = aTri.v[i].y + vtxMdlPos.y;
					if (iNVoS == NVOS_PYZ || iNVoS == NVOS_NYZ)
						for (; i < 3; i++)
							fX[i] = aTri.v[i].y + vtxMdlPos.y, fY[i] = aTri.v[i].z + vtxMdlPos.z;
					if (iNVoS == NVOS_PXZ || iNVoS == NVOS_NXZ)
						for (; i < 3; i++)
							fX[i] = aTri.v[i].x + vtxMdlPos.x, fY[i] = aTri.v[i].z + vtxMdlPos.z;
					int iFlipH = (iNVoS == NVOS_NXY || iNVoS == NVOS_NYZ || iNVoS == NVOS_NXZ);
					int nBad = 0;
					if (fX[0] == fX[1] && fY[0] == fY[1]) nBad++;
					if (fX[1] == fX[2] && fY[1] == fY[2]) nBad++;
					if (fX[0] == fX[2] && fY[0] == fY[2]) nBad++;
					if (nBad > 0)
						continue;
					if (calcTexcAndReverse(
						cxTex, cyTex, 
						fX[0], fY[0], lpVtx3A->texcoord[0][0], lpVtx3A->texcoord[0][1], 
						fX[1], fY[1], lpVtx3B->texcoord[0][0], lpVtx3B->texcoord[0][1],
						fX[2], fY[2], lpVtx3C->texcoord[0][0], lpVtx3C->texcoord[0][1],
						iFlipH, 
						fTexcX, fTexcY, fTexAngle, fTexsX, fTexsY
						)
					) {
						bTexCoordReady = TRUE;
					}
				}
			}
		}
#if 0
		//	int nJoinPoly;
		//	do {
		//		nJoinPoly = 0;
		//		int iPoly1, iPoly2, nMyPolys = arrPolys.GetSize();
		//		for (iPoly1 = 0; iPoly1 + 1 < nMyPolys; iPoly1++) {
		//			for (iPoly2 = iPoly1 + 1; iPoly2 < nMyPolys; iPoly2++) {
		//				if (arrPolyNorms.ElementAt(iPoly1) != arrPolyNorms.ElementAt(iPoly2))
		//					continue;
		//				if (!arrPolys.ElementAt(iPoly1).JoinPoly(arrPolys.ElementAt(iPoly2)))
		//					continue;
		//				arrPolys.RemoveAt(iPoly2);
		//				arrPolyNorms.RemoveAt(iPoly2);
		//				arrPolys.ElementAt(iPoly1).Reduce();
		//				iPoly2--;
		//				nMyPolys--;
		//				nJoinPoly++;
		//			}
		//		}
		//	} while (nJoinPoly != 0);
#endif
		const nMyPolys = arrPolys.GetSize();
		iMax_PolysArray = __max(nMyPolys, iMax_PolysArray);
		for (int iMyPoly = 0; iMyPoly < nMyPolys; iMyPoly++, iIndex++) {
			myPoly_t &rPoly = arrPolys.ElementAt(iMyPoly);
			iMax_VertexArray = __max(iMax_VertexArray, rPoly.arrVerts.GetSize());
			rPoly.Reduce();
			if (rPoly.arrVerts.GetSize() < 3)
				continue;

#if 0
			//	{
			//		for (int i = 0, n = rPoly.arrVerts.GetSize(); i < n; i++) {
			//			tstrTemp.Format(
			//				_T("// %6.0f %6.0f %6.0f"), 
			//				(float)rPoly.arrVerts.GetAt(i).x,
			//				(float)rPoly.arrVerts.GetAt(i).y,
			//				(float)rPoly.arrVerts.GetAt(i).z
			//				);
			//			Write2(fileInto, tstrTemp);
			//		}
			//	}
#endif
			myVtx_t &rNorm = arrPolyNorms.ElementAt(iMyPoly);
			myVtx_t aNormThick(
				rNorm * -dThin
				);
#if 0
			//	{
			//		UINT cxTex, cyTex;
			//		if (rTexMetr.Find(tstrTexName, cxTex, cyTex) && cxTex != 0 && cyTex != 0) {
			//			BOOL bFin = FALSE;
			//			for (int iVt0 = 0; !bFin && iVt0 < rPoly.arrVerts.GetSize() - 1; iVt0++) {
			//				for (int iVt1 = iVt0 + 1; !bFin && iVt1 < rPoly.arrVerts.GetSize(); iVt1++) {
			//					const myVtx_t &v0 = rPoly.arrVerts.ElementAt(iVt0);
			//					const myVtx_t &v1 = rPoly.arrVerts.ElementAt(iVt1);
			//					for (int iRrV0 = iVertBase; !bFin && iRrV0 < iVertBase + nVerts; iRrV0++) {
			//						vertex3_t *lpRrV0 = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iRrV0));
			//						myVtx_t aRrV0(lpRrV0->position[0], lpRrV0->position[1], lpRrV0->position[2]);
			//						if (v0 != aRrV0)
			//							continue;
			//						for (int iRrV1 = iVertBase; !bFin && iRrV1 < iVertBase + nVerts; iRrV1++) {
			//							vertex3_t *lpRrV1 = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iRrV1));
			//							myVtx_t aRrV1(lpRrV1->position[0], lpRrV1->position[1], lpRrV1->position[2]);
			//							if (v1 != aRrV1)
			//								continue;
			//							double fX[2], fY[2];
			//							if (iNVoS == NVOS_XY)
			//								fX[0] = aRrV0.x, fY[0] = aRrV0.y,
			//								fX[1] = aRrV1.x, fY[1] = aRrV1.y;
			//							if (iNVoS == NVOS_YZ)
			//								fX[0] = aRrV0.y, fY[0] = aRrV0.z,
			//								fX[1] = aRrV1.y, fY[1] = aRrV1.z;
			//							if (iNVoS == NVOS_XZ)
			//								fX[0] = aRrV0.x, fY[0] = aRrV0.z,
			//								fX[1] = aRrV1.x, fY[1] = aRrV1.z;
			//							if (fX[0] == fX[1] || fY[0] == fY[1])
			//								continue;
			//							if (comp4(lpRrV0->texcoord[0][0], lpRrV1->texcoord[0][0]) == 0 || comp4(lpRrV0->texcoord[0][1], lpRrV1->texcoord[0][1]) == 0)
			//								continue;
			//							calcTexcAndReverse(
			//								cxTex, cyTex, 
			//								fX[0], fY[0], lpRrV0->texcoord[0][0], lpRrV0->texcoord[0][1], 
			//								fX[1], fY[1], lpRrV1->texcoord[0][0], lpRrV1->texcoord[0][1]
			//								);
			//							bFin = TRUE;
			//						}
			//					}
			//				}
			//			}
			//		}
			//	}
#endif

			tstrLine.Format(
				_T("// Brush #%d"), iIndex
				);
			Write2(fileInto, tstrLine);
			Write2(fileInto, _T("{"));
			// ->
			const nMyTris = rPoly.arrVerts.GetSize();
			ASSERT(nMyTris >= 3);
			int nHatch = (nMyTris > 3) ? 2 : 1;
			for (int iMyTri = 0; iMyTri < nMyTris + nHatch * 2; iMyTri++) {
				myTri_t aTri;
				if (iMyTri < nHatch) {
					const iLip = iMyTri;
					aTri.v[0] = rPoly.arrVerts.ElementAt((2 * iLip + 0) % nMyTris);
					aTri.v[1] = rPoly.arrVerts.ElementAt((2 * iLip + 1) % nMyTris);
					aTri.v[2] = rPoly.arrVerts.ElementAt((2 * iLip + 2) % nMyTris);
				} else if (iMyTri < nHatch * 2) {
					const iLip = iMyTri - nHatch;
					aTri.v[0] = rPoly.arrVerts.ElementAt((nMyTris - 0 - 2 * iLip) % nMyTris) + aNormThick;
					aTri.v[1] = rPoly.arrVerts.ElementAt((nMyTris - 1 - 2 * iLip) % nMyTris) + aNormThick;
					aTri.v[2] = rPoly.arrVerts.ElementAt((nMyTris - 2 - 2 * iLip) % nMyTris) + aNormThick;
				} else {
					int i0 = (iMyTri - nHatch * 2 + 0);
					int i1 = (iMyTri - nHatch * 2 + 1) % nMyTris;
					aTri.v[0] = rPoly.arrVerts.ElementAt(i0);
					aTri.v[1] = rPoly.arrVerts.ElementAt(i1) + aNormThick;
					aTri.v[2] = rPoly.arrVerts.ElementAt(i1);
				}

				tstrPerTri.Empty();
				for (int iMyVert = 0; iMyVert < 3; iMyVert++) {
					aTri.v[iMyVert] += vtxMdlPos;
					double fx = aTri.v[iMyVert].x;
					double fy = aTri.v[iMyVert].y;
					double fz = aTri.v[iMyVert].z;
					tstrTemp.Format(_T("( %.0f %.0f %.0f ) "), (float)fx, (float)fy, (float)fz);
					tstrPerTri += tstrTemp;
				}
				tstrLine = tstrPerTri;
				// ( 640 640 128 ) ( 128 640 128 ) ( 128 128 128 ) test/usweapons 0 0 0.00 1 1 0 0 0 
				{
					tstrTemp.Format(
						_T("%s %.0lf %.0lf %.2lf %lf %lf %ld %ld %ld%s"),
						tstrTexNamec,
						(double)fTexcX, (double)fTexcY, (double)fTexAngle, (double)fTexsX, (double)fTexsY,
						(DWORD)0, (DWORD)0, (DWORD)0,
						tstrTexAdd
					);
				}
				tstrLine+= tstrTemp;
				Write2(fileInto, tstrLine);
			}
			// <-
			Write2(fileInto, _T("}"));
		}
	}
	// 
	if (bAddBBox) {
		tstrLine.Format(
			_T("// Brush #%d"), iIndex
			);
		Write2(fileInto, tstrLine);
		Write2(fileInto, _T("{"));
		// ->
		myVtx_t v8[8] = {
			myVtx_t(lpMdl3->mins[0], lpMdl3->mins[1], lpMdl3->mins[2]),
			myVtx_t(lpMdl3->maxs[0], lpMdl3->mins[1], lpMdl3->mins[2]),
			myVtx_t(lpMdl3->mins[0], lpMdl3->maxs[1], lpMdl3->mins[2]),
			myVtx_t(lpMdl3->maxs[0], lpMdl3->maxs[1], lpMdl3->mins[2]),
			myVtx_t(lpMdl3->mins[0], lpMdl3->mins[1], lpMdl3->maxs[2]),
			myVtx_t(lpMdl3->maxs[0], lpMdl3->mins[1], lpMdl3->maxs[2]),
			myVtx_t(lpMdl3->mins[0], lpMdl3->maxs[1], lpMdl3->maxs[2]),
			myVtx_t(lpMdl3->maxs[0], lpMdl3->maxs[1], lpMdl3->maxs[2]),
		};
		static iTriTbl[6][3] = {
			{4,6,5}, // 上蓋
			{1,3,0}, // 下蓋
			{1,5,3}, // 東
			{2,6,0}, // 西
			{0,4,1}, // 南
			{3,7,2}, // 北
		};
		for (int iMyTri = 0; iMyTri < 6; iMyTri++) {
			myTri_t aTri(
				v8[iTriTbl[iMyTri][0]],
				v8[iTriTbl[iMyTri][1]],
				v8[iTriTbl[iMyTri][2]]
			);
			tstrPerTri.Empty();
			for (int iMyVert = 0; iMyVert < 3; iMyVert++) {
				aTri.v[iMyVert] += vtxMdlPos;
				double fx = aTri.v[iMyVert].x;
				double fy = aTri.v[iMyVert].y;
				double fz = aTri.v[iMyVert].z;
				tstrTemp.Format(_T("( %.0f %.0f %.0f ) "), (float)fx, (float)fy, (float)fz);
				tstrPerTri += tstrTemp;
			}
			tstrLine = tstrPerTri;
			tstrTemp.Format(
				_T("%s %.0lf %.0lf %.2lf %lf %lf %ld %ld %ld"),
				_T("common/trigger"),
				(double)0, (double)0, (double)0, (double)1, (double)1,
				(DWORD)0, (DWORD)0, (DWORD)0
				);
			tstrLine += tstrTemp;
			Write2(fileInto, tstrLine);
		}
		// <-
		Write2(fileInto, _T(" }"));
		iIndex++;
	}
	ASSERT(iModelNo != 0 || (iModelNo == 0 && nBadNorm == 0 && nBadTri == 0));

	return TRUE;
}

BOOL ReparseTextLines2Map(
	CFile &fileInto1,
	MB2MMLOptical &rOpts,
	const CByteArray &arrTexture,
	const CByteArray &arrFace,
	const CByteArray &arrVertex,
	const CByteArray &arrModel,
	const CByteArray &arrEntity,
	const CByteArray &arrMeshVerts,
	const CByteArray &arrTerr,
	const CStringList &l
)
{
	POSITION p = l.GetHeadPosition();
	CString tstr, tstrKey, tstrVal, tstrComment, tstrClassName;
	int iIndent = 0, iEntityIndex = 0, iExtractModelNo = -1;
	myVtx_t vtxMdlOrigin;
	while (p) {
		tstr = l.GetNext(p);
		if (tstr.IsEmpty())
			continue;
		if (iIndent == 0) {
			if (tstr[0] == TCHAR('{')) {
				if (iIndent == 1)
					return FALSE;
				iIndent++;

				tstrComment.Format(
					_T("// Entity #%d"), iEntityIndex
					);
				Write2(fileInto1, tstrComment);
				Write2(fileInto1, tstr);
			}
		} else {
			if (tstr[0] == TCHAR('}')) {
				if (iIndent == 0)
					return FALSE;
				iIndent--;

				BOOL bAddBBox = 0;
				if (rOpts.nMask & MB2MML_F_ENABLE_BBOX_TREATS) {
					for (int iIdx = 0; iIdx < rOpts.arrBBoxTreats.GetSize(); iIdx++) {
						if (rOpts.arrBBoxTreats.ElementAt(iIdx) == tstrClassName) {
							bAddBBox = TRUE; break;
						}
					}
				}

				if (iExtractModelNo != -1) {
					CFindTexMetrExp aTexMetr(rOpts.arrTexMetr);
					WriteModel2(fileInto1, iExtractModelNo, vtxMdlOrigin, bAddBBox, aTexMetr, arrTexture, arrFace, arrVertex, arrModel, arrMeshVerts, arrTerr);
					iExtractModelNo = -1;
				}
				iEntityIndex++;
				vtxMdlOrigin = myVtx_t();
				tstrClassName.Empty();

				Write2(fileInto1, tstr);
			} else {
				if (ParseKeyNVal(tstr, tstrKey, tstrVal)) {
					if (tstrKey == _T("model") && !tstrVal.IsEmpty() && tstrVal[0] == TCHAR('*')) {
						// キタ～
						iExtractModelNo = _ttoi((LPCTSTR)tstrVal + 1);
					} else {
						Write2(fileInto1, tstr);

						if (tstrKey == _T("classname") && tstrVal == _T("worldspawn")) {
							// キタ～
							iExtractModelNo = 0;
						} else if (tstrKey == _T("origin")) {
							// キタ～
							double fx, fy, fz;
							if (_stscanf(tstrVal, _T("%lf %lf %lf"), &fx, &fy, &fz) != 3) {
								fx = fy = fz = 0;
							}
							vtxMdlOrigin = myVtx_t(fx, fy, fz);
						}
						if (tstrKey == _T("classname")) {
							tstrClassName = tstrVal;
						}
					}
				}
			}
		}
	}
	return TRUE;
}

BOOL ReadDirEnt(CFile &file, int iEntryNo, CByteArray &chunk)
{
	TRY
		file.Seek(12 + iEntryNo * 8, CFile::begin);
		DWORD cb = file.GetLength();
		DWORD dw[2] = {0, 0};
		if (file.Read((void *)&dw[0], 8) != 8)
			AfxThrowFileException(CFileException::endOfFile);
		if (dw[0] >= cb || dw[0] + dw[1] > cb)
			AfxThrowFileException(CFileException::endOfFile);
		file.Seek(dw[0], CFile::begin);
		chunk.SetSize(dw[1]);
		if (file.Read((void *)chunk.GetData(), dw[1]) != dw[1])
			AfxThrowFileException(CFileException::endOfFile);
		return TRUE;
	CATCH(CFileException, e)
		chunk.SetSize(0);
	END_CATCH

	return FALSE;
}

CRect calcInnerRect(CRect rc, CSize size)
{
	float fAr = (float)rc.Width() / (float)rc.Height();
	float fAs = (float)size.cx / (float)size.cy;
	if (fAr < fAs) {
		// 横長
		int cx = (int)(rc.Width());
		int cy = (int)(rc.Width() / fAs);
		int ty = rc.CenterPoint().y;
		return CRect(rc.left, ty - cy / 2, rc.right, ty + cy / 2);
	} else {
		// 縦長 or 他
		int cx = (int)(rc.Height() * fAs);
		int cy = (int)(rc.Height());
		int tx = rc.CenterPoint().x;
		return CRect(tx - cx / 2, rc.top, tx + cx / 2, rc.bottom);
	}
}

namespace MohBSPToWinDib
{

struct MYCDFVP
{
	POINTS p[2];
};

struct MYCDFMP
{
	int iIdx;
	POINTS o;

	MYCDFMP() { iIdx = 0; o.x = o.y = 0; }
};

inline bool operator ==(const POINTS &s1, const POINTS &s2)
{
	return (s1.x == s2.x && s1.y == s2.y) ? true : false;
}

inline bool operator !=(const POINTS &s1, const POINTS &s2)
{
	return (s1.x != s2.x || s1.y != s2.y) ? true : false;
}

struct MY2DTRI_T
{
	POINTS p[3];
};

class MY2DPOLY_T
{
public:
	// 
	CArray<POINTS, POINTS> arrVerts;
	// 
	MY2DPOLY_T()
	{
		arrVerts.SetSize(0, 64);
	}
	// 
	MY2DPOLY_T &operator =(const MY2DPOLY_T &s)
	{
		arrVerts.Copy(s.arrVerts);
		return *this;
	}
	// 
	void RemoveAll()
	{
		arrVerts.RemoveAll();
	}
	// 
	BOOL JoinTri(MY2DTRI_T &r)
	{
		int iVert, nVerts = arrVerts.GetSize();
		if (nVerts == 0) {
			arrVerts.SetSize(3);
			arrVerts[0] = r.p[0],
			arrVerts[1] = r.p[1],
			arrVerts[2] = r.p[2];
			return TRUE;
		}
		ASSERT(nVerts >= 3);
		for (iVert = 0; iVert < nVerts; iVert++) {
			const iVert0 = (iVert    ) % nVerts;
			const iVert1 = (iVert + 1) % nVerts;
			const POINTS &rV0 = arrVerts.ElementAt(iVert0);
			const POINTS &rV1 = arrVerts.ElementAt(iVert1);
			ASSERT(rV0 != rV1);
			if (rV0 == r.p[1] && rV1 == r.p[0]) {
				ASSERT(rV0 != r.p[2] && rV1 != r.p[2]);
				arrVerts.InsertAt(iVert + 1, r.p[2]);
				return TRUE;
			}
			if (rV0 == r.p[2] && rV1 == r.p[1]) {
				ASSERT(rV0 != r.p[0] && rV1 != r.p[0]);
				arrVerts.InsertAt(iVert + 1, r.p[0]);
				return TRUE;
			}
			if (rV0 == r.p[0] && rV1 == r.p[2]) {
				ASSERT(rV0 != r.p[1] && rV1 != r.p[1]);
				arrVerts.InsertAt(iVert + 1, r.p[1]);
				return TRUE;
			}
		}
		return FALSE;
	}
	// 
	void Reduce()
	{
		int nRepair;
		do {
			nRepair = 0;
			int i, n = arrVerts.GetSize();
			for (i = 0; i < n && n >= 2; i++) {
				if (arrVerts[(i + 0)] == arrVerts[(i + 1) % n]) {
					arrVerts.RemoveAt((i + 1) % n);
					i--;
					n--;
					nRepair++;
				}
			}
			for (i = 0; i < n && n >= 3; i++) {
				if (arrVerts[(i + 0)] == arrVerts[(i + 2) % n]) {
					arrVerts.RemoveAt((i + 1) % n);
					if (!(i < (i + 1) % n))
						i--;
					n--;
					nRepair++;
					arrVerts.RemoveAt((i + 0));
					i--;
					n--;
					nRepair++;
				}
			}
		} while (nRepair != 0);
	}

};

BOOL ParseKEYnVAL(
	const char *&lpcPos,
	const char *lpcEnd,
	CByteArray &arrKey,
	CByteArray &arrVal
	)
{
	for (; lpcPos != lpcEnd && lpcPos[0] != '\"'; lpcPos++);
	if (lpcPos != lpcEnd) lpcPos++;
	const char *lpcPtr = lpcPos;
	for (; lpcPos != lpcEnd && lpcPos[0] != '\"'; lpcPos++);
	arrKey.SetSize(lpcPos - lpcPtr + 1);
	arrKey[lpcPos - lpcPtr] = 0;
	strncpy((char *)arrKey.GetData(), lpcPtr, lpcPos - lpcPtr);
	if (lpcPos != lpcEnd) lpcPos++;

	for (; lpcPos != lpcEnd && lpcPos[0] != '\"'; lpcPos++);
	if (lpcPos != lpcEnd) lpcPos++;
	lpcPtr = lpcPos;
	for (; lpcPos != lpcEnd && lpcPos[0] != '\"'; lpcPos++);
	arrVal.SetSize(lpcPos - lpcPtr + 1);
	arrVal[lpcPos - lpcPtr] = 0;
	strncpy((char *)arrVal.GetData(), lpcPtr, lpcPos - lpcPtr);
	if (lpcPos != lpcEnd) lpcPos++;

	return TRUE;
}

BOOL CreateDibFrom(
	MB2MMLOptical &rOpts,
	const CByteArray &arrFace,
	const CByteArray &arrVertex,
	const CByteArray &arrModel,
	const CByteArray &arrEntity,
	const CByteArray &arrMeshVerts,
	const CByteArray &arrTerr
)
{
	CArray<MYCDFVP, MYCDFVP &> arrVtx;
	CArray<MYCDFMP, MYCDFMP &> arrMdl;
	CArray<MY2DPOLY_T, MY2DPOLY_T &> arrPolys;
	arrVtx.SetSize(0, 10240);
	arrMdl.SetSize(1, 128);
	arrPolys.SetSize(0, 16);

	LPCSTR lpcPos = (LPCSTR)arrEntity.GetData();
	LPCSTR lpcEnd = lpcPos + arrEntity.GetSize();
	TRY
		int iLv = 0, iEntIdx = 0;
		MYCDFMP aMP;
		aMP.iIdx = 0;
		aMP.o.x = aMP.o.y = 0;
		arrMdl[0] = aMP;
		aMP.iIdx = -1;
		CByteArray arrKey, arrVal;
		while (lpcPos != lpcEnd) {
			LPCSTR lpcParaPos = lpcPos;
			LPCSTR lpcParaEnd = lpcPos;
			for (; lpcParaEnd != lpcEnd && lpcParaEnd[0] != '\n' && lpcParaEnd[0] != '\r'; lpcParaEnd++);

			for (; lpcParaPos != lpcParaEnd && isspace((BYTE)*lpcParaPos); lpcParaPos++);

			while (lpcParaPos != lpcParaEnd) {
				if (lpcParaPos[0] == '/' && &lpcParaPos[1] != lpcParaEnd && lpcParaPos[1] == '/') {
					break;
				}
				if (lpcParaPos[0] == '{') {
					if (iLv == 2)
						AfxThrowArchiveException(CArchiveException::badIndex);
					iLv++;
					lpcParaPos++;
				} else if (lpcParaPos[0] == '}') {
					if (iLv == 0)
						AfxThrowArchiveException(CArchiveException::badIndex);
					iLv--;
					if (iLv == 0) {
						if (!(aMP.iIdx < 0)) {
							arrMdl.Add(aMP);
						}
						aMP.o.x = aMP.o.y = 0;
						aMP.iIdx = -1;

						iEntIdx++;
					}
					lpcParaPos++;
				} else if (lpcParaPos[0] == '\"') {
					if (iLv != 1)
						AfxThrowArchiveException(CArchiveException::badIndex);
					if (ParseKEYnVAL(lpcParaPos, lpcParaEnd, arrKey, arrVal)) {
						if (stricmp("origin", (char *)arrKey.GetData()) == 0) {
							int x, y, z;
							if (sscanf((char *)arrVal.GetData(), "%d %d %d", &x, &y, &z) == 3) {
								aMP.o.x = x;
								aMP.o.y = y;
							}
						} else if (stricmp("model", (char *)arrKey.GetData()) == 0) {
							int i;
							if (sscanf((char *)arrVal.GetData(), "*%d", &i) == 1 && aMP.iIdx < 0) {
								aMP.iIdx = i;
							}
						}
					}
				} else {
					lpcParaPos++;
				}

				for (; lpcParaPos != lpcParaEnd && isspace((BYTE)*lpcParaPos); lpcParaPos++);
			}

			lpcPos = lpcParaEnd;
			if (lpcPos == lpcEnd)
				break;
			if (lpcPos[0] == '\n') {
				lpcPos++;
			} else if (lpcPos[0] == '\r') {
				lpcPos++;
				if (lpcPos != lpcEnd && lpcPos[0] == '\n') {
					lpcPos++;
				}
			}
		}

		const nMdls = arrModel.GetSize() / 40;
		const nVtxs = arrVertex.GetSize() / 44;

		POINT aPtMin = {SHRT_MAX, SHRT_MAX}, aPtMax = {SHRT_MIN, SHRT_MIN};
		for (int iMdlE = 0; iMdlE < arrMdl.GetSize(); iMdlE++) {
			MYCDFMP &r = arrMdl[iMdlE];
			if (r.iIdx < 0 || nMdls <= r.iIdx)
				AfxThrowArchiveException(CArchiveException::badIndex);
			const model3_t *lpMdl3 = reinterpret_cast<const model3_t *>(arrModel.GetData() + 40 * r.iIdx);
			const UINT iSurfMin = lpMdl3->face;
			const UINT iSurfMax = lpMdl3->face + lpMdl3->n_faces;
			UINT iSurf;
			for (iSurf = iSurfMin; iSurf < iSurfMax; iSurf++) {
				const face3_t *lpSurf3 = reinterpret_cast<const face3_t *>(arrFace.GetData() + 108 * iSurf);
				if (lpSurf3->type == 1) {
					const UINT iMVertMin = lpSurf3->meshvert;
					const UINT iMVertMax = lpSurf3->meshvert + lpSurf3->n_meshverts;
					const DWORD *lpMVert = reinterpret_cast<const DWORD *>(arrMeshVerts.GetData());
					UINT iMVert;
					const iVertPos = lpSurf3->vertex;
					ASSERT(iVertPos < nVtxs);
					UINT iIdx;
					POINTS p0, p1;
					arrPolys.RemoveAll();
					for (iMVert = iMVertMin, iIdx = 0; iMVert < iMVertMax; iMVert += 3, iIdx += 3) {
						ASSERT(iMVert + 3 <= iMVertMax);
						const UINT iVertIdx[] = {
							iVertPos + lpMVert[iMVert + 0],
							iVertPos + lpMVert[iMVert + 1],
							iVertPos + lpMVert[iMVert + 2],
						};
						const vertex3_t *lpVtx3[] = {
							reinterpret_cast<const vertex3_t *>(arrVertex.GetData() + 44 * iVertIdx[0]),
							reinterpret_cast<const vertex3_t *>(arrVertex.GetData() + 44 * iVertIdx[1]),
							reinterpret_cast<const vertex3_t *>(arrVertex.GetData() + 44 * iVertIdx[2])
						};
						MY2DTRI_T v3;
						v3.p[0].x = (SHORT)lpVtx3[0]->position[0] + r.o.x;
						v3.p[0].y = (SHORT)lpVtx3[0]->position[1] + r.o.y;
						v3.p[1].x = (SHORT)lpVtx3[1]->position[0] + r.o.x;
						v3.p[1].y = (SHORT)lpVtx3[1]->position[1] + r.o.y;
						v3.p[2].x = (SHORT)lpVtx3[2]->position[0] + r.o.x;
						v3.p[2].y = (SHORT)lpVtx3[2]->position[1] + r.o.y;
						if (v3.p[0] != v3.p[1] && v3.p[1] != v3.p[2] && v3.p[2] != v3.p[0]) {
							for (int x = 0; x < arrPolys.GetSize(); x++) {
								if (arrPolys[x].JoinTri(v3)) break;
							}
							if (x == arrPolys.GetSize()) {
								MY2DPOLY_T o;
								o.JoinTri(v3);
								arrPolys.Add(o);
							}
						}
					}
					for (int iPoly = 0; iPoly < arrPolys.GetSize(); iPoly++) {
						MY2DPOLY_T &o = arrPolys[iPoly];
						o.Reduce();

						for (int iIdx = 0; iIdx < o.arrVerts.GetSize(); iIdx++) {
							POINTS &r = o.arrVerts[iIdx];
							if (iIdx == 0) {
								p0 = r;
							} else {
								MYCDFVP aVP;
								aVP.p[0] = p1;
								aVP.p[1] = r;
								arrVtx.Add(aVP);
							}
							p1 = r;
						}
						MYCDFVP aVP;
						aVP.p[0] = p1;
						aVP.p[1] = p0;
						arrVtx.Add(aVP);
					}
				} else if (lpSurf3->type == 2) {
					const iVertMin = lpSurf3->vertex;
					const iVertMax = iVertMin + lpSurf3->n_vertexes;
					const cyMtx = lpSurf3->size[0];
					const cxMtx = lpSurf3->size[1];
					ASSERT(iVertMin < nVtxs && iVertMax < nVtxs);
					POINTS p0, p1;
					for (int yMtx = 0; yMtx < cyMtx; yMtx++) {
						POINTS p;
						for (int xMtx = 0; xMtx < cxMtx; xMtx++) {
							const iVert = cyMtx * xMtx + yMtx;
							ASSERT(0 <= iVert && iVert < nVtxs);
							vertex3_t *lpVtx3 = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertMin + iVert));
							p.x = (SHORT)lpVtx3->position[0] + r.o.x;
							p.y = (SHORT)lpVtx3->position[1] + r.o.y;
							if (xMtx == 0) {
								p0 = p;
							} else {
								MYCDFVP aVP;
								aVP.p[0] = p1;
								aVP.p[1] = p;
								arrVtx.Add(aVP);
							}
							p1 = p;
						}
						MYCDFVP aVP;
						aVP.p[0] = p1;
						aVP.p[1] = p0;
						arrVtx.Add(aVP);
					}
					for (int xMtx = 0; xMtx < cxMtx; xMtx++) {
						POINTS p;
						for (int yMtx = 0; yMtx < cyMtx; yMtx++) {
							const iVert = cyMtx * xMtx + yMtx;
							ASSERT(0 <= iVert && iVert < nVtxs);
							vertex3_t *lpVtx3 = (vertex3_t *)(((const BYTE *)arrVertex.GetData()) + 44 * (iVertMin + iVert));
							p.x = (SHORT)lpVtx3->position[0] + r.o.x;
							p.y = (SHORT)lpVtx3->position[1] + r.o.y;
							if (xMtx == 0) {
								p0 = p;
							} else {
								MYCDFVP aVP;
								aVP.p[0] = p1;
								aVP.p[1] = p;
								arrVtx.Add(aVP);
							}
							p1 = p;
						}
						MYCDFVP aVP;
						aVP.p[0] = p1;
						aVP.p[1] = p0;
						arrVtx.Add(aVP);
					}
				}
			}
		}

		{
			const nTerrs = arrTerr.GetSize() / 388;
			for (int iTerr = 0; iTerr < nTerrs; iTerr++) {
				const BYTE *lpbEach = arrTerr.GetData() + 388 * iTerr;
				int x = (CHAR)lpbEach[0x24] * (int)64;
				int y = (CHAR)lpbEach[0x25] * (int)64;

				static const POINTS aPts[] = {
					{  0,  0},{512,  0},
					{  0, 64},{512, 64},
					{  0,128},{512,128},
					{  0,192},{512,192},
					{  0,256},{512,256},
					{  0,320},{512,320},
					{  0,384},{512,384},
					{  0,448},{512,448},
					{  0,512},{512,512},

					{  0,  0},{  0,512},
					{ 64,  0},{ 64,512},
					{128,  0},{128,512},
					{192,  0},{192,512},
					{256,  0},{256,512},
					{320,  0},{320,512},
					{384,  0},{384,512},
					{448,  0},{448,512},
					{512,  0},{512,512},

					{384,  0},{512,128},
					{256,  0},{512,256},
					{128,  0},{512,384},
					{  0,  0},{512,512},
					{  0,128},{384,512},
					{  0,256},{256,512},
					{  0,384},{128,512},

					{512,384},{384,512},
					{512,256},{256,512},
					{512,128},{128,512},
					{512,  0},{  0,512},
					{384,  0},{  0,384},
					{256,  0},{  0,256},
					{128,  0},{  0,128},
				};
				const nPts = sizeof(aPts) / sizeof(POINTS);
				for (int iPt = 0; iPt < nPts; iPt += 2) {
					MYCDFVP aVP;
					aVP.p[0].x = x + aPts[iPt + 0].x;
					aVP.p[0].y = y + aPts[iPt + 0].y;
					aVP.p[1].x = x + aPts[iPt + 1].x;
					aVP.p[1].y = y + aPts[iPt + 1].y;
					arrVtx.Add(aVP);
				}
			}
		}

		for (int iMidIdx = 0; iMidIdx < arrVtx.GetSize(); iMidIdx++) {
			for (int i = 0; i < 2; i++)
				aPtMin.x = __min(aPtMin.x, arrVtx[iMidIdx].p[i].x),
				aPtMin.y = __min(aPtMin.y, arrVtx[iMidIdx].p[i].y),
				aPtMax.x = __max(aPtMax.x, arrVtx[iMidIdx].p[i].x),
				aPtMax.y = __max(aPtMax.y, arrVtx[iMidIdx].p[i].y);
		}

		{
			SIZE size = rOpts.sizeMakeBM;
			POINT aPtCent = {size.cx / 2, size.cy / 2};

			POINT aPtInc = {
				((int)(aPtMax.x - aPtMin.x) * 5) / 100,
				((int)(aPtMax.y - aPtMin.y) * 5) / 100
			};
			aPtMin.x -= aPtInc.x;
			aPtMin.y -= aPtInc.y;
			aPtMax.x += aPtInc.x;
			aPtMax.y += aPtInc.y;

			CRect rcIr = calcInnerRect(CRect(0, 0, size.cx, size.cy), CSize((int)aPtMax.x - aPtMin.x, (int)aPtMax.y - aPtMin.y));

			int iXorg = ((int)aPtMax.x + aPtMin.x) / 2;
			int iYorg = ((int)aPtMax.y + aPtMin.y) / 2;
			float fXs = rcIr.Width() / (float)(aPtMax.x - aPtMin.x);
			float fYs = rcIr.Height()/ (float)(aPtMax.y - aPtMin.y);

			BITMAPINFOHEADER aBIH;
			ZeroMemory(&aBIH, sizeof(aBIH));
			aBIH.biSize = sizeof(aBIH),
			aBIH.biWidth = size.cx,
			aBIH.biHeight = size.cy,
			aBIH.biPlanes = 1,
			aBIH.biBitCount = 24,
			aBIH.biCompression = BI_RGB;
			HBITMAP hBM = CreateDIBSection(
				NULL,
				(BITMAPINFO *)&aBIH,
				DIB_RGB_COLORS,
				NULL,
				NULL,
				0
				);
			if (hBM != NULL) {
				CDC aDC;
				if (aDC.CreateCompatibleDC(NULL)) {
					HGDIOBJ hOldBM = aDC.SelectObject(hBM);
					aDC.FillSolidRect(0, 0, size.cx, size.cy, RGB(0, 0, 0));
					aDC.FillSolidRect(rcIr, RGB(30, 30, 30));
					CPen aPen(PS_SOLID, 0, RGB(255, 255, 255));
					HGDIOBJ hOldPen = aDC.SelectObject(aPen.m_hObject);

					for (int i = 0; i < arrVtx.GetSize(); i++) {
						MYCDFVP &r = arrVtx[i];
						aDC.MoveTo(
							(int)(( r.p[0].x - iXorg) * fXs + aPtCent.x),
							(int)((-r.p[0].y + iYorg) * fYs + aPtCent.y)
							);
						aDC.LineTo(
							(int)(( r.p[1].x - iXorg) * fXs + aPtCent.x),
							(int)((-r.p[1].y + iYorg) * fYs + aPtCent.y)
							);
					}

					aDC.SelectObject(hOldPen);
					aDC.SelectObject(hOldBM);
				}
				rOpts.hMadeBM = hBM;
			}
		}
		return TRUE;
	CATCH_ALL(e)

	END_CATCH_ALL

	return FALSE;
}

};

BOOL CreateMapDib(
	MB2MMLOptical &rOpts,
	const CByteArray &arrTexture,
	const CByteArray &arrFace,
	const CByteArray &arrVertex,
	const CByteArray &arrModel,
	const CByteArray &arrEntity,
	const CByteArray &arrMeshVerts
)
{
	SIZE sizeNewBM = rOpts.sizeMakeBM;
	BITMAPINFOHEADER aBih;
	ZeroMemory(&aBih, sizeof(aBih));
	aBih.biSize = sizeof(BITMAPINFOHEADER);
	aBih.biWidth = sizeNewBM.cx;
	aBih.biHeight = sizeNewBM.cy;
	aBih.biPlanes = 1;
	aBih.biBitCount = 24;
	aBih.biCompression = BI_RGB;
	HBITMAP hBM = CreateDIBSection(
		NULL,
		(BITMAPINFO *)&aBih,
		DIB_RGB_COLORS,
		NULL,
		NULL,
		0
		);
	ASSERT(hBM);
	CPoint aPtCenter(sizeNewBM.cx / 2, sizeNewBM.cy / 2);
	{
		float fMinX = FLT_MAX, fMinY = FLT_MAX, fMaxX = FLT_MIN, fMaxY = FLT_MIN;

		UINT iMdl;
		const UINT nMdls = arrModel.GetSize() / 40;
		const UINT nVtxs = arrVertex.GetSize() / 44;
		// 1
		for (iMdl = 0; iMdl < nMdls; iMdl++) {
			const model3_t *lpMdl3 = reinterpret_cast<const model3_t *>(arrModel.GetData() + 40 * iMdl);
			const UINT iSurfMin = lpMdl3->face;
			const UINT iSurfMax = lpMdl3->face + lpMdl3->n_faces;
			UINT iSurf;
			for (iSurf = iSurfMin; iSurf < iSurfMax; iSurf++) {
				const face3_t *lpSurf3 = reinterpret_cast<const face3_t *>(arrFace.GetData() + 108 * iSurf);
				if (lpSurf3->type == 1) {
					const UINT iMVertMin = lpSurf3->meshvert;
					const UINT iMVertMax = lpSurf3->meshvert + lpSurf3->n_meshverts;
					const DWORD *lpMVert = reinterpret_cast<const DWORD *>(arrMeshVerts.GetData());
					UINT iMVert;
					const UINT iVertPos = lpSurf3->vertex;
					ASSERT(iVertPos < nVtxs);
					UINT iIdx;
					for (iMVert = iMVertMin, iIdx = 0; iMVert < iMVertMax; iMVert++, iIdx++) {
						const UINT iMVertIdx = lpMVert[iMVert];
						const UINT iVertIdx = iVertPos + iMVertIdx;
						const vertex3_t *lpVtx3 = reinterpret_cast<const vertex3_t *>(arrVertex.GetData() + 44 * iVertIdx);
						float fX = lpVtx3->position[0], fY = -lpVtx3->position[1];
						if (fX < fMinX)
							fMinX = fX;
						if (fX > fMaxX)
							fMaxX = fX;

						if (fY < fMinY)
							fMinY = fY;
						if (fY > fMaxY)
							fMaxY = fY;
					}
				}
			}
		}

		double fXorg, fYorg, fXs, fYs;
		{
			float fXp = (float)((fMaxX - fMinX) * 0.05f);
			float fYp = (float)((fMaxY - fMinY) * 0.05f);
			fMinX -= fXp;
			fMinY -= fYp;
			fMaxX += fXp;
			fMaxY += fYp;
			float fcX = (fMaxX - fMinX), fcY = (fMaxY - fMinY);
			CRect rcBM(0, 0, sizeNewBM.cx, sizeNewBM.cy);
			CRect rcIr = calcInnerRect(rcBM, CSize((int)fcX, (int)fcY));

			CBrush aBr(RGB(32, 32, 32));
			CDC aDC;
			aDC.CreateCompatibleDC(NULL);
			HGDIOBJ hOldObj = aDC.SelectObject((HGDIOBJ)hBM);
			aDC.FillRect(rcIr, &aBr);
			aDC.SelectObject(hOldObj);

			fXorg = (fMaxX + fMinX) / 2;
			fYorg = (fMaxY + fMinY) / 2;
			fXs = rcIr.Width() / fcX;
			fYs = rcIr.Height() / fcY;
		}
		// 2
		CDC aDC;
		aDC.CreateCompatibleDC(NULL);
		HGDIOBJ hOldObj1 = aDC.SelectObject(hBM);
		CPen aPen(PS_SOLID, 0, RGB(255, 255, 255));
		HGDIOBJ hOldObj2 = aDC.SelectObject(aPen.m_hObject);
		for (iMdl = 0; iMdl < nMdls; iMdl++) {
			const model3_t *lpMdl3 = reinterpret_cast<const model3_t *>(arrModel.GetData() + 40 * iMdl);
			const UINT iSurfMin = lpMdl3->face;
			const UINT iSurfMax = lpMdl3->face + lpMdl3->n_faces;
			UINT iSurf;
			for (iSurf = iSurfMin; iSurf < iSurfMax; iSurf++) {
				const face3_t *lpSurf3 = reinterpret_cast<const face3_t *>(arrFace.GetData() + 108 * iSurf);
				if (lpSurf3->type == 1) {
					const UINT iMVertMin = lpSurf3->meshvert;
					const UINT iMVertMax = lpSurf3->meshvert + lpSurf3->n_meshverts;
					const DWORD *lpMVert = reinterpret_cast<const DWORD *>(arrMeshVerts.GetData());
					UINT iMVert;
					const UINT iVertPos = lpSurf3->vertex;
					UINT iIdx;
					int iOx, iOy;
					for (iMVert = iMVertMin, iIdx = 0; iMVert < iMVertMax; iMVert++, iIdx++) {
						const UINT iMVertIdx = lpMVert[iMVert];
						const UINT iVertIdx = iVertPos + iMVertIdx;
						const vertex3_t *lpVtx3 = reinterpret_cast<const vertex3_t *>(arrVertex.GetData() + 44 * iVertIdx);
						float fX = lpVtx3->position[0], fY = -lpVtx3->position[1];

						int iVx = (int)((fX - fXorg) * fXs) + aPtCenter.x;
						int iVy = (int)((fY - fYorg) * fYs) + aPtCenter.y;

						if (iIdx == 0)
							aDC.MoveTo(iOx = iVx, iOy = iVy);
						else
							aDC.LineTo(iVx, iVy);
					}
					aDC.LineTo(iOx, iOy);
				}
			}
		}
		aDC.SelectObject(hOldObj1);
		aDC.SelectObject(hOldObj2);
	}
	rOpts.hMadeBM = hBM;
	return TRUE;
}

void MyTest()
{

}

};

int _cdecl _MohBSP2Map(LPCTSTR lpszIn1, LPCTSTR lpszInto1, MB2MMLOptical &rOpts)
{
	TRY
#ifdef _DEBUG
		MyTest();
#endif

		CFile fileIn1(lpszIn1, CFile::modeRead | CFile::shareDenyWrite);
		CFile fileInto1(lpszInto1, CFile::modeCreate | CFile::modeNoTruncate | CFile::modeReadWrite | CFile::shareDenyWrite);
		CByteArray arrTexture, arrFace, arrVertex, arrModel, arrEntity, arrMeshVerts, arrTerr;
		if (g_bVerbose) _putts(_T("Reading Textures entry."));
		if (!ReadDirEnt(fileIn1, 0, arrTexture))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Faces entry."));
		if (!ReadDirEnt(fileIn1, 3, arrFace))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Vertices entry."));
		if (!ReadDirEnt(fileIn1, 4, arrVertex))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Mesh vertices entry."));
		if (!ReadDirEnt(fileIn1, 5, arrMeshVerts))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Models entry."));
		if (!ReadDirEnt(fileIn1, 13, arrModel))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Entities entry."));
		if (!ReadDirEnt(fileIn1, 14, arrEntity))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reading Terr entry."));
		if (!ReadDirEnt(fileIn1, 22, arrTerr))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		fileIn1.Close();
		CStringList listEntityLines;
		if (g_bVerbose) _putts(_T("Parsing Entities entry into text lines."));
		if (!ParseChunk2EntityLines(arrEntity, listEntityLines))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (g_bVerbose) _putts(_T("Reparsing Entities text lines into map file format."));
		fileInto1.SetLength(0);
		if (!ReparseTextLines2Map(fileInto1, rOpts, arrTexture, arrFace, arrVertex, arrModel, arrEntity, arrMeshVerts, arrTerr, listEntityLines))
			AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		if (rOpts.nMask & MB2MML_F_MAKE_BMP) {
			if (!MohBSPToWinDib::CreateDibFrom(rOpts, arrFace, arrVertex, arrModel, arrEntity, arrMeshVerts, arrTerr))
				AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
			//if (!CreateMapDib(rOpts, arrTexture, arrFace, arrVertex, arrModel, arrEntity, arrMeshVerts))
			//	AfxThrowArchiveException(CArchiveException::badIndex, fileIn1.GetFilePath());
		}
		return 0;
	CATCH_ALL(e)

	END_CATCH_ALL

	return -1;
}

int _cdecl _MohBSP2Map(LPCTSTR lpszIn1, LPCTSTR lpszInto1)
{
	MB2MMLOptical aOpts;
	return _MohBSP2Map(lpszIn1, lpszInto1, aOpts);
}

#endif // <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
