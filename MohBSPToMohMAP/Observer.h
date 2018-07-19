
#pragma once

#include <map>

struct DataWatcher {
	// 
	typedef std::map<CString, UINT> Name2N;

	// 
	CTime timeStart, timeEnd;
	// 
	UINT nLODt, nEntities;
	// World ...
	DWORD nWorldTris, nWorldBadTris, nWorldBrushes, nWorldJunkTris, nWorldPMesh;
	// Non-world ...
	DWORD nNwTris, nNwBadTris, nNwBrushes, nNwJunkTris, nNwPMesh;
	// 
	Name2N entityStat;
	// 
	Name2N modelStat;

	// 
	bool fProcessingWorld;

	// 
	DataWatcher()
	{
		nLODt =
		nWorldTris = nWorldBadTris = nWorldBrushes = nWorldJunkTris = nWorldPMesh =
		nNwTris = nNwBadTris = nNwBrushes = nNwJunkTris = nNwPMesh = 
		nEntities = 
		0;
	}
	// 
	inline void AddLODt() { nLODt++; }
	// 
	inline void AddTri() { if (fProcessingWorld) nWorldTris++; else nNwTris++; }
	// 
	inline void AddBadTri() { if (fProcessingWorld) nWorldBadTris++; else nNwBadTris++; }
	// 
	inline void AddBrush() { if (fProcessingWorld) nWorldBrushes++; else nNwBrushes++; }
	// 
	inline void AddJunkTri() { if (fProcessingWorld) nWorldJunkTris++; else nNwJunkTris++; }
	// 
	inline void AddPMesh() { if (fProcessingWorld) nWorldPMesh++; else nNwPMesh++; }
	// 
	inline void AddEntityClass(CString str) { Add2Stat(entityStat, str); }
	// 
	inline void AddModelName(CString str) { Add2Stat(modelStat, NormalizeModelName(str)); }

protected:
	// 
	void Add2Stat(Name2N &z, CString strKey)
	{
		nEntities++;
		z[strKey]++;
	}

	// 
	static CString NormalizeModelName(CString str)
	{
		str.Replace("//", "/");
		return str;
	}
};
