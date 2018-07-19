
#pragma once

namespace MB2MM
{
	// 
	class Texc
	{
	public:
		// 
		UINT cxTex, cyTex;
		// 
		double vv[3][2];
		// 
		double tv[3][2];
		// 
		bool fMirror;

		// 
		double tx, ty, angle, tsx, tsy;

		// 
		bool Decode();
	};
};
