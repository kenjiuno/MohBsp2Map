// MohBSPToMohMAP.h : MOHBSPTOMOHMAP アプリケーションのメイン ヘッダー ファイルです。
//

#if !defined(AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_)
#define AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// メイン シンボル

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPApp:
// このクラスの動作の定義に関しては MohBSPToMohMAP.cpp ファイルを参照してください。
//

class CMohBSPToMohMAPApp : public CWinApp
{
public:
	CMohBSPToMohMAPApp();

// オーバーライド
	// ClassWizard は仮想関数のオーバーライドを生成します。
	//{{AFX_VIRTUAL(CMohBSPToMohMAPApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// インプリメンテーション

	//{{AFX_MSG(CMohBSPToMohMAPApp)
		// メモ - ClassWizard はこの位置にメンバ関数を追加または削除します。
		//        この位置に生成されるコードを編集しないでください。
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ は前行の直前に追加の宣言を挿入します。

#endif // !defined(AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_)
