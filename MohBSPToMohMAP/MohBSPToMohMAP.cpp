// MohBSPToMohMAP.cpp : アプリケーション用クラスの定義を行います。
//

#include "stdafx.h"
#include "MohBSPToMohMAP.h"
#include "MohBSPToMohMAPDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPApp

BEGIN_MESSAGE_MAP(CMohBSPToMohMAPApp, CWinApp)
	//{{AFX_MSG_MAP(CMohBSPToMohMAPApp)
		// メモ - ClassWizard はこの位置にマッピング用のマクロを追加または削除します。
		//        この位置に生成されるコードを編集しないでください。
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPApp クラスの構築

CMohBSPToMohMAPApp::CMohBSPToMohMAPApp(): CWinApp(_T("MohBSPToMohMAP 0.0.2"))
{
	// TODO: この位置に構築用のコードを追加してください。
	// ここに InitInstance 中の重要な初期化処理をすべて記述してください。
}

/////////////////////////////////////////////////////////////////////////////
// 唯一の CMohBSPToMohMAPApp オブジェクト

CMohBSPToMohMAPApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPApp クラスの初期化

BOOL CMohBSPToMohMAPApp::InitInstance()
{
	// 標準的な初期化処理
	// もしこれらの機能を使用せず、実行ファイルのサイズを小さくしたけ
	//  れば以下の特定の初期化ルーチンの中から不必要なものを削除して
	//  ください。

	SetRegistryKey(_T("Kentaro-K.21"));

	CMohBSPToMohMAPDlg dlg;
	m_pMainWnd = &dlg;
	int nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{
		// TODO: ダイアログが <OK> で消された時のコードを
		//       記述してください。
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: ダイアログが <ｷｬﾝｾﾙ> で消された時のコードを
		//       記述してください。
	}

	// ダイアログが閉じられてからアプリケーションのメッセージ ポンプを開始するよりは、
	// アプリケーションを終了するために FALSE を返してください。
	return FALSE;
}
