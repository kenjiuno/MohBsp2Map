// MohBSPToMohMAPDlg.cpp : インプリメンテーション ファイル
//

#include "stdafx.h"
#include "MohBSPToMohMAP.h"
#include "MohBSPToMohMAPDlg.h"

#include "MohBSP2MAP.h"
#include "OSP.h"
#include "SplitStr.h"
#include "Str2Clipbrd.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace OSP;

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPDlg 

CMohBSPToMohMAPDlg::CMohBSPToMohMAPDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMohBSPToMohMAPDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CMohBSPToMohMAPDlg)
	m_strFileIn = _T("");
	m_strFileOut = _T("");
	m_strMes = _T("");
	m_iReduction = -1;
	m_nThick = 0;
	m_fUniqueTris = FALSE;
	m_fRemoveSharp = FALSE;
	m_fRecoverTexfc = FALSE;
	m_fLowPriority = FALSE;
	m_fSkew = FALSE;
	//}}AFX_DATA_INIT
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

	m_nThick = 4;
	m_iReduction = 2;
	m_fRecoverTexfc = false;
	m_fLowPriority = true;
	m_fSkew = true;
	m_fRemoveSharp = true;

#ifdef _DEBUG

#if 0
#elif 1
	m_strFileIn = "T:\\intr3.bsp";
	m_strFileOut = "T:\\decomp.map";
#elif 1
	m_fRecoverTexfc = true;
	m_strFileIn = "";
	m_strFileOut = "T:\\decomp.map";
#elif 0
	m_strFileIn = "T:\\exp13.bsp";
	m_strFileOut = "T:\\decomp.map";
#elif 0
	m_strFileIn = "T:\\intr8.bsp";
	m_strFileOut = "T:\\decomp.map";
#elif 0
	m_strFileIn = "T:\\intr6.bsp";
	m_strFileOut = "T:\\decomp.map";
#endif

#endif
}

void CMohBSPToMohMAPDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CMohBSPToMohMAPDlg)
	DDX_Control(pDX, IDC_EDIT_LOG, m_wndEdit);
	DDX_Text(pDX, IDC_EDIT_FILE_IN, m_strFileIn);
	DDX_Text(pDX, IDC_EDIT_FILE_OUT, m_strFileOut);
	DDX_Text(pDX, IDC_EDIT_LOG, m_strMes);
	DDX_CBIndex(pDX, IDC_COMBO_REDUCTION, m_iReduction);
	DDX_Text(pDX, IDC_EDIT_THICK, m_nThick);
	DDV_MinMaxUInt(pDX, m_nThick, 1, 256);
	DDX_Check(pDX, IDC_CHECK_UNIQUE_TRIS, m_fUniqueTris);
	DDX_Check(pDX, IDC_CHECK_REMOVE_SHARP, m_fRemoveSharp);
	DDX_Check(pDX, IDC_CHECK_RECOVER_TEXFC, m_fRecoverTexfc);
	DDX_Check(pDX, IDC_CHECK_LOW_PRIORITY, m_fLowPriority);
	DDX_Check(pDX, IDC_CHECK_SKEW_SINK, m_fSkew);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CMohBSPToMohMAPDlg, CDialog)
	//{{AFX_MSG_MAP(CMohBSPToMohMAPDlg)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_REF_IN, OnButtonRefIn)
	ON_BN_CLICKED(IDC_BUTTON_REF_OUT, OnButtonRefOut)
	ON_BN_CLICKED(IDC_BUTTON_START, OnButtonStart)
	ON_BN_CLICKED(IDC_BUTTON_COPY, OnButtonCopy)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPDlg 

BOOL CMohBSPToMohMAPDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	SetIcon(m_hIcon, TRUE);
	SetIcon(m_hIcon, FALSE);

	CFont f;
	f.CreateStockObject(DEFAULT_GUI_FONT);
	LOGFONT lf;
	f.GetLogFont(&lf);
	lf.lfPitchAndFamily = FIXED_PITCH;
	_tcscpy(lf.lfFaceName, _T(""));
	m_font.CreateFontIndirect(&lf);
	if (m_font.m_hObject) m_wndEdit.SetFont(&m_font);
	
	return TRUE;
}

void CMohBSPToMohMAPDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this);

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

HCURSOR CMohBSPToMohMAPDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CMohBSPToMohMAPDlg::OnOK()
{

}

void CMohBSPToMohMAPDlg::OnCancel()
{
	CDialog::OnCancel();
}

const DWORD nDefSaveFlags = 0
	|OFN_ENABLESIZING
	|OFN_EXPLORER
	|OFN_HIDEREADONLY
	|OFN_LONGNAMES
	|OFN_OVERWRITEPROMPT
	|OFN_PATHMUSTEXIST
	;
const DWORD nDefOpenFlags = 0
	|OFN_ENABLESIZING
	|OFN_EXPLORER
	|OFN_FILEMUSTEXIST
	|OFN_HIDEREADONLY
	|OFN_LONGNAMES
	|OFN_PATHMUSTEXIST
	;

void CMohBSPToMohMAPDlg::OnButtonRefIn() 
{
	if (!UpdateData()) return;

	CFileDialog wndDlg(true, _T("bsp"), m_strFileIn, nDefOpenFlags, _T("MOHAA BSP (*.bsp)|*.bsp|All (*.*)|*.*||"));
	if (wndDlg.DoModal() != IDOK) return;
	m_strFileIn = wndDlg.GetPathName();

	UpdateData(false);
}

void CMohBSPToMohMAPDlg::OnButtonRefOut() 
{
	if (!UpdateData()) return;

	CFileDialog wndDlg(false, _T("map"), m_strFileOut, nDefSaveFlags, _T("MOHAA MAP (*.map)|*.map|All (*.*)|*.*||"));
	if (wndDlg.DoModal() != IDOK) return;
	m_strFileOut = wndDlg.GetPathName();

	UpdateData(false);
}

namespace
{
	// 
	CString MkN2NStat(DataWatcher::Name2N &stat)
	{
		typedef std::multimap<UINT, CString, greater<UINT> > N2Name;

		CString str, strText;

		N2Name n2n;
		DataWatcher::Name2N::iterator
			iter1Pos = stat.begin(),
			iter1End = stat.end();
		for (; iter1Pos != iter1End; iter1Pos++) n2n.insert(N2Name::value_type(iter1Pos->second, iter1Pos->first));
		N2Name::iterator
			iter2Pos = n2n.begin(),
			iter2End = n2n.end();
		for (; iter2Pos != iter2End; iter2Pos++) {
			if (iter2Pos->second.IsEmpty())
				continue;
			str.Format(
				"%s\t%3u\n"
				, (LPCTSTR)iter2Pos->second
				, iter2Pos->first
				);
			strText += str;
		}
		return strText;
	}
	// 
	CString FormatLocalDateTime(CTime timeThen)
	{
		tm osTime;
		CTime t(mktime(timeThen.GetLocalTm(&osTime)));
		return t.Format(_T("%Y-%m-%d %H:%M:%S"));
	}
	// 
	CString FormatCostTime(CTimeSpan t)
	{
		CString str;
		if (1 <= t.GetDays()) return "XX:XX:XX";

		str.Format(
			"%02u:%02u:%02u"
			, t.GetHours()
			, t.GetMinutes()
			, t.GetSeconds()
			);
		return str;
	}
	// 
	class ProcessPriorityEnforcer
	{
		// 
		DWORD nSaved;

	public:
		// 
		ProcessPriorityEnforcer()
		{
			nSaved = 0;
		}
		// 
		~ProcessPriorityEnforcer()
		{
			Close();
		}
		// 
		void Close()
		{
			if (nSaved != 0) {
				VERIFY(SetPriorityClass(GetCurrentProcess(), nSaved));
				nSaved = 0;
			}
		}
		// 
		void Enforce(DWORD nNew)
		{
			DWORD n = GetPriorityClass(GetCurrentProcess());
			if (n == 0) return;
			Close();
			nSaved = n;

			VERIFY(SetPriorityClass(GetCurrentProcess(), nNew));
		}

	};
};

void CMohBSPToMohMAPDlg::OnButtonStart() 
{
	if (!UpdateData()) return;

	CWaitCursor wc;

	MB2MM::Knowledge kb;
	MB2MM::Optz optz;
	optz.nThick = m_nThick;
	optz.nReduction = m_iReduction;
	optz.fUniqueTris = m_fUniqueTris ? true : false;
	optz.fRemoveSharp = m_fRemoveSharp ? true : false;
	optz.fRecoverTexfc = m_fRecoverTexfc ? true : false;
	optz.fSkew = m_fSkew ? true : false;
	DataWatcher mo;
	MB2MM::Decompiler z(kb, mo, optz);
	CString strDir = OSP_GetDir(OSP_GetModuleFileName(NULL));
	kb.LoadBBoxFiles(strDir);
	kb.LoadTexMetrFiles(strDir);
	ProcessPriorityEnforcer ef;
	if (m_fLowPriority) ef.Enforce(IDLE_PRIORITY_CLASS);
	mo.timeStart = CTime::GetCurrentTime();
	bool fOk = z.Decompile(m_strFileIn, m_strFileOut);
	mo.timeEnd = CTime::GetCurrentTime();

	CString strText;
	CString str;
	strText.Format(
		"Result ... %s\n"
		"\n"
		"Working Time\n"
		"From %s\n"
		"To   %s\n"
		"\n"
		"Approximate time spent: %s\n"
		"\n"
		"=== World Entity Stat. ===\n"
		"Tris           %7u\n"
		"Bad Tris       %7u\n"
		"Brushes        %7u\n"
		"Junk Tris      %7u\n"
		"Patch Meshes   %7u\n"
		"\n"
		"=== Non-world Entities Stat. ===\n"
		"Tris           %7u\n"
		"Bad Tris       %7u\n"
		"Brushes        %7u\n"
		"Junk Tris      %7u\n"
		"Patch Meshes   %7u\n"
		"\n"
		"=== Other Stat. ===\n"
		"LOD Terrains   %7u\n"
		//
		, fOk ? _T("Ok") : _T("ERR")
		, (LPCTSTR)FormatLocalDateTime(mo.timeStart)
		, (LPCTSTR)FormatLocalDateTime(mo.timeEnd)
		, (LPCTSTR)FormatCostTime(mo.timeEnd - mo.timeStart)
		//
		, mo.nWorldTris
		, mo.nWorldBadTris
		, mo.nWorldBrushes
		, mo.nWorldJunkTris
		, mo.nWorldPMesh
		//
		, mo.nNwTris
		, mo.nNwBadTris
		, mo.nNwBrushes
		, mo.nNwJunkTris
		, mo.nNwPMesh
		//
		, mo.nLODt
		);
	strText += "\n" "=== Entity Usage Stat. ===\n";
	strText += MkN2NStat(mo.entityStat);
	strText += "\n" "=== Model Name Usage Stat. ===\n";
	strText += MkN2NStat(mo.modelStat);

	strText.Replace("\n", "\r\n");

	m_wndEdit.SetWindowText(strText);

	AfxMessageBox(fOk
		? "Ok"
		: "ERR"
		, 0 |(fOk ? MB_ICONINFORMATION : MB_ICONEXCLAMATION)
		);
	if (!fOk) return;
}

void CMohBSPToMohMAPDlg::OnButtonCopy() 
{
	if (!UpdateData()) return;

	CStr2Clipbrd().SetClipboardTextData2(m_strMes);
}
