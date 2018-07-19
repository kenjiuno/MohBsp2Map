// MohBSPToMohMAPDlg.h : ヘッダー ファイル
//

#if !defined(AFX_MOHBSPTOMOHMAPDLG_H__8FE2724D_72D8_43D1_BE0F_D4C52E65AF91__INCLUDED_)
#define AFX_MOHBSPTOMOHMAPDLG_H__8FE2724D_72D8_43D1_BE0F_D4C52E65AF91__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPDlg ダイアログ

class CMohBSPToMohMAPDlg : public CDialog
{
public:
	// 
	CMohBSPToMohMAPDlg(CWnd* pParent = NULL);

protected:
	// 
	HICON m_hIcon;
	// 
	CFont m_font;

public:
	//{{AFX_DATA(CMohBSPToMohMAPDlg)
	enum { IDD = IDD_MOHBSPTOMOHMAP_DIALOG };
	CEdit	m_wndEdit;
	CString	m_strFileIn;
	CString	m_strFileOut;
	CString	m_strMes;
	int		m_iReduction;
	UINT	m_nThick;
	BOOL	m_fUniqueTris;
	BOOL	m_fRemoveSharp;
	BOOL	m_fRecoverTexfc;
	BOOL	m_fLowPriority;
	BOOL	m_fSkew;
	//}}AFX_DATA

	//{{AFX_VIRTUAL(CMohBSPToMohMAPDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);
	//}}AFX_VIRTUAL

protected:
	//{{AFX_MSG(CMohBSPToMohMAPDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnOK();
	afx_msg void OnCancel();
	afx_msg void OnButtonRefIn();
	afx_msg void OnButtonRefOut();
	afx_msg void OnButtonStart();
	afx_msg void OnButtonCopy();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

private:

};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ は前行の直前に追加の宣言を挿入します。

#endif // !defined(AFX_MOHBSPTOMOHMAPDLG_H__8FE2724D_72D8_43D1_BE0F_D4C52E65AF91__INCLUDED_)
