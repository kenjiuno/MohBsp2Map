// MohBSPToMohMAP.h : MOHBSPTOMOHMAP �A�v���P�[�V�����̃��C�� �w�b�_�[ �t�@�C���ł��B
//

#if !defined(AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_)
#define AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// ���C�� �V���{��

/////////////////////////////////////////////////////////////////////////////
// CMohBSPToMohMAPApp:
// ���̃N���X�̓���̒�`�Ɋւ��Ă� MohBSPToMohMAP.cpp �t�@�C�����Q�Ƃ��Ă��������B
//

class CMohBSPToMohMAPApp : public CWinApp
{
public:
	CMohBSPToMohMAPApp();

// �I�[�o�[���C�h
	// ClassWizard �͉��z�֐��̃I�[�o�[���C�h�𐶐����܂��B
	//{{AFX_VIRTUAL(CMohBSPToMohMAPApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// �C���v�������e�[�V����

	//{{AFX_MSG(CMohBSPToMohMAPApp)
		// ���� - ClassWizard �͂��̈ʒu�Ƀ����o�֐���ǉ��܂��͍폜���܂��B
		//        ���̈ʒu�ɐ��������R�[�h��ҏW���Ȃ��ł��������B
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ �͑O�s�̒��O�ɒǉ��̐錾��}�����܂��B

#endif // !defined(AFX_MOHBSPTOMOHMAP_H__17DCF008_0B64_4704_8604_299BD79D6948__INCLUDED_)
