
// MatrixPropertiesDlg.h: файл заголовка
//

#pragma once



// Диалоговое окно CMatrixPropertiesDlg
class CMatrixPropertiesDlg : public CDialogEx
{
// Создание
public:
	CMatrixPropertiesDlg(CWnd* pParent = nullptr);	// стандартный конструктор

	int svd_hestenes(int, int, float*, float*, float *, float*);
	float det(float** T, UINT32 N);
	bool inverse(float** matrix, float** result, int size);
	void MatrMultiply(int n, int m, float** matrix, float* vektor, float* res);
	void RedrawAll(float mn, float mx);
	void Mashtab(float arr[], int dim, float* mmin, float* mmax);

// Данные диалогового окна
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MATRIXPROPERTIES_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// поддержка DDX/DDV


// Реализация
protected:
	HICON m_hIcon;

	// Созданные функции схемы сообщений
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedExit();
	afx_msg void OnBnClickedSvd();
	afx_msg void OnBnClickedRadioRandCond();
	afx_msg void OnBnClickedRadioChooseCond();
	afx_msg void OnBnClickedDrawGraph();

	CWnd* PicWnd;				//области рисования
	CDC* PicDc;
	CRect Pic;

	CPen osi_pen;				//ручки
	CPen setka_pen;
	CPen obratnaya_pen;
	CPen psevdo_pen;

	int max_cond = 1.E+6;
	int step = max_cond / 30;
	static const int kolvo_iter = 15;

	float xp = 0, yp = 0,				//коэфициенты пересчета
		xmin = -max_cond*0.1, xmax = max_cond,				//максисимальное и минимальное значение х 
		ymin = -0.001, ymax = 0.006,				//максисимальное и минимальное значение y
		mn = 0, mx = 0;					//коэффициенты масштабирования

	float Pi = 3.141592653589;

	int dim;
	float minimum;
	float maximum;
	float cond;
	
	CString InformationToScreen;

	CButton m_choose_cond;
	CButton m_rand_cond;
	CButton m_obr_graph;
	CButton m_psevdo_graph;

	CEdit status;
};
