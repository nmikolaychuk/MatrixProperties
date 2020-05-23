
// MatrixPropertiesDlg.cpp: файл реализации
//

#include "pch.h"
#include "framework.h"
#include "MatrixProperties.h"
#include "MatrixPropertiesDlg.h"
#include "afxdialogex.h"

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <random>
#include <cassert>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define DOTS(x,y) (xp*((x)-xmin)),(yp*((y)-ymax))							// макрос перевода координат. график сигнала

using namespace std;

// Диалоговое окно CMatrixPropertiesDlg

struct PRNG
{
	std::mt19937 engine;
};

void initGenerator(PRNG& generator)
{
	// Создаём псевдо-устройство для получения случайного зерна.
	std::random_device device;
	// Получаем случайное зерно последовательности
	generator.engine.seed(device());
}

CMatrixPropertiesDlg::CMatrixPropertiesDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_MATRIXPROPERTIES_DIALOG, pParent)
	, dim(7)
	, InformationToScreen(_T(""))
	, minimum(-50)
	, maximum(50)
	, cond(1000)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMatrixPropertiesDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_DIM, dim);
	DDX_Text(pDX, IDC_INFORMATION, InformationToScreen);
	DDX_Text(pDX, IDC_KOEF_MIN, minimum);
	DDX_Text(pDX, IDC_KOEF_MAX, maximum);
	DDX_Text(pDX, IDC_COND, cond);
	DDX_Control(pDX, IDC_RADIO_CHOOSE_COND, m_choose_cond);
	DDX_Control(pDX, IDC_RADIO_RAND_COND, m_rand_cond);
	DDX_Control(pDX, IDC_COND, status);
	DDX_Control(pDX, IDC_RADIO_OBR_GRAPH, m_obr_graph);
	DDX_Control(pDX, IDC_RADIO_PSEVDO_GRAPH, m_psevdo_graph);
}

BEGIN_MESSAGE_MAP(CMatrixPropertiesDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_EXIT, &CMatrixPropertiesDlg::OnBnClickedExit)
	ON_BN_CLICKED(IDC_SVD, &CMatrixPropertiesDlg::OnBnClickedSvd)
	ON_BN_CLICKED(IDC_RADIO_RAND_COND, &CMatrixPropertiesDlg::OnBnClickedRadioRandCond)
	ON_BN_CLICKED(IDC_RADIO_CHOOSE_COND, &CMatrixPropertiesDlg::OnBnClickedRadioChooseCond)
	ON_BN_CLICKED(IDC_DRAW_GRAPH, &CMatrixPropertiesDlg::OnBnClickedDrawGraph)
END_MESSAGE_MAP()


// Обработчики сообщений CMatrixPropertiesDlg

BOOL CMatrixPropertiesDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Задает значок для этого диалогового окна.  Среда делает это автоматически,
	//  если главное окно приложения не является диалоговым
	SetIcon(m_hIcon, TRUE);			// Крупный значок
	SetIcon(m_hIcon, FALSE);		// Мелкий значок

	// TODO: добавьте дополнительную инициализацию

	PicWnd = GetDlgItem(IDC_GRAPH);			//связываем с ID окон
	PicDc = PicWnd->GetDC();
	PicWnd->GetClientRect(&Pic);

	// перья
	setka_pen.CreatePen(		//для сетки
		PS_DOT,					//пунктирная
		-1,						//толщина 1 пиксель
		RGB(0, 0, 0));			//цвет  черный

	osi_pen.CreatePen(			//координатные оси
		PS_SOLID,				//сплошная линия
		2,						//толщина 2 пикселя
		RGB(0, 0, 0));			//цвет черный

	obratnaya_pen.CreatePen(			//график
		PS_SOLID,				//сплошная линия
		-1,						//толщина -1 пикселя
		RGB(255, 0, 0));			//цвет синий

	psevdo_pen.CreatePen(			//график
		PS_SOLID,				//сплошная линия
		-1,						//толщина -1 пикселя
		RGB(0, 0, 255));			//цвет синий

	UpdateData(FALSE);
	return TRUE;  // возврат значения TRUE, если фокус не передан элементу управления
}

// При добавлении кнопки свертывания в диалоговое окно нужно воспользоваться приведенным ниже кодом,
//  чтобы нарисовать значок.  Для приложений MFC, использующих модель документов или представлений,
//  это автоматически выполняется рабочей областью.

void CMatrixPropertiesDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // контекст устройства для рисования

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Выравнивание значка по центру клиентского прямоугольника
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Нарисуйте значок
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
		RedrawAll(ymin, ymax);
	}
}

// Система вызывает эту функцию для получения отображения курсора при перемещении
//  свернутого окна.
HCURSOR CMatrixPropertiesDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CMatrixPropertiesDlg::Mashtab(float arr[], int dim, float* mmin, float* mmax)		//определяем функцию масштабирования
{
	*mmin = *mmax = arr[0];

	for (int i = 0; i < dim; i++)
	{
		if (*mmin > arr[i]) *mmin = arr[i];
		if (*mmax < arr[i]) *mmax = arr[i];
	}
}

void CMatrixPropertiesDlg::RedrawAll(float ymin, float ymax)
{
	PicDc->FillSolidRect(&Pic, RGB(250, 250, 250));			//закрашиваю фон 
	PicDc->SelectObject(&osi_pen);		//выбираем перо

	float window_signal_width = Pic.Width();
	float window_signal_height = Pic.Height();
	xp = (window_signal_width / (xmax - xmin));			//Коэффициенты пересчёта координат по Х
	yp = -(window_signal_height / (ymax - ymin));			//Коэффициенты пересчёта координат по У

	//создаём Ось Y
	PicDc->MoveTo(DOTS(0, ymax));
	PicDc->LineTo(DOTS(0, ymin));
	//создаём Ось Х
	PicDc->MoveTo(DOTS(xmin, 0));
	PicDc->LineTo(DOTS(xmax, 0));

	PicDc->SelectObject(&setka_pen);

	//отрисовка сетки по y
	for (float x = 0; x <= xmax; x += max_cond / 10)
	{
		if (x != 0) {
			PicDc->MoveTo(DOTS(x, ymax));
			PicDc->LineTo(DOTS(x, ymin));
		}
	}
	//отрисовка сетки по x
	for (float y = 0; y < ymax; y += ymax / 5)
	{
		if (y != 0) {
			PicDc->MoveTo(DOTS(xmin, y));
			PicDc->LineTo(DOTS(xmax, y));
		}
	}


	//подпись точек на оси
	CFont font;
	font.CreateFontW(14.5, 0, 0, 0, FW_REGULAR, 0, 0, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS || CLIP_LH_ANGLES, DEFAULT_QUALITY, DEFAULT_PITCH, _T("Century Gothic"));
	PicDc->SelectObject(font);

	//подпись осей
	PicDc->TextOutW(DOTS(0.02 * max_cond, 0.98 * ymax), _T("p")); //Y
	PicDc->TextOutW(DOTS(xmax - max_cond / 15, -ymin - 0.1 * ymax), _T("Cond")); //X

	//по Y с шагом 5
	for (float i = 0; i <= ymax; i += ymax / 5)
	{
		CString str;
		if (i != 0)
		{
			str.Format(_T("%.3f"), i);
			PicDc->TextOutW(DOTS(xmin + max_cond / 50, i + 0.03 * ymax), str);
		}
	}
	//по X с шагом 0.5
	for (float j = 0; j <= xmax; j += max_cond / 10)
	{
		CString str;
		if (j != 0) {
			str.Format(_T("%.e"), j);
			PicDc->TextOutW(DOTS(j - max_cond / 100, -ymin - 0.22 * ymax), str);
		}
	}
}

int CMatrixPropertiesDlg::svd_hestenes(int m_m, int n_n, float* a, float* u, float* v, float* sigma)
{
	float thr = 1.E-4f, nul = 1.E-16f;
	int n, m, i, j, l, k, lort, iter, in, ll, kk;
	float alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
	n = n_n;
	m = m_m;
	for (i = 0; i < n; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
			if (i == j) v[in + j] = 1.;
			else v[in + j] = 0.;
	}
	for (i = 0; i < m; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
		{
			u[in + j] = a[in + j];
		}
	}

	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l < n - 1; l++)
			for (k = l + 1; k < n; k++)
			{
				alfa = 0.; betta = 0.; hamma = 0.;
				for (i = 0; i < m; i++)
				{
					in = i * n;
					ll = in + l;
					kk = in + k;
					alfa += u[ll] * u[ll];
					betta += u[kk] * u[kk];
					hamma += u[ll] * u[kk];
				}

				if (sqrt(alfa * betta) < nul)	continue;
				if (fabs(hamma) / sqrt(alfa * betta) < thr) continue;

				lort = 1;
				eta = (betta - alfa) / (2.f * hamma);
				t = (float)((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta * eta)));
				cos0 = (float)(1. / sqrt(1. + t * t));
				sin0 = t * cos0;

				for (i = 0; i < m; i++)
				{
					in = i * n;
					buf = u[in + l] * cos0 - u[in + k] * sin0;
					u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
					u[in + l] = buf;

					if (i >= n) continue;
					buf = v[in + l] * cos0 - v[in + k] * sin0;
					v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
					v[in + l] = buf;
				}
			}

		if (!lort) break;
	}

	for (i = 0; i < n; i++)
	{
		s = 0.;
		for (j = 0; j < m; j++)	s += u[j * n + i] * u[j * n + i];
		s = (float)sqrt(s);
		sigma[i] = s;
		if (s < nul)	continue;
		for (j = 0; j < m; j++)	u[j * n + i] /= s;
	}
	//======= Sortirovka ==============
	for (i = 0; i < n - 1; i++)
		for (j = i; j < n; j++)
			if (sigma[i] < sigma[j])
			{
				s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
				for (k = 0; k < m; k++)
				{
					s = u[i + k * n]; u[i + k * n] = u[j + k * n]; u[j + k * n] = s;
				}
				for (k = 0; k < n; k++)
				{
					s = v[i + k * n]; v[i + k * n] = v[j + k * n]; v[j + k * n] = s;
				}
			}

	return iter;
}

// Генерирует число с плавающей точкой в диапазоне [minValue, maxValue)
float getRandomFloat(PRNG& generator, float minValue, float maxValue)
{
	// Проверяем корректность аргументов
	assert(minValue < maxValue);

	// Создаём распределение
	std::uniform_real_distribution<float> distribution(minValue, maxValue);

	// Вычисляем псевдослучайное число: вызовем распределение как функцию,
	//  передав генератор произвольных целых чисел как аргумент.
	return distribution(generator.engine);
}

// ----------------------------- Определитель --------------------------
// Вычисляет определитель матрицы T размерностью N
float CMatrixPropertiesDlg::det(float** T, UINT32 N)
{
	float det__;
	int sub_j, s;
	float** subT;    // Субматрица как набор ссылок на исходную матрицу
	switch (N)
	{
	case 1:
		return T[0][0];
	case 2:
		return T[0][0] * T[1][1] - T[0][1] * T[1][0];
	default:
		if (N < 1)  // Некорректная матрица
		{
			MessageBox(L"Некорректная матрица!", L"Ошибка", MB_ICONERROR | MB_OK);
		}
		subT = new float* [N - 1];  // Массив ссылок на столбцы субматрицы
		det__ = 0;
		s = 1;        // Знак минора
		for (int i = 0; i < N; i++)  // Разложение по первому столбцу
		{
			sub_j = 0;
			for (int j = 0; j < N; j++)// Заполнение субматрицы ссылками на исходные столбцы
				if (i != j)      // исключить i строку
					subT[sub_j++] = T[j] + 1;  // здесь + 1 исключает первый столбец

			det__ = det__ + s * T[i][0] * det(subT, N - 1);
			s = -s;
		};
		delete[] subT;
		return det__;
	};
}

// Функция, производящая обращение матрицы.
// Принимает:
//     matrix - матрица для обращения
//     result - матрица достаточного размера для вмещения результата
//     size   - размерность матрицы
// Возвращает:
//     true в случае успешного обращения, false в противном случае
bool CMatrixPropertiesDlg::inverse(float** matrix, float** result, int size)
{
	// Изначально результирующая матрица является единичной
	// Заполняем единичную матрицу
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
			result[i][j] = 0.0;

		result[i][i] = 1.0;
	}

	// Копия исходной матрицы
	float** copy = new float* [size]();

	// Заполняем копию исходной матрицы
	for (int i = 0; i < size; ++i)
	{
		copy[i] = new float[size];

		for (int j = 0; j < size; ++j)
			copy[i][j] = matrix[i][j];
	}

	// Проходим по строкам матрицы (назовём их исходными)
	// сверху вниз. На данном этапе происходит прямой ход
	// и исходная матрица превращается в верхнюю треугольную
	for (int k = 0; k < size; ++k)
	{
		// Если элемент на главной диагонали в исходной
		// строке - нуль, то ищем строку, где элемент
		// того же столбца не нулевой, и меняем строки
		// местами
		if (fabs(copy[k][k]) < 1e-8)
		{
			// Ключ, говорязий о том, что был произведён обмен строк
			bool changed = false;

			// Идём по строкам, расположенным ниже исходной
			for (int i = k + 1; i < size; ++i)
			{
				// Если нашли строку, где в том же столбце
				// имеется ненулевой элемент
				if (fabs(copy[i][k]) > 1e-8)
				{
					// Меняем найденную и исходную строки местами
					// как в исходной матрице, так и в единичной
					std::swap(copy[k], copy[i]);
					std::swap(result[k], result[i]);

					// Взводим ключ - сообщаем о произведённом обмене строк
					changed = true;

					break;
				}
			}

			// Если обмен строк произведён не был - матрица не может быть
			// обращена
			if (!changed)
			{
				// Чистим память
				for (int i = 0; i < size; ++i)
					delete[] copy[i];

				delete[] copy;

				// Сообщаем о неудаче обращения
				return false;
			}
		}

		// Запоминаем делитель - диагональный элемент
		float div = copy[k][k];

		// Все элементы исходной строки делим на диагональный
		// элемент как в исходной матрице, так и в единичной
		for (int j = 0; j < size; ++j)
		{
			copy[k][j] /= div;
			result[k][j] /= div;
		}

		// Идём по строкам, которые расположены ниже исходной
		for (int i = k + 1; i < size; ++i)
		{
			// Запоминаем множитель - элемент очередной строки,
			// расположенный под диагональным элементом исходной
			// строки
			float multi = copy[i][k];

			// Отнимаем от очередной строки исходную, умноженную
			// на сохранённый ранее множитель как в исходной,
			// так и в единичной матрице
			for (int j = 0; j < size; ++j)
			{
				copy[i][j] -= multi * copy[k][j];
				result[i][j] -= multi * result[k][j];
			}
		}
	}

	// Проходим по вернхней треугольной матрице, полученной
	// на прямом ходе, снизу вверх
	// На данном этапе происходит обратный ход, и из исходной
	// матрицы окончательно формируется единичная, а из единичной -
	// обратная
	for (int k = size - 1; k > 0; --k)
	{
		// Идём по строкам, которые расположены выше исходной
		for (int i = k - 1; i + 1 > 0; --i)
		{
			// Запоминаем множитель - элемент очередной строки,
			// расположенный над диагональным элементом исходной
			// строки
			float multi = copy[i][k];

			// Отнимаем от очередной строки исходную, умноженную
			// на сохранённый ранее множитель как в исходной,
			// так и в единичной матрице
			for (int j = 0; j < size; ++j)
			{
				copy[i][j] -= multi * copy[k][j];
				result[i][j] -= multi * result[k][j];
			}
		}
	}

	// Чистим память
	for (int i = 0; i < size; ++i)
		delete[] copy[i];

	delete[] copy;

	// Сообщаем об успехе обращения
	return true;
}

void CMatrixPropertiesDlg::OnBnClickedExit()
{
	// TODO: добавьте свой код обработчика уведомлений
	CDialogEx::OnCancel();
}


void CMatrixPropertiesDlg::OnBnClickedSvd()
{
	// TODO: добавьте свой код обработчика уведомлений
	UpdateData(TRUE);

	PRNG generator;
	initGenerator(generator);

	if (m_rand_cond.GetCheck() == BST_CHECKED)
	{
		float** A = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			A[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				A[i][j] = getRandomFloat(generator, minimum, maximum);
			}
		}

		float** Aobr = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Aobr[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Aobr[i][j] = 0;
			}
		}

		ofstream out("result_rand_cond.txt");

		out << "\t\t\t\t\t\t      УЧЕБНАЯ ПРАКТИКА\n";
		out << "\t\t\t     Задание 12. Исследование свойств обратной и псевдообратной матриц.";

		out << "\n\nМатрица коэффициентов А:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << A[i][j] << "\t";
			}
			out << "\n";
		}

		float determinant = 0;
		determinant = det(A, dim);
		out << "\nОпределитель матрицы: " << determinant;

		out << "\n\n******************************************************** ОБРАТНАЯ МАТРИЦА (КЛАССИЧЕСКИЙ МЕТОД) *******************************************************";

		if (determinant != 0)
		{
			inverse(A, Aobr, dim);

			out << "\n\nОбратная матрица к матрице А:\n\n";
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					out << setprecision(3) << Aobr[i][j] << "\t";
				}
				out << "\n";
			}
		}

		//решение методом псевдообратной матрицы Мура-Пенроуза

		float* MatrA = new float[dim * dim];
		float* MatrU = new float[dim * dim];
		float* MatrV = new float[dim * dim];
		float* VectorSigm = new float[dim];
		float* MatrSigm = new float[dim * dim];
		float* TranspMatrSigm = new float[dim * dim];

		float* Buffer = new float[dim * dim];
		float* Buffer1 = new float[dim];
		float* Buffer2 = new float[dim * dim];
		float* Buffer3 = new float[dim * dim];

		for (int i = 0; i < dim; i++)
		{
			VectorSigm[i] = 0;
			Buffer1[i] = 0;
			for (int j = 0; j < dim; j++)
			{
				MatrA[i * dim + j] = A[i][j];
				MatrU[i * dim + j] = 0;
				MatrV[i * dim + j] = 0;
				MatrSigm[i * dim + j] = 0;
				TranspMatrSigm[i * dim + j] = 0;
				Buffer[i * dim + j] = 0;
				Buffer2[i * dim + j] = 0;
				Buffer3[i * dim + j] = 0;
			}
		}

		out << "\n************************************************************** ПСЕВДООБРАТНАЯ МАТРИЦА (СВД) **************************************************************";

		svd_hestenes(dim, dim, MatrA, MatrU, MatrV, VectorSigm);

		float porog = VectorSigm[0] / 10;

		//транспонируем матрицу U
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				Buffer[i * dim + j] = MatrU[i + j * dim];
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
					MatrSigm[i * dim + j] = VectorSigm[i];
				else
					MatrSigm[i * dim + j] = 0;
			}
		}

		//транспонируем матрицу Sigma
		for (int i = 0; i < dim; i++)
		{
			if (VectorSigm[i] <= porog)
				Buffer1[i] = 0;
			else
				Buffer1[i] = 1 / VectorSigm[i];
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
					TranspMatrSigm[i * dim + j] = Buffer1[i];
				else
					TranspMatrSigm[i * dim + j] = 0;
			}
		}

		//умножаем V на транспонированную Sigma
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					Buffer2[i * dim + j] += MatrV[i * dim + k] * TranspMatrSigm[k * dim + j];
				}
			}
		}

		//умножаем результат предыдущего действия на транспонированную матрицу U (определение псевдообратной матрицы A-крест)
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					Buffer3[i * dim + j] += Buffer2[i * dim + k] * Buffer[k * dim + j];
				}
			}
		}

		out << "\n\nВектор Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			out << setprecision(3) << VectorSigm[i] << "\t";
		}

		out << "\n\nМатрица Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << MatrSigm[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\nТранспонированная матрица Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << TranspMatrSigm[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\n\nМатрица U:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << MatrU[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\nМатрица V:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << MatrV[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\nПроизведение VSigma^-1:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Buffer2[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\nПсевдообратная матрица Мура-Пенроуза (A+):\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Buffer3[i * dim + j] << "\t";
			}
			out << "\n";
		}

		//число обусловленности

		float cond_rand = VectorSigm[0] / VectorSigm[dim - 1];

		out << "\nЧисло обусловленности матрицы А: " << cond_rand;

		//проверка свойств

		//AA^-1 = E      и      A^-1A = E
		float** edinichnaya_1 = new float* [dim];
		float** edinichnaya_2 = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			edinichnaya_1[i] = new float[dim];
			edinichnaya_2[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				edinichnaya_1[i][j] = 0;
				edinichnaya_2[i][j] = 0;
				for (int k = 0; k < dim; k++)
				{
					edinichnaya_1[i][j] += A[i][k] * Aobr[k][j];
					edinichnaya_2[i][j] += Aobr[i][k] * A[k][j];
				}
			}
		}

		/*for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (edinichnaya_1[i][j] < 0.001)
				{
					edinichnaya_1[i][j] = 0;
				}
				if (edinichnaya_2[i][j] < 0.001)
				{
					edinichnaya_2[i][j] = 0;
				}
			}
		}*/

		out << "\n\n******************************************************************** ПРОВЕРКА СВОЙСТВ ************************************************************************";

		out << "\n\nСвойства обратной матрицы:";

		out << "\n\n1 свойство: AA^-1 = E и A^-1A = E:\n\n";
		out << "AA^-1:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << edinichnaya_1[i][j] << "\t";
			}
			out << "\n";
		}

		out << "A^-1A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << edinichnaya_2[i][j] << "\t";
			}
			out << "\n";
		}
		out << "\t\t\t\t\t\t\t\t\t1 свойство: выполняется!\n";

		float determ_obr = 0;
		determ_obr = det(Aobr, dim);

		out << "\n2 свойство: det(A^-1) = 1/det(A)\n\n";
		
		out << "det(A^-1) = " << determ_obr << " = 1/det(A)\n";

		out << "\t\t\t\t\t\t\t\t\t2 свойство: выполняется!";

		out << "\n\n3 свойство: (A^-1)T = (AT)^-1\n\n";
		out << "(A^-1)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Aobr[j][i] << "\t";
			}
			out << "\n";
		}

		out << "(AT)^-1\n";
		float** Atransp = new float* [dim];
		float** Atransp_obr = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Atransp[i] = new float[dim];
			Atransp_obr[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Atransp[i][j] = A[j][i];
				Atransp_obr[i][j] = 0;
			}
		}
		inverse(Atransp, Atransp_obr, dim);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Atransp_obr[i][j] << "\t";
			}
			out << "\n";
		}

		out << "\t\t\t\t\t\t\t\t\t3 свойство: выполняется!\n";

		out << "\n4 свойство: (A^-1)^-1 = A\n\n";

		float** Aobr_obr = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Aobr_obr[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Aobr_obr[i][j] = 0;
			}
		}
		out << "(A^-1)^-1:\n";
		inverse(Aobr, Aobr_obr, dim);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Aobr_obr[i][j] << "\t";
			}
			out << "\n";
		}

		out << "\t\t\t\t\t\t\t\t\t4 свойство: выполняется!\n";

		out << "\nСвойства псевдообратной матрицы Мура-Пенроуза:";

		//AA+A = A, A+AA+ = A+

		float** svoistvo1 = new float* [dim];
		float** svoistvo2 = new float* [dim];
		float** svoistvo3_1 = new float* [dim];
		float** svoistvo3_2 = new float* [dim];

		float** svoistvo1_res = new float* [dim];
		float** svoistvo2_res = new float* [dim];
		float** svoistvo3_1_res = new float* [dim];
		float** svoistvo3_2_res = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			svoistvo1[i] = new float[dim];
			svoistvo2[i] = new float[dim];
			svoistvo3_1[i] = new float[dim];
			svoistvo3_2[i] = new float[dim];

			svoistvo1_res[i] = new float[dim];
			svoistvo2_res[i] = new float[dim];
			svoistvo3_1_res[i] = new float[dim];
			svoistvo3_2_res[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				svoistvo1[i][j] = 0;
				svoistvo2[i][j] = 0;
				svoistvo3_1[i][j] = 0;
				svoistvo3_2[i][j] = 0;

				svoistvo1_res[i][j] = 0;
				svoistvo2_res[i][j] = 0;
				svoistvo3_1_res[i][j] = 0;
				svoistvo3_2_res[i][j] = 0;
				for (int k = 0; k < dim; k++)
				{
					svoistvo1[i][j] += A[i][k] * Buffer3[k * dim + j];
					svoistvo2[i][j] += Buffer3[i * dim + k] * A[k][j];
					svoistvo3_1[i][j] += A[i][k] * Buffer3[k * dim + j];
					svoistvo3_2[i][j] += Buffer3[i * dim + k] * A[k][j];
				}
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				svoistvo3_1_res[i][j] = svoistvo3_1[j][i];
				svoistvo3_2_res[i][j] = svoistvo3_2[j][i];
				for (int k = 0; k < dim; k++)
				{
					svoistvo1_res[i][j] += svoistvo1[i][k] * A[k][j];
					svoistvo2_res[i][j] += svoistvo2[i][k] * Buffer3[k * dim + j];
				}
			}
		}

		out << "\n\n1 свойство: AA+A=A\n\n";
		out << "AA+A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo1_res[i][j] << "\t";
			}
			out << "\n";
		}

		out << "A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << A[i][j] << "\t";
			}
			out << "\n";
		}

		out << "\t\t\t\t\t\t\t\t\t1 свойство: выполняется!\n";

		out << "\n2 свойство: A+AA+=A+\n\n";
		out << "A+AA+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo2_res[i][j] << "\t";
			}
			out << "\n";
		}

		out << "A+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << Buffer3[i * dim + j] << "\t";
			}
			out << "\n";
		}

		out << "\t\t\t\t\t\t\t\t\t2 свойство: выполняется!\n";

		out << "\n3 свойство: (AA+)T = AA+ и (A+A)T = A+A\n\n";
		out << "(AA+)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo3_1_res[i][j] << "\t";
			}
			out << "\n";
		}

		out << "AA+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo3_1[i][j] << "\t";
			}
			out << "\n";
		}

		out << "(A+A)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo3_2_res[i][j] << "\t";
			}
			out << "\n";
		}

		out << "A+A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << setprecision(3) << svoistvo3_2[i][j] << "\t";
			}
			out << "\n";
		}

		out << "\t\t\t\t\t\t\t\t\t3 свойство: выполняется!\n";

		out.close();

		CFile file;
		file.Open(L"result_rand_cond.txt", CFile::modeRead);
		CStringA str;
		LPSTR pBuf = str.GetBuffer(file.GetLength() + 1);
		file.Read(pBuf, file.GetLength() + 1);
		pBuf[file.GetLength()] = NULL;
		CStringA decodedText = str;
		file.Close();
		str.ReleaseBuffer();

		InformationToScreen = "";
		InformationToScreen = str;

		UpdateData(FALSE);


		//открытие файла
		//ShellExecute(NULL, L"open", L"result.txt", NULL, NULL, SW_SHOWNORMAL);

		//освобождение памяти
		for (int i = 0; i < dim; ++i)
		{
			delete[] A[i];
		}
		delete[] A;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Aobr[i];
		}
		delete[] Aobr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Aobr_obr[i];
		}
		delete[] Aobr_obr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Atransp[i];
		}
		delete[] Atransp;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Atransp_obr[i];
		}
		delete[] Atransp_obr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] edinichnaya_1[i];
			delete[] edinichnaya_2[i];
		}
		delete[] edinichnaya_1;
		delete[] edinichnaya_2;

		for (int i = 0; i < dim; ++i)
		{
			delete[] svoistvo1[i];
			delete[] svoistvo2[i];
			delete[] svoistvo3_1[i];
			delete[] svoistvo3_2[i];

			delete[] svoistvo1_res[i];
			delete[] svoistvo2_res[i];
			delete[] svoistvo3_1_res[i];
			delete[] svoistvo3_2_res[i];
		}
		delete[] svoistvo1;
		delete[] svoistvo2;
		delete[] svoistvo3_1;
		delete[] svoistvo3_2;

		delete[] svoistvo1_res;
		delete[] svoistvo2_res;
		delete[] svoistvo3_1_res;
		delete[] svoistvo3_2_res;

		delete[] MatrA;
		delete[] MatrU;
		delete[] MatrSigm;
		delete[] MatrV;
		delete[] Buffer;
		delete[] Buffer1;
		delete[] Buffer2;
		delete[] Buffer3;
	}

	if (m_choose_cond.GetCheck() == BST_CHECKED)
	{
		float** A = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			A[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				A[i][j] = getRandomFloat(generator, minimum, maximum);
			}
		}

		ofstream out("result_choose_cond.txt");

		out << "\t\t\t\t\t\t      УЧЕБНАЯ ПРАКТИКА\n";
		out << "\t\t\t     Задание 12. Исследование свойств обратной и псевдообратной матриц.";

		float* MatrA = new float[dim * dim];
		float* MatrU = new float[dim * dim];
		float* MatrV = new float[dim * dim];

		float* FinalMatrU = new float[dim * dim];
		float* FinalMatrV = new float[dim * dim];
		float* FinalVectorSigm = new float[dim];
		float* FinalMatrSigm = new float[dim * dim];

		float* VectorSigm = new float[dim];
		float* NewVectorSigm = new float[dim];
		float* NewMatrSigm = new float[dim * dim];
		float* Inner = new float[dim - 2];
		float* NewMatrA = new float[dim * dim];
		float* NewMatrA_res = new float[dim * dim];
		float* TranspMatrSigm = new float[dim * dim];

		float* Buffer = new float[dim * dim];
		float* Buffer1 = new float[dim];
		float* Buffer2 = new float[dim * dim];
		float* Buffer3 = new float[dim * dim];

		for (int i = 0; i < dim - 2; i++)
		{
			Inner[i] = 0;
		}

		for (int i = 0; i < dim; i++)
		{
			VectorSigm[i] = 0;
			NewVectorSigm[i] = 0;
			FinalVectorSigm[i] = 0;
			Buffer1[i] = 0;
			for (int j = 0; j < dim; j++)
			{
				MatrA[i * dim + j] = A[i][j];
				MatrU[i * dim + j] = 0;
				MatrV[i * dim + j] = 0;

				FinalMatrU[i * dim + j] = 0;
				FinalMatrV[i * dim + j] = 0;
				FinalMatrSigm[i * dim + j] = 0;

				NewMatrSigm[i * dim + j] = 0;
				TranspMatrSigm[i * dim + j] = 0;
				NewMatrA[i * dim + j] = 0;
				NewMatrA_res[i * dim + j] = 0;
				Buffer[i * dim + j] = 0;
				Buffer2[i * dim + j] = 0;
				Buffer3[i * dim + j] = 0;
			}
		}

		svd_hestenes(dim, dim, MatrA, MatrU, MatrV, VectorSigm);

		for (int i = 0; i < dim; i++)
		{
			NewVectorSigm[i] = VectorSigm[i];
		}

		NewVectorSigm[dim - 1] = NewVectorSigm[0] / cond;

		for (int i = 0; i < dim - 2; i++)
		{
			Inner[i] = getRandomFloat(generator, NewVectorSigm[dim - 1], NewVectorSigm[0]);
		}

		for (int j = 0; j < dim - 3; j++)
		{
			int max_element = j;
			for (int i = j + 1; i < dim - 2; i++)
			{
				if (Inner[i] > Inner[max_element])
				{
					max_element = i;
				}
			}
			swap(Inner[j], Inner[max_element]);
		}

		for (int i = 1; i < dim - 1; i++)
		{
			NewVectorSigm[i] = Inner[i - 1];
		}

		//собираем обратно матрицу A

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
					NewMatrSigm[i * dim + j] = NewVectorSigm[i];
				else
					NewMatrSigm[i * dim + j] = 0;
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					NewMatrA[i * dim + j] += MatrU[i * dim + k] * NewMatrSigm[k * dim + j];
				}
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					NewMatrA_res[i * dim + j] += NewMatrA[i * dim + k] * MatrV[k * dim + j];
				}
			}
		}

		out << "\n\nМатрица коэффициентов А:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << NewMatrA_res[i * dim + j] << " ";
			}
			out << "\n";
		}

		float** Aobr = new float* [dim];
		float** ANew = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Aobr[i] = new float[dim];
			ANew[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Aobr[i][j] = 0;
				ANew[i][j] = 0;
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				ANew[i][j] = NewMatrA_res[i * dim + j];
			}
		}

		float determinant = det(ANew, dim);
		out << "\nОпределитель матрицы: " << determinant;

		out << "\n\n******************************************************** ОБРАТНАЯ МАТРИЦА (КЛАССИЧЕСКИЙ МЕТОД) *******************************************************";

		if (determinant != 0)
		{
			inverse(ANew, Aobr, dim);

			out << "\n\nОбратная матрица к матрице А:\n\n";
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					out << Aobr[i][j] << " ";
				}
				out << "\n";
			}
		}

		out << "\n************************************************************** ПСЕВДООБРАТНАЯ МАТРИЦА (СВД) **************************************************************";

		svd_hestenes(dim, dim, NewMatrA_res, FinalMatrU, FinalMatrV, FinalVectorSigm);

		float porog = FinalVectorSigm[0] / 10;

		//транспонируем матрицу U
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				Buffer[i * dim + j] = FinalMatrU[i + j * dim];
			}
		}

		//транспонируем вектор Sigma
		for (int i = 0; i < dim; i++)
		{
			if (FinalVectorSigm[i] <= porog)
				Buffer1[i] = 0;
			else
				Buffer1[i] = 1 / FinalVectorSigm[i];
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
					TranspMatrSigm[i * dim + j] = Buffer1[i];
				else
					TranspMatrSigm[i * dim + j] = 0;
			}
		}

		//умножаем V на транспонированную Sigma
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					Buffer2[i * dim + j] += FinalMatrV[i * dim + k] * TranspMatrSigm[k * dim + j];
				}
			}
		}

		//умножаем результат предыдущего действия на транспонированную матрицу U (определение псевдообратной матрицы A-крест)
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					Buffer3[i * dim + j] += Buffer2[i * dim + k] * Buffer[k * dim + j];
				}
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
					FinalMatrSigm[i * dim + j] = FinalVectorSigm[i];
				else
					FinalMatrSigm[i * dim + j] = 0;
			}
		}

		out << "\n\nВектор Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			out << FinalVectorSigm[i] << " ";
		}

		out << "\n\nМатрица Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << FinalMatrSigm[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\nТранспонированная матрица Сигма:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << TranspMatrSigm[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\n\nМатрица U:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << FinalMatrU[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\nМатрица V:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << FinalMatrV[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\nПроизведение VSigma^-1:\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Buffer2[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\nПсевдообратная матрица Мура-Пенроуза (A+):\n\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Buffer3[i * dim + j] << " ";
			}
			out << "\n";
		}

		//число обусловленности

		float cond_rand = FinalVectorSigm[0] / FinalVectorSigm[dim - 1];

		out << "\nЧисло обусловленности матрицы А: " << cond_rand;

		//проверка свойств

		//AA^-1 = E      и      A^-1A = E
		float** edinichnaya_1 = new float* [dim];
		float** edinichnaya_2 = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			edinichnaya_1[i] = new float[dim];
			edinichnaya_2[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				edinichnaya_1[i][j] = 0;
				edinichnaya_2[i][j] = 0;
				for (int k = 0; k < dim; k++)
				{
					edinichnaya_1[i][j] += ANew[i][k] * Aobr[k][j];
					edinichnaya_2[i][j] += Aobr[i][k] * ANew[k][j];
				}
			}
		}

		//for (int i = 0; i < dim; i++)
		//{
		//	for (int j = 0; j < dim; j++)
		//	{
		//		if (edinichnaya_1[i][j] < 0.001)
		//		{
		//			edinichnaya_1[i][j] = 0;
		//		}
		//		if (edinichnaya_2[i][j] < 0.001)
		//		{
		//			edinichnaya_2[i][j] = 0;
		//		}
		//	}
		//}

		out << "\n\n******************************************************************** ПРОВЕРКА СВОЙСТВ ************************************************************************";

		out << "\n\nСвойства обратной матрицы:";

		out << "\n\n1 свойство: AA^-1 = E и A^-1A = E:\n\n";
		out << "AA^-1:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << edinichnaya_1[i][j] << " ";
			}
			out << "\n";
		}

		out << "A^-1A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << edinichnaya_2[i][j] << " ";
			}
			out << "\n";
		}

		out << "\n2 свойство: det(A^-1) = 1/det(A)\n\n";
		out << "det(A^-1) = " << det(Aobr, dim) << " = 1/det(A)";

		out << "\n\n3 свойство: (A^-1)T = (AT)^-1\n\n";
		out << "(A^-1)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Aobr[j][i] << " ";
			}
			out << "\n";
		}

		out << "(AT)^-1\n";
		float** Atransp = new float* [dim];
		float** Atransp_obr = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Atransp[i] = new float[dim];
			Atransp_obr[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Atransp[i][j] = ANew[j][i];
				Atransp_obr[i][j] = 0;
			}
		}
		inverse(Atransp, Atransp_obr, dim);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Atransp_obr[i][j] << " ";
			}
			out << "\n";
		}

		out << "\n4 свойство: (A^-1)^-1 = A\n\n";

		float** Aobr_obr = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			Aobr_obr[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				Aobr_obr[i][j] = 0;
			}
		}
		out << "(A^-1)^-1:\n";
		inverse(Aobr, Aobr_obr, dim);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Aobr_obr[i][j] << " ";
			}
			out << "\n";
		}

		out << "\nСвойства псевдообратной матрицы Мура-Пенроуза:";

		//AA+A = A, A+AA+ = A+

		float** svoistvo1 = new float* [dim];
		float** svoistvo2 = new float* [dim];
		float** svoistvo3_1 = new float* [dim];
		float** svoistvo3_2 = new float* [dim];

		float** svoistvo1_res = new float* [dim];
		float** svoistvo2_res = new float* [dim];
		float** svoistvo3_1_res = new float* [dim];
		float** svoistvo3_2_res = new float* [dim];
		for (int i = 0; i < dim; i++)
		{
			svoistvo1[i] = new float[dim];
			svoistvo2[i] = new float[dim];
			svoistvo3_1[i] = new float[dim];
			svoistvo3_2[i] = new float[dim];

			svoistvo1_res[i] = new float[dim];
			svoistvo2_res[i] = new float[dim];
			svoistvo3_1_res[i] = new float[dim];
			svoistvo3_2_res[i] = new float[dim];
			for (int j = 0; j < dim; j++)
			{
				svoistvo1[i][j] = 0;
				svoistvo2[i][j] = 0;
				svoistvo3_1[i][j] = 0;
				svoistvo3_2[i][j] = 0;

				svoistvo1_res[i][j] = 0;
				svoistvo2_res[i][j] = 0;
				svoistvo3_1_res[i][j] = 0;
				svoistvo3_2_res[i][j] = 0;
				for (int k = 0; k < dim; k++)
				{
					svoistvo1[i][j] += ANew[i][k] * Buffer3[k * dim + j];
					svoistvo2[i][j] += Buffer3[i * dim + k] * ANew[k][j];
					svoistvo3_1[i][j] += ANew[i][k] * Buffer3[k * dim + j];
					svoistvo3_2[i][j] += Buffer3[i * dim + k] * ANew[k][j];
				}
			}
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				svoistvo3_1_res[i][j] = svoistvo3_1[j][i];
				svoistvo3_2_res[i][j] = svoistvo3_2[j][i];
				for (int k = 0; k < dim; k++)
				{
					svoistvo1_res[i][j] += svoistvo1[i][k] * ANew[k][j];
					svoistvo2_res[i][j] += svoistvo2[i][k] * Buffer3[k * dim + j];
				}
			}
		}

		out << "\n\n1 свойство: AA+A=A\n\n";
		out << "AA+A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo1_res[i][j] << " ";
			}
			out << "\n";
		}

		out << "A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << ANew[i][j] << " ";
			}
			out << "\n";
		}

		out << "\n2 свойство: A+AA+=A+\n\n";
		out << "A+AA+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo2_res[i][j] << " ";
			}
			out << "\n";
		}

		out << "A+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << Buffer3[i * dim + j] << " ";
			}
			out << "\n";
		}

		out << "\n3 свойство: (AA+)T = AA+ и (A+A)T = A+A\n\n";
		out << "(AA+)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo3_1_res[i][j] << " ";
			}
			out << "\n";
		}

		out << "AA+:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo3_1[i][j] << " ";
			}
			out << "\n";
		}

		out << "(A+A)T:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo3_2_res[i][j] << " ";
			}
			out << "\n";
		}

		out << "A+A:\n";
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				out << svoistvo3_2[i][j] << " ";
			}
			out << "\n";
		}

		out.close();

		CFile file;
		file.Open(L"result_choose_cond.txt", CFile::modeRead);
		CStringA str;
		LPSTR pBuf = str.GetBuffer(file.GetLength() + 1);
		file.Read(pBuf, file.GetLength() + 1);
		pBuf[file.GetLength()] = NULL;
		CStringA decodedText = str;
		file.Close();
		str.ReleaseBuffer();

		InformationToScreen = "";
		InformationToScreen = str;

		UpdateData(FALSE);

		//освобождение памяти
		for (int i = 0; i < dim; ++i)
		{
			delete[] A[i];
		}
		delete[] A;

		for (int i = 0; i < dim; ++i)
		{
			delete[] ANew[i];
		}
		delete[] ANew;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Aobr[i];
		}
		delete[] Aobr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Aobr_obr[i];
		}
		delete[] Aobr_obr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Atransp[i];
		}
		delete[] Atransp;

		for (int i = 0; i < dim; ++i)
		{
			delete[] Atransp_obr[i];
		}
		delete[] Atransp_obr;

		for (int i = 0; i < dim; ++i)
		{
			delete[] edinichnaya_1[i];
			delete[] edinichnaya_2[i];
		}
		delete[] edinichnaya_1;
		delete[] edinichnaya_2;

		for (int i = 0; i < dim; ++i)
		{
			delete[] svoistvo1[i];
			delete[] svoistvo2[i];
			delete[] svoistvo3_1[i];
			delete[] svoistvo3_2[i];

			delete[] svoistvo1_res[i];
			delete[] svoistvo2_res[i];
			delete[] svoistvo3_1_res[i];
			delete[] svoistvo3_2_res[i];
		}
		delete[] svoistvo1;
		delete[] svoistvo2;
		delete[] svoistvo3_1;
		delete[] svoistvo3_2;

		delete[] svoistvo1_res;
		delete[] svoistvo2_res;
		delete[] svoistvo3_1_res;
		delete[] svoistvo3_2_res;

		delete[] MatrA;
		delete[] MatrU;
		delete[] MatrV;
		delete[] VectorSigm;
		delete[] NewVectorSigm;
		delete[] NewMatrSigm;
		delete[] Inner;
		delete[] NewMatrA;
		delete[] NewMatrA_res;
		delete[] TranspMatrSigm;
		delete[] Buffer;
		delete[] Buffer1;
		delete[] Buffer2;
		delete[] Buffer3;
	}

	if (m_rand_cond.GetCheck() == BST_UNCHECKED && m_choose_cond.GetCheck() == BST_UNCHECKED)
	{
		MessageBox(L"Сначала укажите тип решения задачи (Случайная или заданная обусловленность матрицы А)", L"Внимание", MB_ICONEXCLAMATION | MB_OK);
	}
}


void CMatrixPropertiesDlg::OnBnClickedRadioRandCond()
{
	// TODO: добавьте свой код обработчика уведомлений
	if (m_rand_cond.GetCheck() == BST_CHECKED)
	{
		status.EnableWindow(FALSE);
	}
}


void CMatrixPropertiesDlg::OnBnClickedRadioChooseCond()
{
	// TODO: добавьте свой код обработчика уведомлений
	if (m_choose_cond.GetCheck() == BST_CHECKED)
	{
		status.EnableWindow(TRUE);
	}
}


void CMatrixPropertiesDlg::OnBnClickedDrawGraph()
{
	// TODO: добавьте свой код обработчика уведомлений
	UpdateData(TRUE);

	PRNG generator;
	initGenerator(generator);

	if (m_obr_graph.GetCheck() == BST_UNCHECKED && m_psevdo_graph.GetCheck() == BST_UNCHECKED)
	{
		MessageBox(L"Необходимо указать для какой матрицы выявить зависимость (обратная или псевдообратная)", L"Внимание", MB_ICONEXCLAMATION | MB_OK);
	}

	float** SumElementObrMatr = new float* [kolvo_iter];
	float** SumElementPsevdoMatr = new float* [kolvo_iter];

	for (int i = 0; i < kolvo_iter; i++)
	{
		SumElementObrMatr[i] = new float[max_cond / step];
		SumElementPsevdoMatr[i] = new float[max_cond / step];
		for (int j = 0; j < max_cond / step; j++)
		{
			SumElementObrMatr[i][j] = 0;
			SumElementPsevdoMatr[i][j] = 0;
		}
	}

	for (int kolvo = 0; kolvo < kolvo_iter; kolvo++)
	{
		//основной цикл для построения зависимости

		for (int iter_cond = 500; iter_cond < max_cond; iter_cond += step)
		{
			float** A = new float* [dim];
			for (int i = 0; i < dim; i++)
			{
				A[i] = new float[dim];
				for (int j = 0; j < dim; j++)
				{
					A[i][j] = getRandomFloat(generator, minimum, maximum);
				}
			}

			float* MatrA = new float[dim * dim];
			float* MatrU = new float[dim * dim];
			float* MatrV = new float[dim * dim];

			float* FinalMatrU = new float[dim * dim];
			float* FinalMatrV = new float[dim * dim];
			float* FinalVectorSigm = new float[dim];
			float* FinalMatrSigm = new float[dim * dim];

			float* VectorSigm = new float[dim];
			float* NewVectorSigm = new float[dim];
			float* NewMatrSigm = new float[dim * dim];
			float* Inner = new float[dim - 2];
			float* NewMatrA = new float[dim * dim];
			float* NewMatrA_res = new float[dim * dim];
			float* TranspMatrSigm = new float[dim * dim];

			float* Buffer = new float[dim * dim];
			float* Buffer1 = new float[dim];
			float* Buffer2 = new float[dim * dim];
			float* Buffer3 = new float[dim * dim];

			for (int i = 0; i < dim - 2; i++)
			{
				Inner[i] = 0;
			}

			for (int i = 0; i < dim; i++)
			{
				VectorSigm[i] = 0;
				NewVectorSigm[i] = 0;
				FinalVectorSigm[i] = 0;
				Buffer1[i] = 0;
				for (int j = 0; j < dim; j++)
				{
					MatrA[i * dim + j] = A[i][j];
					MatrU[i * dim + j] = 0;
					MatrV[i * dim + j] = 0;

					FinalMatrU[i * dim + j] = 0;
					FinalMatrV[i * dim + j] = 0;
					FinalMatrSigm[i * dim + j] = 0;

					NewMatrSigm[i * dim + j] = 0;
					TranspMatrSigm[i * dim + j] = 0;
					NewMatrA[i * dim + j] = 0;
					NewMatrA_res[i * dim + j] = 0;
					Buffer[i * dim + j] = 0;
					Buffer2[i * dim + j] = 0;
					Buffer3[i * dim + j] = 0;
				}
			}

			svd_hestenes(dim, dim, MatrA, MatrU, MatrV, VectorSigm);

			for (int i = 0; i < dim; i++)
			{
				NewVectorSigm[i] = VectorSigm[i];
			}

			NewVectorSigm[dim - 1] = NewVectorSigm[0] / iter_cond;

			for (int i = 0; i < dim - 2; i++)
			{
				Inner[i] = getRandomFloat(generator, NewVectorSigm[dim - 1], NewVectorSigm[0]);
			}

			for (int j = 0; j < dim - 3; j++)
			{
				int max_element = j;
				for (int i = j + 1; i < dim - 2; i++)
				{
					if (Inner[i] > Inner[max_element])
					{
						max_element = i;
					}
				}
				swap(Inner[j], Inner[max_element]);
			}

			for (int i = 1; i < dim - 1; i++)
			{
				NewVectorSigm[i] = Inner[i - 1];
			}

			//собираем обратно матрицу A

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					if (i == j)
						NewMatrSigm[i * dim + j] = NewVectorSigm[i];
					else
						NewMatrSigm[i * dim + j] = 0;
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					for (int k = 0; k < dim; k++)
					{
						NewMatrA[i * dim + j] += MatrU[i * dim + k] * NewMatrSigm[k * dim + j];
					}
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					for (int k = 0; k < dim; k++)
					{
						NewMatrA_res[i * dim + j] += NewMatrA[i * dim + k] * MatrV[k * dim + j];
					}
				}
			}

			float** Aobr = new float* [dim];
			float** ANew = new float* [dim];
			for (int i = 0; i < dim; i++)
			{
				Aobr[i] = new float[dim];
				ANew[i] = new float[dim];
				for (int j = 0; j < dim; j++)
				{
					Aobr[i][j] = 0;
					ANew[i][j] = 0;
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					ANew[i][j] = NewMatrA_res[i * dim + j];
				}
			}

			float determinant = det(ANew, 3);

			if (determinant != 0)
			{
				inverse(ANew, Aobr, dim);
			}

			svd_hestenes(dim, dim, NewMatrA_res, FinalMatrU, FinalMatrV, FinalVectorSigm);

			float porog = FinalVectorSigm[0] / 20;

			//транспонируем матрицу U
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					Buffer[i * dim + j] = FinalMatrU[i + j * dim];
				}
			}

			//транспонируем вектор Sigma
			for (int i = 0; i < dim; i++)
			{
				if (FinalVectorSigm[i] <= porog)
					Buffer1[i] = 0;
				else
					Buffer1[i] = 1 / FinalVectorSigm[i];
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					if (i == j)
						TranspMatrSigm[i * dim + j] = Buffer1[i];
					else
						TranspMatrSigm[i * dim + j] = 0;
				}
			}

			//умножаем V на транспонированную Sigma
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					for (int k = 0; k < dim; k++)
					{
						Buffer2[i * dim + j] += FinalMatrV[i * dim + k] * TranspMatrSigm[k * dim + j];
					}
				}
			}

			//умножаем результат предыдущего действия на транспонированную матрицу U (определение псевдообратной матрицы A-крест)
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					for (int k = 0; k < dim; k++)
					{
						Buffer3[i * dim + j] += Buffer2[i * dim + k] * Buffer[k * dim + j];
					}
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					if (i == j)
						FinalMatrSigm[i * dim + j] = FinalVectorSigm[i];
					else
						FinalMatrSigm[i * dim + j] = 0;
				}
			}

			//проверка свойств


			//AA^-1 = E      и      A^-1A = E
			float** edinichnaya_1 = new float* [dim];
			float** E = new float* [dim];
			for (int i = 0; i < dim; i++)
			{
				edinichnaya_1[i] = new float[dim];
				E[i] = new float[dim];
				for (int j = 0; j < dim; j++)
				{
					edinichnaya_1[i][j] = 0;
					E[i][j] = 0;
					for (int k = 0; k < dim; k++)
					{
						edinichnaya_1[i][j] += ANew[i][k] * Aobr[k][j];
					}
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					if (i == j)
						E[i][j] = 1;
					else
						E[i][j] = 0;
				}
			}



			float sum_e = 0;
			float sum_aobr = 0;

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					sum_e += E[i][j];
					sum_aobr += edinichnaya_1[i][j];
				}
			}

			SumElementObrMatr[kolvo][iter_cond / step] = abs(sum_e - sum_aobr);

			//AA+A = A, A+AA+ = A+

			float** svoistvo1 = new float* [dim];
			float** svoistvo1_res = new float* [dim];

			for (int i = 0; i < dim; i++)
			{
				svoistvo1[i] = new float[dim];
				svoistvo1_res[i] = new float[dim];
				for (int j = 0; j < dim; j++)
				{
					svoistvo1[i][j] = 0;
					svoistvo1_res[i][j] = 0;
					for (int k = 0; k < dim; k++)
					{
						svoistvo1[i][j] += ANew[i][k] * Buffer3[k * dim + j];
					}
				}
			}

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					for (int k = 0; k < dim; k++)
					{
						svoistvo1_res[i][j] += svoistvo1[i][k] * ANew[k][j];
					}
				}
			}

			float sum_anew = 0;
			float sum_proizv = 0;

			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					sum_anew += ANew[i][j];
					sum_proizv += svoistvo1_res[i][j];
				}
			}

			SumElementPsevdoMatr[kolvo][iter_cond / step] = abs(sum_anew - sum_proizv);

			//освобождение памяти

			for (int i = 0; i < dim; ++i)
			{
				delete[] ANew[i];
			}
			delete[] ANew;

			for (int i = 0; i < dim; ++i)
			{
				delete[] Aobr[i];
			}
			delete[] Aobr;

			for (int i = 0; i < dim; ++i)
			{
				delete[] edinichnaya_1[i];
				delete[] E[i];
			}
			delete[] edinichnaya_1;
			delete[] E;

			for (int i = 0; i < dim; ++i)
			{
				delete[] svoistvo1[i];
				delete[] svoistvo1_res[i];
			}
			delete[] svoistvo1;
			delete[] svoistvo1_res;

			for (int i = 0; i < dim; ++i)
			{
				delete[] A[i];
			}
			delete[] A;

			delete[] MatrA;
			delete[] MatrU;
			delete[] MatrV;
			delete[] VectorSigm;
			delete[] NewVectorSigm;
			delete[] NewMatrSigm;
			delete[] Inner;
			delete[] NewMatrA;
			delete[] NewMatrA_res;
			delete[] TranspMatrSigm;
			delete[] Buffer;
			delete[] Buffer1;
			delete[] Buffer2;
			delete[] Buffer3;
		}
	}

	float* average_obr = new float[max_cond / step];
	float* average_psevdo = new float[max_cond / step];
	for (int i = 0; i < max_cond / step; i++)
	{
		average_obr[i] = 0;
		average_psevdo[i] = 0;
	}

	for (int i = 0; i < max_cond / step; i++)
	{
		for (int j = 0; j < kolvo_iter; j++)
		{
			average_obr[i] += SumElementObrMatr[j][i];
			average_psevdo[i] += SumElementPsevdoMatr[j][i];
		}
		average_obr[i] /= kolvo_iter;
		average_psevdo[i] /= kolvo_iter;
	}

	if (m_obr_graph.GetCheck() == BST_CHECKED)
	{
		Mashtab(average_obr, max_cond / step, &ymin, &ymax);
		RedrawAll(-0.17 * ymax, ymax);
		PicDc->SelectObject(&obratnaya_pen);
		PicDc->MoveTo(DOTS(500, average_obr[1]));

		for (int i = 500; i < max_cond; i += step)
		{
			PicDc->LineTo(DOTS(i, average_obr[i / step]));
		}
	}

	if (m_psevdo_graph.GetCheck() == BST_CHECKED)
	{
		Mashtab(average_psevdo, max_cond / step, &ymin, &ymax);
		RedrawAll(-0.15 * ymax, ymax * 0.9);
		PicDc->SelectObject(&psevdo_pen);
		PicDc->MoveTo(DOTS(500, average_psevdo[1]));

		for (int i = 500; i < max_cond; i += step)
		{
			PicDc->LineTo(DOTS(i, average_psevdo[i / step]));
		}
	}

	ofstream out("graph_value.txt");

	out << "Для обратной матрицы:\n\n";
	for (int i = 0; i < max_cond / step; i++)
	{
		out << average_obr[i] << ", ";
	}
	out << "\n\nДля псевдообратной матрицы:\n\n";
	for (int i = 0; i < max_cond / step; i++)
	{
		out << average_psevdo[i] << ", ";
	}
	out.close();

	for (int i = 0; i < kolvo_iter; i++)
	{
		delete[] SumElementObrMatr[i];
		delete[] SumElementPsevdoMatr[i];
	}
	delete[] SumElementObrMatr;
	delete[] SumElementPsevdoMatr;

	delete[] average_obr;
	delete[] average_psevdo;
}
