
// MatrixProperties.h: главный файл заголовка для приложения PROJECT_NAME
//

#pragma once

#ifndef __AFXWIN_H__
	#error "включить pch.h до включения этого файла в PCH"
#endif

#include "resource.h"		// основные символы


// CMatrixPropertiesApp:
// Сведения о реализации этого класса: MatrixProperties.cpp
//

class CMatrixPropertiesApp : public CWinApp
{
public:
	CMatrixPropertiesApp();

// Переопределение
public:
	virtual BOOL InitInstance();

// Реализация

	DECLARE_MESSAGE_MAP()
};

extern CMatrixPropertiesApp theApp;
