//-----------------------------------------------------------------------------
// Created on: 25 September 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/, http://quaoar.su/
//-----------------------------------------------------------------------------

#ifndef analysis_situs_h
#define analysis_situs_h

#define ASitus_NotUsed(x)

#ifdef ASitus_EXPORTS
  #define ASitus_EXPORT __declspec(dllexport)
#else
  #define ASitus_EXPORT __declspec(dllimport)
#endif

// Active Data (API) includes
#include <ActAPI_Common.h>

// VTK includes
#include <vtkAutoInit.h>

#endif
