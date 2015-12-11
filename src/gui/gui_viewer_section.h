//-----------------------------------------------------------------------------
// Created on: 09 December 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://quaoar.su/blog/
//-----------------------------------------------------------------------------

#ifndef gui_viewer_section_h
#define gui_viewer_section_h

// A-Situs includes
#include <analysis_situs.h>

// A-Situs (GUI) includes
#include <gui_viewer.h>

//! Viewer for a single skinning section.
class gui_viewer_section : public gui_viewer
{
  Q_OBJECT

public:

  gui_viewer_section(QWidget* parent = NULL);
  virtual ~gui_viewer_section();

public:

  void Repaint();

public slots:

  void onResetView();

};

#endif