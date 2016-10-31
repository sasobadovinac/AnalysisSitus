//-----------------------------------------------------------------------------
// Created on: 14 December 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

#ifndef gui_viewer_tess_h
#define gui_viewer_tess_h

// A-Situs includes
#include <analysis_situs.h>

// Common includes
#include <asiData_Model.h>

// GUI includes
#include <gui_viewer.h>

// Visualization includes
#include <visu_interactor_style_pick.h>
#include <visu_pick_callback.h>
#include <visu_rotation_callback.h>

// VTK includes
#include <vtkOrientationMarkerWidget.h>

//! Viewer for mesh.
class gui_viewer_mesh : public gui_viewer
{
  Q_OBJECT

public:

  gui_viewer_mesh(QWidget* parent = NULL);
  virtual ~gui_viewer_mesh();

public:

  void Repaint();

public slots:

  void onResetView();
  void onMeshNodePicked();

private:

  //! Default interactor style.
  vtkSmartPointer<visu_interactor_style_pick> m_interactorStyleDefault;

  //! Pick callback.
  vtkSmartPointer<visu_pick_callback> m_pickCallback;

  //! Rotation callback.
  vtkSmartPointer<visu_rotation_callback> m_rotoCallback;

  //! Orientation marker widget.
  vtkSmartPointer<vtkOrientationMarkerWidget> m_axesWidget;

};

#endif