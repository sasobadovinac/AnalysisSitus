//-----------------------------------------------------------------------------
// Created on: 02 December 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

#ifndef visu_face_domain_pipeline_h
#define visu_face_domain_pipeline_h

// A-Situs includes
#include <visu_data_provider.h>
#include <visu_pipeline.h>

// VTK includes
#include <vtkExtractSelection.h>
#include <vtkGeometryFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>

//-----------------------------------------------------------------------------

DEFINE_STANDARD_HANDLE(visu_face_domain_pipeline, visu_pipeline)

//! Visualization pipeline for face domain.
class visu_face_domain_pipeline : public visu_pipeline
{
public:

  // OCCT RTTI
  DEFINE_STANDARD_RTTI_INLINE(visu_face_domain_pipeline, visu_pipeline)

public:

  visu_face_domain_pipeline(const bool isDefaultColorScheme = true);

public:

  virtual void
    SetInput(const Handle(visu_data_provider)& DP);

public:

  inline void ForceExecution() { m_bForced = true; }

public:

  void
    SetSelectedCells(const TColStd_PackedMapOfInteger& mask);

private:

  virtual void callback_add_to_renderer      (vtkRenderer* theRenderer);
  virtual void callback_remove_from_renderer (vtkRenderer* theRenderer);
  virtual void callback_update               ();

private:

  //! Copying prohibited.
  visu_face_domain_pipeline(const visu_face_domain_pipeline&);

  //! Assignment prohibited.
  visu_face_domain_pipeline& operator=(const visu_face_domain_pipeline&);

private:

  double computeTipSize(const TopoDS_Face& F) const;

private:

  bool                                 m_bDefaultColorScheme; //!< Indicates whether to use a default color scheme.
  bool                                 m_bForced;             //!< Forced update.
  bool                                 m_bMapperColorsSet;    //!< Boolean flag indicating whether lookup table is set.
  vtkSmartPointer<vtkIdTypeArray>      m_selected;            //!< Poles selected for visualization.
  vtkSmartPointer<vtkSelectionNode>    m_selectionNode;       //!< VTK selection node.
  vtkSmartPointer<vtkSelection>        m_selection;           //!< VTK selection.
  vtkSmartPointer<vtkExtractSelection> m_extractSelection;    //!< VTK selection extractor.
  vtkSmartPointer<vtkGeometryFilter>   m_toPolyData;          //!< VTK geometry filter.

};

#endif