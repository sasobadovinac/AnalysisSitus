//-----------------------------------------------------------------------------
// Created on: 26 November 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <visu_mesh_prs.h>

// Visualization includes
#include <visu_display_mode.h>
#include <visu_mesh_contour_pipeline.h>
#include <visu_mesh_data_provider.h>
#include <visu_mesh_pipeline.h>
#include <visu_utils.h>

// GUI includes
#include <gui_common.h>

// VTK includes
#include <vtkActor.h>
#include <vtkMapper.h>
#include <vtkProperty.h>

//! Creates a Presentation object for the passed Mesh Node.
//! \param theNode [in] Mesh Node to create a Presentation for.
visu_mesh_prs::visu_mesh_prs(const Handle(ActAPI_INode)& theNode) : visu_prs(theNode)
{
  // Create Data Provider
  Handle(visu_mesh_data_provider)
    DP = new visu_mesh_data_provider( theNode->GetId(),
                                      ActAPI_ParameterStream() << theNode->Parameter(asiData_TessNode::PID_Mesh) );

  // Pipeline for mesh
  this->addPipeline(Pipeline_Mesh, new visu_mesh_pipeline);
  this->assignDataProvider(Pipeline_Mesh, DP);

  // Pipeline for mesh contour
  this->addPipeline(Pipeline_MeshContour, new visu_mesh_contour_pipeline);
  this->assignDataProvider(Pipeline_MeshContour, DP);
  //
  this->GetPipeline(Pipeline_MeshContour)->Actor()->GetProperty()->SetOpacity(0.5);
  this->GetPipeline(Pipeline_MeshContour)->Actor()->SetPickable(0);

  // We use CONTOUR mesh pipeline along with an ordinary one. Thus it is
  // really necessary to resolve coincident primitives to avoid blinking
  // on mesh edges
  vtkMapper::SetResolveCoincidentTopology(1);

  /* =====================================
   *  Prepare a pipeline for highlighting
   * ===================================== */

  Handle(visu_mesh_contour_pipeline) aHiliPipeline = new visu_mesh_contour_pipeline();

  // Set color, opacity and line width
  double aHiliColor[3];
  visu_utils::DefaultDetectionColor(aHiliColor[0], aHiliColor[1], aHiliColor[2]);
  aHiliPipeline->Actor()->GetProperty()->SetColor(aHiliColor[0], aHiliColor[1], aHiliColor[2]);
  aHiliPipeline->Actor()->GetProperty()->SetOpacity(0.75);
  aHiliPipeline->Actor()->GetProperty()->SetLineWidth( visu_utils::DefaultDetectionLineWidth() );

  // Set picking pipeline
  this->installDetectPipeline( aHiliPipeline, DP->Clone() );
}

//! Factory method for Node's Presentation.
//! \param theNode [in] Mesh Node to create a Presentation for.
//! \return new Presentation instance.
Handle(visu_prs) visu_mesh_prs::Instance(const Handle(ActAPI_INode)& theNode)
{
  return new visu_mesh_prs(theNode);
}

//! Returns true if the Presentation is visible, false -- otherwise.
//! \return true/false.
bool visu_mesh_prs::IsVisible() const
{
  return true; // TODO: make visibility controllable
}

//! Sets SHADING visualization mode.
void visu_mesh_prs::doShading() const
{
  Handle(visu_mesh_pipeline) aMeshPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_Mesh) );
  Handle(visu_mesh_pipeline) aMeshContourPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_MeshContour) );

  aMeshPL->Actor()->GetProperty()->EdgeVisibilityOff();
  aMeshPL->Actor()->GetProperty()->SetRepresentationToSurface();
  aMeshContourPL->Actor()->VisibilityOn();
}

//! Sets WIREFRAME visualization mode.
void visu_mesh_prs::doWireframe() const
{
  Handle(visu_mesh_pipeline) aMeshPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_Mesh) );
  Handle(visu_mesh_pipeline) aMeshContourPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_MeshContour) );

  aMeshPL->Actor()->GetProperty()->EdgeVisibilityOn();
  aMeshPL->Actor()->GetProperty()->SetRepresentationToWireframe();
  aMeshContourPL->Actor()->VisibilityOff();
}

//! Sets custom color for the Mesh.
//! \param theColor [in] color to set.
void visu_mesh_prs::doColor(const QColor& theColor) const
{
  if ( !theColor.isValid() )
    return;

  Handle(visu_mesh_pipeline) aMeshPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_Mesh) );

  aMeshPL->Mapper()->ScalarVisibilityOff();
  aMeshPL->Actor()->GetProperty()->SetColor( theColor.redF(),
                                             theColor.greenF(),
                                             theColor.blueF() );
}

//! Unsets custom color for the Mesh.
void visu_mesh_prs::doUnColor() const
{
  Handle(visu_mesh_pipeline) aMeshPL =
    Handle(visu_mesh_pipeline)::DownCast( this->GetPipeline(Pipeline_Mesh) );

  aMeshPL->Mapper()->ScalarVisibilityOn();
}

//! Callback for initialization of Presentation pipelines.
void visu_mesh_prs::beforeInitPipelines()
{
  // Do nothing...
}

//! Callback for initialization of Presentation pipelines.
void visu_mesh_prs::afterInitPipelines()
{
  // Access a dedicated pipeline for highlighting
  const Handle(visu_mesh_contour_pipeline)& detect_pl =
    Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );

  // Initialize the pipeline's input
  detect_pl->SetInput( this->dataProvider(Pipeline_Mesh) );
}

//! Callback for updating of Presentation pipelines invoked before the
//! kernel update routine starts.
void visu_mesh_prs::beforeUpdatePipelines() const
{
  Handle(asiData_TessNode) Mesh_Node = Handle(asiData_TessNode)::DownCast( this->GetNode() );

  visu_display_mode aDMode = (visu_display_mode) Mesh_Node->GetDisplayMode();
  if ( aDMode == DisplayMode_Undefined || aDMode == DisplayMode_Shading )
    this->doShading();
  else if ( aDMode == DisplayMode_Wireframe )
    this->doWireframe();
}

//! Callback for updating of Presentation pipelines invoked after the
//! kernel update routine completes.
void visu_mesh_prs::afterUpdatePipelines() const
{
  /* ======================================
   *  Update highlighting pipeline as well
   * ====================================== */

  // Access a dedicated pipeline for highlighting
  const Handle(visu_mesh_contour_pipeline)& detect_pl =
    Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );

  // IMPORTANT: We update our highlighting pipeline here just to make things
  // faster. The better place to do that is "highlight" method, because
  // we do not really need to build highlighting pipeline just after
  // the Nodal Presentation is created. Logically, we would better to prepare
  // this pipeline only on actual pick request from user. However, in the
  // latter case the reactivity of application might significantly slow down
  detect_pl->Update();

  /* =================
   *  Actualize color
   * ================= */

  Handle(asiData_TessNode) Mesh_Node = Handle(asiData_TessNode)::DownCast( this->GetNode() );
  if ( Mesh_Node->HasColor() )
  {
    QColor aColor = gui_common::IntToColor( Mesh_Node->GetColor() );
    this->doColor(aColor);
  }
  else
    this->doUnColor();
}

//! Callback for highlighting.
//! \param theRenderer  [in] renderer.
//! \param thePickRes   [in] picking results.
//! \param theSelNature [in] selection kind.
void visu_mesh_prs::highlight(vtkRenderer*                 ASitus_NotUsed(theRenderer),
                              const visu_pick_result&      ASitus_NotUsed(thePickRes),
                              const visu_selection_nature& theSelNature) const
{
  //---------------------------------------------------------------------------
  // Update highlighting pipelines
  //---------------------------------------------------------------------------

  // Access pipeline for highlighting
  Handle(visu_mesh_contour_pipeline) hili_pl;
  //
  if ( theSelNature == SelectionNature_Pick )
    hili_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetPickPipeline() );
  else
    hili_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );

  if ( !hili_pl )
    return;

  // Set selection mask...
  //hili_pl->SetSelectedNode(cellIds);
  //hili_pl->ForceExecution();
  hili_pl->SetInput( this->dataProviderDetect() );

  // ... and visibility
  //hili_pl->Actor()->SetVisibility( !cellIds.IsEmpty() );
}

//! Callback for un-highlighting.
//! \param theRenderer  [in] renderer.
//! \param theSelNature [in] selection kind.
void visu_mesh_prs::unHighlight(vtkRenderer*                 ASitus_NotUsed(theRenderer),
                                const visu_selection_nature& theSelNature) const
{
  // Access pipeline for highlighting
  Handle(visu_mesh_contour_pipeline) hili_pl;
  //
  if ( theSelNature == SelectionNature_Pick )
    hili_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetPickPipeline() );
  else
    hili_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );

  if ( !hili_pl )
    return;

  // ... and visibility
  hili_pl->Actor()->SetVisibility(0);
}

//! Callback for rendering.
//! \param theRenderer [in] renderer.
void visu_mesh_prs::renderPipelines(vtkRenderer* theRenderer) const
{
  //---------------------------------------------------------------------------
  // Highlighting
  //---------------------------------------------------------------------------

  const Handle(visu_mesh_contour_pipeline)&
    detect_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );

  // Picking pipeline must be added to renderer the LAST (!). Otherwise
  // we experience some strange coloring bug because of their coincidence
  /* (1) */ detect_pl->AddToRenderer(theRenderer);
}

//! Callback for de-rendering.
//! \param theRenderer [in] renderer.
void visu_mesh_prs::deRenderPipelines(vtkRenderer* theRenderer) const
{
  //---------------------------------------------------------------------------
  // Highlighting
  //---------------------------------------------------------------------------

  Handle(visu_mesh_contour_pipeline)
    detect_pl = Handle(visu_mesh_contour_pipeline)::DownCast( this->GetDetectPipeline() );
  //
  detect_pl->RemoveFromRenderer(theRenderer);
}