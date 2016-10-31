//-----------------------------------------------------------------------------
// Created on: 06 April 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <visu_points_source.h>

// Visualization includes
#include <visu_utils.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

//-----------------------------------------------------------------------------
// Construction
//-----------------------------------------------------------------------------

vtkStandardNewMacro(visu_points_source);

//! Default constructor.
visu_points_source::visu_points_source() : vtkPolyDataAlgorithm ()
{
  this->SetNumberOfInputPorts(0); // Connected directly to our own Data Provider
                                  // which has nothing to do with VTK pipeline
}

//! Destructor.
visu_points_source::~visu_points_source()
{
}

//-----------------------------------------------------------------------------
// Kernel methods
//-----------------------------------------------------------------------------

//! Sets input points to visualize.
//! \param points [in] points to visualize.
void visu_points_source::SetInputPoints(const Handle(asiAlgo_PointCloud)& points)
{
  m_points = points;
  //
  this->Modified();
}

//-----------------------------------------------------------------------------

//! This is called by the superclass. Creates VTK polygonal data set
//! from the input arrays.
//! \param request      [in] information about data object.
//! \param inputVector  [in] the input data. As a data source is the start
//!                          stage of the VTK pipeline, the Input Vector is
//!                          empty and not used (no input port).
//! \param outputVector [in] the pointer to output data, that is filled
//!                          in this method.
//! \return status.
int visu_points_source::RequestData(vtkInformation*        request,
                                    vtkInformationVector** inputVector,
                                    vtkInformationVector*  outputVector)
{
  if ( m_points.IsNull() )
  {
    vtkErrorMacro( << "Invalid domain: NULL point cloud" );
    return 0;
  }

  Handle(TColStd_HArray1OfReal) coords = m_points->GetPoints();
  //
  if ( coords.IsNull() )
  {
    vtkErrorMacro( << "Invalid domain: NULL point cloud" );
    return 0;
  }

  /* ==============================
   *  Prepare involved collections
   * ============================== */

  vtkPolyData* polyOutput = vtkPolyData::GetData(outputVector);
  polyOutput->Allocate();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  polyOutput->SetPoints(points);

  //---------------------------------------------------------------------------

  for ( int i = coords->Lower(); i <= coords->Upper() - 2; i += 3 )
  {
    gp_Pnt P( coords->Value(i), coords->Value(i + 1), coords->Value(i + 2) );
    //
    this->registerVertex( this->registerGridPoint(P, polyOutput), polyOutput );
  }

  return Superclass::RequestData(request, inputVector, outputVector);
}

//! Adds the given point to the passed polygonal data set.
//! \param point    [in]     point to add.
//! \param polyData [in/out] polygonal data set being populated.
//! \return ID of the just added VTK point.
vtkIdType visu_points_source::registerGridPoint(const gp_Pnt& point,
                                                vtkPolyData*  polyData)
{
  // Access necessary arrays
  vtkPoints* points = polyData->GetPoints();

  // Push the point into VTK data set
  vtkIdType pid = points->InsertNextPoint( point.X(),
                                           point.Y(),
                                           point.Z() );

  return pid;
}

//! Adds a vertex cell into the polygonal data set.
//! \param n        [in]     index of the corresponding geometric point.
//! \param polyData [in/out] polygonal data set being populated.
//! \return ID of the just added VTK cell.
vtkIdType visu_points_source::registerVertex(const vtkIdType n,
                                             vtkPolyData*    polyData)
{
  std::vector<vtkIdType> nodes;
  nodes.push_back(n);
  //
  vtkIdType cellID = polyData->InsertNextCell(VTK_VERTEX, 1, &nodes[0]);
  //
  return cellID;
}