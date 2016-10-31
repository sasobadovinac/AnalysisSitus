//-----------------------------------------------------------------------------
// Created on: 16 April 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <visu_iv_point_set_2d_data_provider.h>

// Common includes
#include <common_facilities.h>

//-----------------------------------------------------------------------------

//! Constructor.
//! \param pointSet [in] point set to source data from.
visu_iv_point_set_2d_data_provider::visu_iv_point_set_2d_data_provider(const Handle(ActAPI_INode)& pointSet)
: visu_points_data_provider(pointSet)
{}

//-----------------------------------------------------------------------------

//! \return point cloud to visualize.
Handle(asiAlgo_PointCloud) visu_iv_point_set_2d_data_provider::GetPoints() const
{
  Handle(asiData_IVPointSet2dNode)
    points_n = Handle(asiData_IVPointSet2dNode)::DownCast(m_node);
  //
  if ( points_n.IsNull() || !points_n->IsWellFormed() )
    return NULL;

  Handle(TColStd_HArray1OfReal) coords2d = points_n->GetPoints();
  if ( coords2d.IsNull() )
    return NULL;

  // Re-pack to 3d
  const int nPts = coords2d->Length() / 2;
  Handle(TColStd_HArray1OfReal) coords3d = new TColStd_HArray1OfReal(0, nPts*3 - 1);
  //
  for ( int coordidx3d = 0, coordidx2d = 0; coordidx3d < coords3d->Length(); coordidx3d += 3, coordidx2d += 2 )
  {
    coords3d->ChangeValue(coordidx3d)     = coords2d->Value(coordidx2d);
    coords3d->ChangeValue(coordidx3d + 1) = coords2d->Value(coordidx2d + 1);
    coords3d->ChangeValue(coordidx3d + 2) = 0.0;
  }

  return new asiAlgo_PointCloud(coords3d);
}

//-----------------------------------------------------------------------------

//! Enumerates Data Parameters playing as sources for DOMAIN -> VTK
//! translation process.
//! \return source Parameters.
Handle(ActAPI_HParameterList) visu_iv_point_set_2d_data_provider::translationSources() const
{
  // Resulting Parameters
  ActParamStream out;

  Handle(asiData_IVPointSet2dNode)
    points_n = Handle(asiData_IVPointSet2dNode)::DownCast(m_node);
  //
  if ( points_n.IsNull() || !points_n->IsWellFormed() )
    return out;

  // Register Parameter as sensitive
  out << points_n->Parameter(asiData_IVPointSet2dNode::PID_Geometry);

  return out;
}