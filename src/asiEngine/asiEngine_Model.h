//-----------------------------------------------------------------------------
// Created on: 26 November 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

#ifndef asiEngine_Model_h
#define asiEngine_Model_h

// A-Situs includes
#include <asiEngine.h>

// A-Situs (data) includes
#include <asiData_BoundaryEdgesNode.h>
#include <asiData_IVCurveNode.h>
#include <asiData_IVCurvesNode.h>
#include <asiData_IVNode.h>
#include <asiData_IVPointSet2dNode.h>
#include <asiData_IVPointSetNode.h>
#include <asiData_IVPoints2dNode.h>
#include <asiData_IVPointsNode.h>
#include <asiData_IVSurfaceNode.h>
#include <asiData_IVSurfacesNode.h>
#include <asiData_IVTessItemNode.h>
#include <asiData_IVTessNode.h>
#include <asiData_IVTextItemNode.h>
#include <asiData_IVTextNode.h>
#include <asiData_IVTopoItemNode.h>
#include <asiData_IVTopoNode.h>
#include <asiData_Partition.h>
#include <asiData_PartNode.h>
#include <asiData_RootNode.h>

// Active Data includes
#include <ActData_BaseModel.h>

//-----------------------------------------------------------------------------

DEFINE_STANDARD_HANDLE(asiEngine_Model, ActData_BaseModel)

//! Standard implementation of Data Model in Analysis Situs. This
//! standard implementation answers many typical needs related to analysis.
//! However, in the prototyping scope you normally need something more.
//! If the basic data model is not enough, you may either extend it with your
//! custom Data Nodes or implement a brand new one. In the latter case you
//! cannot use some handy services from the Engine package (e.g.
//! initialization of Part Nodes).
class asiEngine_Model : public ActData_BaseModel
{
public:

  // OCCT RTTI
  DEFINE_STANDARD_RTTI_INLINE(asiEngine_Model, ActData_BaseModel)

public:

  asiEngine_EXPORT
    asiEngine_Model();

//-----------------------------------------------------------------------------
// Populate and clear:
public:

  asiEngine_EXPORT virtual void
    Populate();

  asiEngine_EXPORT virtual void
    Clear();

//-----------------------------------------------------------------------------
// Accessors to Nodes:
public:

  asiEngine_EXPORT Handle(asiData_PartNode)
    GetPartNode() const;

  asiEngine_EXPORT Handle(asiData_IVNode)
    GetIVNode() const;

//-----------------------------------------------------------------------------
// Overridden:
public:

  //! Create a cloned instance of Data Model.
  //! \return cloned instance.
  virtual Handle(ActAPI_IModel) Clone() const
  {
    return ActData_BaseModel::CloneInstance<asiEngine_Model>();
  }

//-----------------------------------------------------------------------------
// Partitions:
public:

  //! Accessor for a Partition instance dedicated to root Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_RootNode>) GetRootPartition() const
  {
    return Handle(asiData_Partition<asiData_RootNode>)::DownCast( this->Partition(Partition_Root) );
  }

//-----------------------------------------------------------------------------

  //! Accessor for a Partition instance dedicated to Geometry Part Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_PartNode>) GetPartPartition() const
  {
    return Handle(asiData_Partition<asiData_PartNode>)::DownCast( this->Partition(Partition_GeomPart) );
  }

  //! Accessor for a Partition instance dedicated to Geometry Face Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_FaceNode>) GetGeomFacePartition() const
  {
    return Handle(asiData_Partition<asiData_FaceNode>)::DownCast( this->Partition(Partition_GeomFace) );
  }

  //! Accessor for a Partition instance dedicated to Geometry Surface Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_SurfNode>) GetGeomSurfacePartition() const
  {
    return Handle(asiData_Partition<asiData_SurfNode>)::DownCast( this->Partition(Partition_GeomSurface) );
  }

  //! Accessor for a Partition instance dedicated to Geometry Edge Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_EdgeNode>) GetGeomEdgePartition() const
  {
    return Handle(asiData_Partition<asiData_EdgeNode>)::DownCast( this->Partition(Partition_GeomEdge) );
  }

  //! Accessor for a Partition instance dedicated to Geometry Curve Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_CurveNode>) GetGeomCurvePartition() const
  {
    return Handle(asiData_Partition<asiData_CurveNode>)::DownCast( this->Partition(Partition_GeomCurve) );
  }

  //! Accessor for a Partition instance dedicated to Geometry Boundary Edges Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_BoundaryEdgesNode>) GetGeomBoundaryEdgesPartition() const
  {
    return Handle(asiData_Partition<asiData_BoundaryEdgesNode>)::DownCast( this->Partition(Partition_GeomBoundaryEdges) );
  }

  //! Accessor for a Partition instance dedicated to Contour Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_ContourNode>) GetGeomContourPartition() const
  {
    return Handle(asiData_Partition<asiData_ContourNode>)::DownCast( this->Partition(Partition_GeomContour) );
  }

//-----------------------------------------------------------------------------

  //! Accessor for a Partition instance dedicated to IV Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVNode>) GetIVPartition() const
  {
    return Handle(asiData_Partition<asiData_IVNode>)::DownCast( this->Partition(Partition_IV) );
  }

  //! Accessor for a Partition instance dedicated to IV Points 2D Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVPoints2dNode>) GetIVPoints2dPartition() const
  {
    return Handle(asiData_Partition<asiData_IVPoints2dNode>)::DownCast( this->Partition(Partition_IV_Points2d) );
  }

  //! Accessor for a Partition instance dedicated to IV Points Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVPointsNode>) GetIVPointsPartition() const
  {
    return Handle(asiData_Partition<asiData_IVPointsNode>)::DownCast( this->Partition(Partition_IV_Points) );
  }

  //! Accessor for a Partition instance dedicated to IV Point Set 2D Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVPointSet2dNode>) GetIVPointSet2dPartition() const
  {
    return Handle(asiData_Partition<asiData_IVPointSet2dNode>)::DownCast( this->Partition(Partition_IV_PointSet2d) );
  }

  //! Accessor for a Partition instance dedicated to IV Point Set Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVPointSetNode>) GetIVPointSetPartition() const
  {
    return Handle(asiData_Partition<asiData_IVPointSetNode>)::DownCast( this->Partition(Partition_IV_PointSet) );
  }

  //! Accessor for a Partition instance dedicated to IV Curves Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVCurvesNode>) GetIVCurvesPartition() const
  {
    return Handle(asiData_Partition<asiData_IVCurvesNode>)::DownCast( this->Partition(Partition_IV_Curves) );
  }

  //! Accessor for a Partition instance dedicated to IV Curve Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVCurveNode>) GetIVCurvePartition() const
  {
    return Handle(asiData_Partition<asiData_IVCurveNode>)::DownCast( this->Partition(Partition_IV_Curve) );
  }

  //! Accessor for a Partition instance dedicated to IV Surfaces Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVSurfacesNode>) GetIVSurfacesPartition() const
  {
    return Handle(asiData_Partition<asiData_IVSurfacesNode>)::DownCast( this->Partition(Partition_IV_Surfaces) );
  }

  //! Accessor for a Partition instance dedicated to IV Surface Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVSurfaceNode>) GetIVSurfacePartition() const
  {
    return Handle(asiData_Partition<asiData_IVSurfaceNode>)::DownCast( this->Partition(Partition_IV_Surface) );
  }

  //! Accessor for a Partition instance dedicated to IV Topo Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTopoNode>) GetIVTopoPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTopoNode>)::DownCast( this->Partition(Partition_IV_Topo) );
  }

  //! Accessor for a Partition instance dedicated to IV Topo Item Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTopoItemNode>) GetIVTopoItemPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTopoItemNode>)::DownCast( this->Partition(Partition_IV_TopoItem) );
  }

  //! Accessor for a Partition instance dedicated to IV Tessellation Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTessNode>) GetIVTessPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTessNode>)::DownCast( this->Partition(Partition_IV_Tess) );
  }

  //! Accessor for a Partition instance dedicated to IV Tessellation Item Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTessItemNode>) GetIVTessItemPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTessItemNode>)::DownCast( this->Partition(Partition_IV_TessItem) );
  }

  //! Accessor for a Partition instance dedicated to IV Text Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTextNode>) GetIVTextPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTextNode>)::DownCast( this->Partition(Partition_IV_Text) );
  }

  //! Accessor for a Partition instance dedicated to IV Text Item Nodes.
  //! \return requested Partition.
  Handle(asiData_Partition<asiData_IVTextItemNode>) GetIVTextItemPartition() const
  {
    return Handle(asiData_Partition<asiData_IVTextItemNode>)::DownCast( this->Partition(Partition_IV_TextItem) );
  }

protected:

  asiEngine_EXPORT virtual void
    initPartitions();

  asiEngine_EXPORT virtual void
    initFunctionDrivers();

protected:

  asiEngine_EXPORT virtual Handle(ActAPI_IPartition)
    getVariablePartition(const VariableType& varType) const;

  asiEngine_EXPORT virtual Handle(ActAPI_INode)
    getRootNode() const;

  asiEngine_EXPORT virtual void
    invariantCopyRefs(ActAPI_FuncGUIDStream&         FuncGUIDs,
                      ActAPI_ParameterLocatorStream& Refs) const;

protected:

  asiEngine_EXPORT virtual int
    actualVersionApp();

  asiEngine_EXPORT virtual Handle(ActData_CAFConverter)
    converterApp();

protected:

  //! IDs of the registered Partitions.
  enum PartitionId
  {
    Partition_RealVar = 1, // Good practice to start with 1 when working with CAF
  //---------------------------------------------------------------------------
    Partition_Root,
  //---------------------------------------------------------------------------
    Partition_GeomPart,
    Partition_GeomFace,
    Partition_GeomSurface,
    Partition_GeomEdge,
    Partition_GeomCurve,
    Partition_GeomBoundaryEdges,
    Partition_GeomContour,
  //---------------------------------------------------------------------------
    Partition_IV,
    Partition_IV_Points2d,
    Partition_IV_Points,
    Partition_IV_PointSet2d,
    Partition_IV_PointSet,
    Partition_IV_Curves,
    Partition_IV_Curve,
    Partition_IV_Surfaces,
    Partition_IV_Surface,
    Partition_IV_Topo,
    Partition_IV_TopoItem,
    Partition_IV_Tess,
    Partition_IV_TessItem,
    Partition_IV_Text,
    Partition_IV_TextItem,
  //---------------------------------------------------------------------------
    Partition_LAST
  };

};

#endif