//-----------------------------------------------------------------------------
// Created on: 11 October 2018
// Created by: Edouard MEI
//-----------------------------------------------------------------------------
// Copyright (c) 2018 OPEN CASCADE SAS
//-----------------------------------------------------------------------------

#ifndef algoRepair_FindOpenEdges_h
#define algoRepair_FindOpenEdges_h

// algoFeat includes
#include <asiAlgo_AAG.h>

// Active Data includes
#include <ActAPI_IAlgorithm.h>

// OCCT includes
#include <NCollection_Vector.hxx>
#include <TopTools_HSequenceOfShape.hxx>

//-----------------------------------------------------------------------------

namespace cadpro {

//! \ingroup CADPRO_REPAIR
//!
//! This class detects open edges for the requested CAD model.
class algoRepair_FindOpenEdges : public ActAPI_IAlgorithm
{
public:

  // OCCT RTTI
  DEFINE_STANDARD_RTTI_INLINE(algoRepair_FindOpenEdges, ActAPI_IAlgorithm)

public:

  //! Constructor.
  //! \param[in] shape    working shape.
  //! \param[in] progress Progress Notifier.
  //! \param[in] plotter  Imperative Plotter.
  asiAlgo_EXPORT
    algoRepair_FindOpenEdges(const TopoDS_Shape&  shape,
                             ActAPI_ProgressEntry progress = NULL,
                             ActAPI_PlotterEntry  plotter = NULL);

  //! Find open edges.
  //! \return true in case of success, false -- otherwise.
  asiAlgo_EXPORT bool
    Perform();

public:

  //! Initializes open edges detector.
  void Init(const TopoDS_Shape& shape)
  {
    // Clean up and reinitialize inputs.
    m_shape = shape;
    m_aag = new asiAlgo_AAG(m_shape);

    // Clean up result.
    m_result.Clear();
  }

  //! \return detected open edges.
  const TopTools_IndexedMapOfShape& GetResultEdges() const
  {
    return m_result.edges;
  }

  //! \return indices of the detected open edges.
  const TColStd_PackedMapOfInteger& GetResultIndices() const
  {
    return m_result.ids;
  }

  //! Cleans the result.
  void ClearResult()
  {
    m_result.edges.Clear();
    m_result.ids.Clear();
  }

  //! Gets number of boundary wires.
  //! \return number of boundary wires.
  const int& GetNbBoundaryWires() const
  {
    return m_nbBoundaryWires;
  }

  //! Gets boundary contours.
  //! \return vector of boundary contours.
  const NCollection_Vector<TColStd_PackedMapOfInteger>& GetBoundaryContoursIds() const
  {
    return m_BoundaryContoursIds;
  }

  //! Gets boundary contours.
  //! \return vector of boundary contours.
  const NCollection_Vector<TopTools_ListOfShape>& GetBoundaryContours() const
  {
    return m_BoundaryContours;
  }

protected:

  //! Extracts IDs for the given topological edges.
  //! \param[in]  edges   edges to get the IDs for.
  //! \param[out] indices indices of the edges of interest.
  void getIds(const TopTools_IndexedMapOfShape& edges,
              TColStd_PackedMapOfInteger&       indices) const;

  //! Constructs loops from open edges.
  //! Loop is unordered set of edges which make a closed contour.
  void ConstructLoops();

protected:

  //-------------------------------------------------------------------------//
  // IN
  //-------------------------------------------------------------------------//

  TopoDS_Shape        m_shape; //!< Master CAD.
  Handle(asiAlgo_AAG) m_aag;   //!< Master AAG.

  //-------------------------------------------------------------------------//
  // OUT
  //-------------------------------------------------------------------------//

  struct
  {
    TopTools_IndexedMapOfShape edges; //!< Detected feature faces.
    TColStd_PackedMapOfInteger ids;   //!< Indices of the detected feature faces.

    void Clear()
    {
      edges.Clear();
      ids.Clear();
    }
  } m_result;

private:

  int m_nbBoundaryWires; //!< Number of border wires.
  NCollection_Vector<TColStd_PackedMapOfInteger> m_BoundaryContoursIds; //!< Border contours ids.
  NCollection_Vector<TopTools_ListOfShape> m_BoundaryContours; //!< Border contours.
};

};

#endif // algoRepair_FindOpenEdges_h
