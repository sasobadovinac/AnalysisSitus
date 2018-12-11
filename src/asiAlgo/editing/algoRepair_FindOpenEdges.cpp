//-----------------------------------------------------------------------------
// Created on: 11 October 2018
// Created by: Edouard MEI
//-----------------------------------------------------------------------------
// Copyright (c) 2016-2018 OPEN CASCADE SAS
//-----------------------------------------------------------------------------

// Own include
#include <algoRepair_FindOpenEdges.h>

// OCCT includes
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <NCollection_Vector.hxx>

//-----------------------------------------------------------------------------

cadpro::algoRepair_FindOpenEdges::algoRepair_FindOpenEdges(const TopoDS_Shape&  shape,
                                                           ActAPI_ProgressEntry progress,
                                                           ActAPI_PlotterEntry  plotter)
: ActAPI_IAlgorithm(progress, plotter)
{
  this->Init(shape);
}

//-----------------------------------------------------------------------------

bool cadpro::algoRepair_FindOpenEdges::Perform()
{
  // Build map of edges to extract open ("naked") ones.
  TopTools_IndexedDataMapOfShapeListOfShape edgesFaces;
  TopExp::MapShapesAndAncestors(m_shape, TopAbs_EDGE, TopAbs_FACE, edgesFaces);

  // Find open edges.
  Handle(TopTools_HSequenceOfShape) openEdgesSeq = new TopTools_HSequenceOfShape;
  for (int k = 1; k <= edgesFaces.Extent(); ++k)
  {
    // Progress.
    if (this->Progress().IsCancelling())
      return false;

    const TopTools_ListOfShape& faces = edgesFaces(k);
    if (faces.Extent() > 1)
      continue; // Skip non-open edges.

    const TopoDS_Edge& edge = TopoDS::Edge(edgesFaces.FindKey(k));
    if (BRep_Tool::Degenerated(edge))
      continue; // Skip degenerated.

    m_result.edges.Add(edge);
    openEdgesSeq->Append(edge);
  }

  // TODO: change to a) something easier b) that which doesn't change edges

  /*
  // Compose boundary wires from the naked edges.
  Handle(TopTools_HSequenceOfShape) boundaryWiresSeq = new TopTools_HSequenceOfShape;
  ShapeAnalysis_FreeBounds::ConnectEdgesToWires(openEdgesSeq, 1e-3, false, boundaryWiresSeq);
  m_nbBoundaryWires = boundaryWiresSeq->Size();
  */
  /*
  //  Check if edges are remained unchaged during 'wire-isation'

  int nbEdgesInWires = 0;
  int nbMatches = 0;

  for (TopTools_HSequenceOfShape::Iterator itW(*boundaryWiresSeq); itW.More(); itW.Next())
  {
    for (TopExp_Explorer expE(itW.Value(), TopAbs_EDGE); expE.More(); expE.Next())
    {
      TopoDS_Edge E = TopoDS::Edge(expE.Current());
      ++nbEdgesInWires;

      for (TopTools_HSequenceOfShape::Iterator itE(*openEdgesSeq); itE.More(); itE.Next())
      { 
        TopoDS_Edge E0 = TopoDS::Edge(itE.Value());
        
        if (E.IsSame(E0))
          ++nbMatches;
      }
    }
  }
  if (nbMatches == m_result.edges.Size())
  {
    cout << "mm nice";
  }

  //  End of Check
  //  result: 10 of 36 remains unchanged
  */


  /* no more boundary wires, onlu boundary contours
  m_BoundaryWires.Clear();
  
  TopTools_HSequenceOfShape::Iterator it((*boundaryWiresSeq), false);
  for (int ii = boundaryWiresSeq->Size()-1; it.More(); it.Next(), --ii)
    m_BoundaryWires.SetValue(ii, it.Value());
  */

  // Recover IDs.
  this->getIds(m_result.edges, m_result.ids);
  
  ConstructLoops();
  m_nbBoundaryWires = m_BoundaryContours.Size();

  return true;
}

//-----------------------------------------------------------------------------

void cadpro::algoRepair_FindOpenEdges::getIds(const TopTools_IndexedMapOfShape& edges,
                                              TColStd_PackedMapOfInteger&       indices) const
{
  const TopTools_IndexedMapOfShape& allEdges = m_aag->RequestMapOfEdges();
  for (int eidx = 1; eidx <= edges.Extent(); ++eidx)
  {
    const TopoDS_Shape& edgeOfInterest = edges(eidx);
    indices.Add(allEdges.FindIndex(edgeOfInterest));
  }
}

void cadpro::algoRepair_FindOpenEdges::ConstructLoops(void)
{
  // TODO: redo algoritm to make sure we take degenerated edges into account. Their existance are ignored for now

  m_BoundaryContours.Clear();
  m_BoundaryContoursIds.Clear();

  //! List of all vertices, included into edges
  TopTools_IndexedDataMapOfShapeListOfShape vertexEdges;
  TopExp::MapShapesAndAncestors(m_shape, TopAbs_VERTEX, TopAbs_EDGE, vertexEdges);

  //! A "closed" list of already used edges
  TColStd_PackedMapOfInteger processedEdges;

  //! A loop. An unordered set of edges which forms a closed contour.  
  TColStd_PackedMapOfInteger currentLoop;

  //! A way to get edge by it's index
  const TopTools_IndexedMapOfShape& allEdges = m_aag->RequestMapOfEdges();

  //! Last used vx
  TopoDS_Vertex currentVx;

  //! A current edge. Start from a random one
  TColStd_PackedMapOfInteger::Iterator currentIterator(m_result.ids);
  
  //! Algo
  int currentEdgeId = -1;
  bool done = false;
  while (!done)
  {
    // if we have no loops now - take the first non-processed edge from heap
    if (currentEdgeId < 0)
    {
      // clear the current loop
      currentLoop.Clear();

      // clear last vx
      currentVx.Nullify();

      // find first non-processed edge
      for (; currentIterator.More(); currentIterator.Next())
      {
        const int id = currentIterator.Key();
        if (!processedEdges.Contains(id))
        {
          currentEdgeId = id;
          break;
        }
      }

      // if no more "free" edges - we have all we can and can return found edge-loops
      if (currentEdgeId < 0)
      {
        done = true;
        continue;
        //>> TODO: no more edges, drop the algo and return result
      }
    }

    // get the current edge from id
    TopoDS_Edge currentEdge = TopoDS::Edge(allEdges.FindKey(currentEdgeId));
    
    // get the vertices of current edge
    TopoDS_Vertex vertices[2];
    TopExp::Vertices(currentEdge, vertices[0], vertices[1]);
    
    // no vx yet? get the random one
    if (currentVx.IsNull())
      currentVx = vertices[1];
    else
      currentVx = (currentVx.IsSame(vertices[1])) ? vertices[0] : vertices[1];

    int linkedIds[2]; // ids of open edges matching vertex is belong
    int linkedIndex;  // 
    int otherId;      // id of other edge, insted of current edge, current vertex is belong

    linkedIndex = 0;

    // get the edges which use that vertex 
    const TopTools_ListOfShape& edges = vertexEdges.FindFromKey(currentVx);
    TopTools_ListOfShape::Iterator it(edges);

    // loop throu edges, count opened
    for (; it.More(); it.Next())
    {
      const TopoDS_Shape& edgeOfInterest = it.Value();
      const int index = allEdges.FindIndex(edgeOfInterest);
      if (m_result.ids.Contains(index))
      {
        if (linkedIndex < 2)
          linkedIds[linkedIndex] = index;          
        ++linkedIndex;
      }
    }
    
    // if vertex used in more than 2 open edges - thats not a contour, drop it and start anew
    if (linkedIndex != 2)
    {
      processedEdges.Add(currentEdgeId);
      currentEdgeId = -1;
      continue;
    }
    
    // get the id of the other edge
    otherId = (currentEdgeId == linkedIds[0])? linkedIds[1]: linkedIds[0];
    processedEdges.Add(currentEdgeId);

    // final check - if other edge is already processed:
    // 1. we already used it in this contour. which means contour is closed
    // 2. we processed it another time and definetly something wrong now
    if (processedEdges.Contains(otherId))
    {
      if (currentLoop.Contains(otherId))
      {
        currentLoop.Add(currentEdgeId);
        m_BoundaryContoursIds.Append(currentLoop);
      }
      currentEdgeId = -1;
      continue;
    }

    // add current edge to the loop and mark as processed, set other edge as current
    currentLoop.Add(currentEdgeId);

    // select a next edge as current one
    currentEdgeId = otherId;
  }

  m_BoundaryContours.Clear();
  NCollection_Vector<TColStd_PackedMapOfInteger>::Iterator itVec(m_BoundaryContoursIds);
  for (; itVec.More(); itVec.Next())
  {
    TopTools_ListOfShape lst;
    TColStd_PackedMapOfInteger::Iterator itPack(itVec.Value());
    for (; itPack.More(); itPack.Next())
    {
      TopoDS_Edge e = TopoDS::Edge(allEdges.FindKey(itPack.Key()));
      lst.Append(e);
    }
    m_BoundaryContours.Append(lst);
  }

  return;
}