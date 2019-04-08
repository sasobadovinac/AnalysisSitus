//-----------------------------------------------------------------------------
// Created on: 24 August 2017
//-----------------------------------------------------------------------------
// Copyright (c) 2017, Sergey Slyadnev
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder(s) nor the
//      names of all contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//-----------------------------------------------------------------------------

// cmdEngine includes
#include <cmdEngine.h>

// asiTcl includes
#include <asiTcl_PluginMacro.h>

// asiEngine includes
#include <asiEngine_Curve.h>
#include <asiEngine_Editing.h>
#include <asiEngine_Part.h>
#include <asiEngine_TolerantShapes.h>

// asiAlgo includes
#include <asiAlgo_CheckDihedralAngle.h>
#include <asiAlgo_CheckValidity.h>
#include <asiAlgo_CompleteEdgeLoop.h>
#include <asiAlgo_MeshCheckInter.h>
#include <asiAlgo_MeshConvert.h>
#include <asiAlgo_RecognizeBlends.h>
#include <asiAlgo_Timer.h>
#include <asiAlgo_Utils.h>

// asiVisu includes
#include <asiVisu_Utils.h>

// asiUI includes
#include <asiUI_Plot2d.h>
#include <asiUI_DialogEuler.h>
#include <asiUI_IV.h>

#ifdef USE_MOBIUS
  #include <mobius/cascade_BSplineCurve.h>
  #include <mobius/cascade_BSplineSurface.h>
#endif

// OCCT includes
#include <BRep_Builder.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

//-----------------------------------------------------------------------------

int ENGINE_Explode(const Handle(asiTcl_Interp)& interp,
                   int                          argc,
                   const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Process input flags.
  const bool isNoLocation = interp->HasKeyword(argc, argv, "noloc");
  const bool isVertex     = interp->HasKeyword(argc, argv, "vertex");
  const bool isEdge       = interp->HasKeyword(argc, argv, "edge");
  const bool isWire       = interp->HasKeyword(argc, argv, "wire");
  const bool isFace       = interp->HasKeyword(argc, argv, "face");
  const bool isShell      = interp->HasKeyword(argc, argv, "shell");
  const bool isSolid      = interp->HasKeyword(argc, argv, "solid");

  if ( isNoLocation && (isVertex ||
                        isEdge   ||
                        isFace   ||
                        isShell  ||
                        isSolid) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot use '-noloc' flag with subshape qualifiers.");
    return TCL_ERROR;
  }

  // Disable visualization for efficiency.
  Handle(asiUI_IV)
    IV = Handle(asiUI_IV)::DownCast( interp->GetPlotter().Plotter() );
  //
  if ( !IV.IsNull() )
  {
    IV->BROWSER_OFF();
    IV->VISUALIZATION_OFF();
  }

  // Get Part Node and shape.
  Handle(asiData_PartNode) partNode = cmdEngine::model->GetPartNode();
  //
  TopoDS_Shape partShape = partNode->GetShape();

  if ( isNoLocation || (argc == 1) )
  {
    // Explode.
    for ( TopoDS_Iterator it(partShape, true, !isNoLocation); it.More(); it.Next() )
    {
      const TopoDS_Shape& subShape = it.Value();

      // Generate name.
      TCollection_AsciiString name = asiAlgo_Utils::ShapeTypeStr(subShape).c_str();

      // Draw imperatively (populates Data Model).
      interp->GetPlotter().DRAW_SHAPE(subShape, name);
    }
  }
  else
  {
    // Get qualifier of sub-shape.
    TopAbs_ShapeEnum subshapeType;

    if ( isVertex )
      subshapeType = TopAbs_VERTEX;
    //
    else if ( isEdge )
      subshapeType = TopAbs_EDGE;
    //
    else if ( isWire )
      subshapeType = TopAbs_WIRE;
    //
    else if ( isFace )
      subshapeType = TopAbs_FACE;
    //
    else if ( isShell )
      subshapeType = TopAbs_SHELL;
    //
    else if ( isSolid )
      subshapeType = TopAbs_SOLID;
    //
    else
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Unexpected sub-shape qualifier.");
      return TCL_ERROR;
    }

    // Explode.
    for ( TopExp_Explorer exp(partShape, subshapeType); exp.More(); exp.Next() )
    {
      const TopoDS_Shape& subShape = exp.Current();

      // Generate name.
      TCollection_AsciiString name = asiAlgo_Utils::ShapeTypeStr(subShape).c_str();

      // Draw imperatively (populates Data Model).
      interp->GetPlotter().DRAW_SHAPE(subShape, name);
    }
  }

  // Enable back UI updates.
  if ( !IV.IsNull() )
  {
    IV->BROWSER_ON();
    IV->VISUALIZATION_ON();
  }

  // Update browser.
  if ( !cmdEngine::cf.IsNull() && cmdEngine::cf->ObjectBrowser )
    cmdEngine::cf->ObjectBrowser->Populate();

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_PrintSummary(const Handle(asiTcl_Interp)& interp,
                        int                          argc,
                        const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node and shape.
  Handle(asiData_PartNode) partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part is not initialized.");
    return TCL_OK;
  }
  //
  TopoDS_Shape partShape = partNode->GetShape();

  // Build explicit topology graph
  asiAlgo_TopoGraph topograph(partShape);

  int numCompounds  = 0,
      numCompSolids = 0,
      numSolids     = 0,
      numShells     = 0,
      numFaces      = 0,
      numWires      = 0,
      numEdges      = 0,
      numVertices   = 0;
  //
  topograph.CalculateSummary(numCompounds,
                             numCompSolids,
                             numSolids,
                             numShells,
                             numFaces,
                             numWires,
                             numEdges,
                             numVertices);

  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "=============================================");
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Summary of part's contents");
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "---------------------------------------------");
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. compounds:  %1" << numCompounds);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. compsolids: %1" << numCompSolids);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. solids:     %1" << numSolids);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. shells:     %1" << numShells);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. faces:      %1" << numFaces);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. wires:      %1" << numWires);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. edges:      %1" << numEdges);
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\tNum. vertices:   %1" << numVertices);

  // Check naming.
  Handle(asiAlgo_Naming) naming = partNode->GetNaming();
  //
  if ( !naming.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "\t... Naming service is active.");
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_GetSummary(const Handle(asiTcl_Interp)& interp,
                      int                          argc,
                      const char**                 argv)
{
  if ( argc != 9 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node and shape.
  Handle(asiData_PartNode) partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part is not initialized.");
    return TCL_OK;
  }
  //
  TopoDS_Shape partShape = partNode->GetShape();

  // Build explicit topology graph
  asiAlgo_TopoGraph topograph(partShape);

  int numCompounds  = 0,
      numCompSolids = 0,
      numSolids     = 0,
      numShells     = 0,
      numFaces      = 0,
      numWires      = 0,
      numEdges      = 0,
      numVertices   = 0;
  //
  topograph.CalculateSummary(numCompounds,
                             numCompSolids,
                             numSolids,
                             numShells,
                             numFaces,
                             numWires,
                             numEdges,
                             numVertices);

  // Set Tcl variables.
  int varIdx = 0;
  //
  interp->SetVar<int>(argv[++varIdx], numCompounds);
  interp->SetVar<int>(argv[++varIdx], numCompSolids);
  interp->SetVar<int>(argv[++varIdx], numSolids);
  interp->SetVar<int>(argv[++varIdx], numShells);
  interp->SetVar<int>(argv[++varIdx], numFaces);
  interp->SetVar<int>(argv[++varIdx], numWires);
  interp->SetVar<int>(argv[++varIdx], numEdges);
  interp->SetVar<int>(argv[++varIdx], numVertices);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_FaceAddr(const Handle(asiTcl_Interp)& interp,
                    int                          argc,
                    const char**                 argv)
{
  if ( argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  const int fidx = atoi(argv[1]);
  //
  if ( fidx < 1 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Face index should be 1-based.");
    return TCL_OK;
  }

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Get face.
  const TopoDS_Face& face = part_n->GetAAG()->GetFace(fidx);

  // Output.
  interp->GetProgress().SendLogMessage( LogInfo(Normal) << "Face %1: %2"
                                                        << fidx
                                                        << asiAlgo_Utils::ShapeAddr(face).c_str() );
  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckDist(const Handle(asiTcl_Interp)& interp,
                     int                          argc,
                     const char**                 argv)
{
  if ( argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Get topological item to imprint.
  Handle(asiData_IVTopoItemNode)
    topoItem_n = Handle(asiData_IVTopoItemNode)::DownCast( cmdEngine::model->FindNodeByName(argv[1]) );
  //
  if ( topoItem_n.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot find topological object with name %1." << argv[1]);
    return TCL_OK;
  }

  TIMER_NEW
  TIMER_GO

  BRepExtrema_DistShapeShape distSS(part_n->GetShape(),
                                    topoItem_n->GetShape(),
                                    Extrema_ExtFlag_MIN);

  TIMER_FINISH
  TIMER_COUT_RESULT_NOTIFIER(interp->GetProgress(), "Shape-shape distance")

  if ( !distSS.IsDone() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Extrema not done.");
    return TCL_OK;
  }

  // Get number of solutions.
  const int numSol = distSS.NbSolution();
  //
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Number of extrema solutions: %1." << numSol);

  for ( int k = 1; k <= numSol; ++k )
  {
    const gp_Pnt& P1 = distSS.PointOnShape1(k);
    const gp_Pnt& P2 = distSS.PointOnShape2(k);

    interp->GetPlotter().DRAW_POINT(P1, Color_Red, "dist_P1");
    interp->GetPlotter().DRAW_POINT(P2, Color_Red, "dist_P2");
    interp->GetPlotter().DRAW_LINK(P1, P2, Color_Red, "dist_P1P2");
  }

  interp->GetProgress().SendLogMessage( LogInfo(Normal) << "Distance: %1."
                                                        << distSS.Value() );

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckCurvature(const Handle(asiTcl_Interp)& interp,
                          int                          argc,
                          const char**                 argv)
{
  if ( argc > 6 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  /* ================
   *  Prepare inputs
   * ================ */

  // Find control keys.
  bool noPlot           = false;
  bool noAlongCurvature = false;
  //
  for ( int i = 1; i < argc; ++i )
  {
    if ( interp->IsKeyword(argv[i], "noplot") )
      noPlot = true;

    if ( interp->IsKeyword(argv[i], "noalong") )
      noAlongCurvature = true;
  }

  // Get numerical values.
  int    numPts      = 100;
  double scaleFactor = 1.0;
  double amplFactor  = 1.0;
  //
  if ( argc >= 2 )
    numPts = atoi(argv[1]);
  if ( argc >= 3 )
    scaleFactor = atof(argv[2]);
  if ( argc >= 4 )
    amplFactor = atof(argv[3]);

  // Get Part Node to access the selected edge.
  Handle(asiData_PartNode) partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape                      partShape = partNode->GetShape();
  const TopTools_IndexedMapOfShape& subShapes = partNode->GetAAG()->RequestMapOfSubShapes();

  // Curve Node is expected.
  Handle(asiData_CurveNode) curveNode = partNode->GetCurveRepresentation();
  //
  if ( curveNode.IsNull() || !curveNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Curve Node is null or ill-defined.");
    return TCL_OK;
  }

  // Get ID of the selected edge.
  const int edgeIdx = curveNode->GetSelectedEdge();
  //
  if ( edgeIdx <= 0 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Please, select edge first.");
    return TCL_OK;
  }

  // Get host curve of the selected edge.
  const TopoDS_Shape& edgeShape = subShapes(edgeIdx);
  //
  if ( edgeShape.ShapeType() != TopAbs_EDGE )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Unexpected topological type of the selected edge.");
    return TCL_OK;
  }
  //
  double f, l;
  Handle(Geom_Curve) curve = BRep_Tool::Curve( TopoDS::Edge(edgeShape), f, l );

  /* ==========================
   *  Evaluate curvature combs
   * ========================== */

  // Create curvature combs.
  asiEngine_Curve CurveAPI( cmdEngine::model, interp->GetProgress(), interp->GetPlotter() );
  //
  Handle(asiData_CurvatureCombsNode) combsNode;
  //
  std::vector<gp_Pnt> points;
  std::vector<double> params;
  std::vector<double> curvatures;
  std::vector<gp_Vec> combs;
  std::vector<bool>   combsOk;
  //
  cmdEngine::model->OpenCommand();
  {
    // Calculate curvature field.
    if ( !asiAlgo_Utils::CalculateCurvatureCombs(curve,
                                                 f,
                                                 l,
                                                 numPts,
                                                 amplFactor,
                                                 points,
                                                 params,
                                                 curvatures,
                                                 combs,
                                                 combsOk) )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot calculate curvature field.");
      return TCL_OK;
    }

    // Create persistent object.
    combsNode = CurveAPI.CreateOrUpdateCurvatureCombs(curveNode,
                                                      scaleFactor,
                                                      points,
                                                      combsOk,
                                                      params,
                                                      curvatures,
                                                      combs);
  }
  cmdEngine::model->CommitCommand();

  // Actualize.
  if ( cmdEngine::cf )
  {
    cmdEngine::cf->ObjectBrowser->Populate();
    cmdEngine::cf->ViewerPart->PrsMgr()->Actualize(combsNode);
  }

  /* ======================
   *  Build curvature plot
   * ====================== */

  if ( cmdEngine::cf && !noPlot )
  {
    // Open curvature plot.
    asiUI_Plot2d*
      cPlot = new asiUI_Plot2d( interp->GetProgress(),
                                interp->GetPlotter() );
    //
    cPlot->Render(params, curvatures, "Parameter", "Curvature", "Curvature Plot");
  }

  /* ========================================
   *  Check along-curvature (if faces exist)
   * ======================================== */

  if ( !noAlongCurvature )
  {
    // Get owner faces.
    TopTools_IndexedDataMapOfShapeListOfShape M;
    TopExp::MapShapesAndAncestors(partShape, TopAbs_EDGE, TopAbs_FACE, M);
    //
    const TopTools_ListOfShape& edgeFaces = M.FindFromKey(edgeShape);
    //
    for ( TopTools_ListIteratorOfListOfShape fit(edgeFaces); fit.More(); fit.Next() )
    {
      const TopoDS_Face& faceShape = TopoDS::Face( fit.Value() );

      // Evaluate along curvature.
      double k;
      if ( !asiAlgo_Utils::EvaluateAlongCurvature(faceShape, TopoDS::Edge(edgeShape), k) )
      {
        interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot evaluate along-curvature.");
        continue;
      }
      else
      {
        interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Along curvature: %1." << k);

        if ( Abs(k) < 1.e-5 )
          interp->GetProgress().SendLogMessage( LogInfo(Normal) << "Curvature radius is infinite." );
        else
          interp->GetProgress().SendLogMessage( LogInfo(Normal) << "Curvature radius: %1." << Abs(1.0 / k) );
      }
    }
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckContinuity(const Handle(asiTcl_Interp)& interp,
                           int                          argc,
                           const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  /* ================
   *  Prepare inputs
   * ================ */

  // Get Part Node to access the selected face.
  Handle(asiData_PartNode) partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape                      partShape = partNode->GetShape();
  const TopTools_IndexedMapOfShape& subShapes = partNode->GetAAG()->RequestMapOfSubShapes();

  // Surf Node is expected.
  Handle(asiData_SurfNode) surfNode = partNode->GetSurfaceRepresentation();
  //
  if ( surfNode.IsNull() || !surfNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Surface Node is null or ill-defined.");
    return TCL_OK;
  }

  // Get ID of the selected face.
  const int faceIdx = surfNode->GetSelectedFace();
  //
  if ( faceIdx <= 0 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Please, select face first.");
    return TCL_OK;
  }

  // Get host surface of the selected face.
  const TopoDS_Shape& faceShape = subShapes(faceIdx);
  //
  if ( faceShape.ShapeType() != TopAbs_FACE )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Unexpected topological type of the selected face.");
    return TCL_OK;
  }
  //
  Handle(Geom_Surface)
    surface = BRep_Tool::Surface( TopoDS::Face(faceShape) );

  // Only B-surfaces are allowed.
  Handle(Geom_BSplineSurface)
    bsurf = Handle(Geom_BSplineSurface)::DownCast(surface);
  //
  if ( bsurf.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Only B-surfaces are allowed in this function.");
    return TCL_OK;
  }

  /* =====================
   *  Evaluate continuity
   * ===================== */

  const TColStd_Array1OfReal&    uKnots = bsurf->UKnots();
  const TColStd_Array1OfReal&    vKnots = bsurf->VKnots();
  const TColStd_Array1OfInteger& uMults = bsurf->UMultiplicities();
  const TColStd_Array1OfInteger& vMults = bsurf->VMultiplicities();

  // Draw U defects.
  for ( int i = 1; i <= uKnots.Length(); ++i )
  {
    const double u     = uKnots(i);
    const int    uMult = uMults(i);

    if ( uMult > 1 )
    {
      Handle(Geom_Curve) iso = bsurf->UIso(u);

      TCollection_AsciiString isoName("mult=");
      isoName += uMult;
      isoName += " u=";
      isoName += u;

      interp->GetPlotter().REDRAW_CURVE(isoName, iso, Color_Red);
    }
  }

  // Draw V defects.
  for ( int i = 1; i <= vKnots.Length(); ++i )
  {
    const double v     = vKnots(i);
    const int    vMult = vMults(i);

    if ( vMult > 1 )
    {
      Handle(Geom_Curve) iso = bsurf->VIso(v);

      TCollection_AsciiString isoName("mult=");
      isoName += vMult;
      isoName += " v=";
      isoName += v;

      interp->GetPlotter().REDRAW_CURVE(isoName, iso, Color_Red);
    }
  }

  // TODO

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckEuler(const Handle(asiTcl_Interp)& interp,
                      int                          argc,
                      const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  if ( cmdEngine::cf && (argc == 1) )
  {
    asiUI_DialogEuler* pEuler = new asiUI_DialogEuler( cmdEngine::model,
                                                       interp->GetProgress() );
    //
    pEuler->show();
  }
  else
  {
    const int genus = atoi(argv[1]);

    // Calculate the Euler-Poincare property for the active part.
    const bool
      isOk = asiEngine_Editing( cmdEngine::model,
                                interp->GetProgress(),
                                interp->GetPlotter() ).CheckEulerPoincare(genus);
    //
    if ( !isOk )
      interp->GetProgress().SendLogMessage(LogWarn(Normal) << "Euler-Poincare property does not hold.");

    // Append result to the interpreter.
    *interp << isOk;
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_DumpTopographDot(const Handle(asiTcl_Interp)& interp,
                            int                          argc,
                            const char**                 argv)
{
  if ( argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Access or create topograph.
  Handle(asiAlgo_TopoGraph) topograph;
  //
  if ( part_n->HasNaming() )
    topograph = part_n->GetNaming()->GetTopoGraph(); // With names.
  else
    topograph = new asiAlgo_TopoGraph( part_n->GetShape() ); // Without names.

  // Dump to buffer.
  std::stringstream buff;
  topograph->Dump(buff);

  // Dump to file.
  std::ofstream FILE(argv[1]);
  if ( !FILE.is_open() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot open file for writing.");
    return TCL_OK;
  }
  //
  FILE << buff.str().c_str();
  //
  FILE.close();

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_EvalCurve(const Handle(asiTcl_Interp)& interp,
                     int                          argc,
                     const char**                 argv)
{
  if ( argc < 4 || argc > 5 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Find Curve Node by name.
  Handle(ActAPI_INode) node = cmdEngine::model->FindNodeByName(argv[1]);
  //
  if ( node.IsNull() || !node->IsKind( STANDARD_TYPE(asiData_IVCurveNode) ) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Node '%1' is not a curve."
                                                        << argv[1]);
    return TCL_OK;
  }
  //
  Handle(asiData_IVCurveNode)
    curveNode = Handle(asiData_IVCurveNode)::DownCast(node);

  // Get curve.
  double f, l;
  Handle(Geom_Curve) occtCurve = curveNode->GetCurve(f, l);
  //
  if ( occtCurve.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "The curve in question is NULL.");
    return TCL_OK;
  }

  // Get parameter value.
  const double u = atof(argv[2]);

  // Get order.
  const int order = atoi(argv[3]);
  //
  if ( order < 0 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Negative order is meaningless.");
    return TCL_OK;
  }

  // Check whether Mobius evaluation is requested.
  bool isMobius = interp->HasKeyword(argc, argv, "mobius");

  // Evaluate curve.
  TCollection_AsciiString Message("\n========================================");
  Message                      += "\n Curve evaluation";
  Message                      += "\n----------------------------------------";
  //
  if ( !isMobius )
  {
    Message += "\nMethod: OpenCascade";

    if ( order == 0 )
    {
      gp_Pnt eval_P;
      occtCurve->D0(u, eval_P);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 1 )
    {
      gp_Pnt eval_P;
      gp_Vec eval_D1;
      occtCurve->D1(u, eval_P, eval_D1);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1", eval_P, eval_D1, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 2 )
    {
      gp_Pnt eval_P;
      gp_Vec eval_D1, eval_D2;
      occtCurve->D2(u, eval_P, eval_D1, eval_D2);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1", eval_P, eval_D1, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2", eval_P, eval_D2, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 3 )
    {
      gp_Pnt eval_P;
      gp_Vec eval_D1, eval_D2, eval_D3;
      occtCurve->D3(u, eval_P, eval_D1, eval_D2, eval_D3);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1", eval_P, eval_D1, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2", eval_P, eval_D2, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D3", eval_P, eval_D3, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD3(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.Z() ).c_str();
      Message += ")";
    }
    else
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Order %1 is not supported for OpenCascade method."
                                                          << order);
      return TCL_OK;
    }
  }
  else
  {
#ifndef USE_MOBIUS
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Mobius module is disabled.");
    return TCL_ERROR;
#else
    Message += "\nMethod: Mobius";

    Handle(Geom_BSplineCurve)
      occtBCurve = Handle(Geom_BSplineCurve)::DownCast(occtCurve);
    //
    if ( occtBCurve.IsNull() )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "The curve in question is not a B-curve.");
      return TCL_OK;
    }

    // Convert to Mobius curve.
    mobius::cascade_BSplineCurve converter(occtBCurve);
    converter.DirectConvert();
    //
    const mobius::ptr<mobius::bcurve>&
      mobCurve = converter.GetMobiusCurve();

    // Evaluate.
    if ( order == 0 )
    {
      mobius::xyz eval_P;
      mobCurve->Eval(u, eval_P);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                        gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                        Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 1 )
    {
      mobius::xyz eval_P;
      mobCurve->Eval(u, eval_P);
      //
      mobius::xyz eval_D1;
      mobCurve->Eval_Dk(u, 1, eval_D1);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                        gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                        Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D1.X(), eval_D1.Y(), eval_D1.Z() ),
                                            Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 2 )
    {
      mobius::xyz eval_P;
      mobCurve->Eval(u, eval_P);
      //
      mobius::xyz eval_D1;
      mobCurve->Eval_Dk(u, 1, eval_D1);
      //
      mobius::xyz eval_D2;
      mobCurve->Eval_Dk(u, 2, eval_D2);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                        gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                        Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D1.X(), eval_D1.Y(), eval_D1.Z() ),
                                            Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D2.X(), eval_D2.Y(), eval_D2.Z() ),
                                            Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 3 )
    {
      mobius::xyz eval_P;
      mobCurve->Eval(u, eval_P);
      //
      mobius::xyz eval_D1;
      mobCurve->Eval_Dk(u, 1, eval_D1);
      //
      mobius::xyz eval_D2;
      mobCurve->Eval_Dk(u, 2, eval_D2);
      //
      mobius::xyz eval_D3;
      mobCurve->Eval_Dk(u, 3, eval_D3);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                        gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                        Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D1.X(), eval_D1.Y(), eval_D1.Z() ),
                                            Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D2.X(), eval_D2.Y(), eval_D2.Z() ),
                                            Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D3",
                                            gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                            gp_Vec( eval_D3.X(), eval_D3.Y(), eval_D3.Z() ),
                                            Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD3(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D3.Z() ).c_str();
      Message += ")";
    }
    else
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Order %1 is not supported for Mobius method."
                                                          << order);
      return TCL_OK;
    }
#endif
  }
  Message += "\n----------------------------------------";
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << Message);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_EvalSurf(const Handle(asiTcl_Interp)& interp,
                    int                          argc,
                    const char**                 argv)
{
  if ( argc < 5 || argc > 6 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Find Surface Node by name.
  Handle(ActAPI_INode) node = cmdEngine::model->FindNodeByName(argv[1]);
  //
  if ( node.IsNull() || !node->IsKind( STANDARD_TYPE(asiData_IVSurfaceNode) ) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Node '%1' is not a surface."
                                                        << argv[1]);
    return TCL_OK;
  }
  //
  Handle(asiData_IVSurfaceNode)
    surfNode = Handle(asiData_IVSurfaceNode)::DownCast(node);

  // Get surface.
  Handle(Geom_Surface) occtSurface = surfNode->GetSurface();
  //
  if ( occtSurface.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "The surface in question is NULL.");
    return TCL_OK;
  }

  // Get parameter values.
  const double u = atof(argv[2]);
  const double v = atof(argv[3]);

  // Get order.
  const int order = atoi(argv[4]);
  //
  if ( order < 0 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Negative order is meaningless.");
    return TCL_OK;
  }

  // Check whether Mobius evaluation is requested.
  bool isMobius = interp->HasKeyword(argc, argv, "mobius");

  // Evaluate surface.
  TCollection_AsciiString Message("\n========================================");
  Message                      += "\n Surface evaluation";
  Message                      += "\n----------------------------------------";
  //
  if ( !isMobius )
  {
    Message += "\nMethod: OpenCascade";

    if ( order == 0 )
    {
      gp_Pnt eval_P;
      occtSurface->D0(u, v, eval_P);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 1 )
    {
      gp_Pnt eval_P;
      gp_Vec eval_D1u, eval_D1v;
      occtSurface->D1(u, v, eval_P, eval_D1u, eval_D1v);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1u", eval_P, eval_D1u, Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1v", eval_P, eval_D1v, Color_Green);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 2 )
    {
      gp_Pnt eval_P;
      gp_Vec eval_D1u, eval_D1v, eval_D2u, eval_D2v, eval_D2uv;
      occtSurface->D2(u, v, eval_P, eval_D1u, eval_D1v, eval_D2u, eval_D2v, eval_D2uv);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P", eval_P, Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1u", eval_P, eval_D1u, Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1v", eval_P, eval_D1v, Color_Green);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2u", eval_P, eval_D2u, Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2v", eval_P, eval_D2v, Color_Green);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2uv", eval_P, eval_D2uv, Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2uv(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.Z() ).c_str();
      Message += ")";
    }
    else
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Order %1 is not supported for OpenCascade method."
                                                          << order);
      return TCL_OK;
    }
  }
  else
  {
#ifndef USE_MOBIUS
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Mobius module is disabled.");
    return TCL_ERROR;
#else
    Message += "\nMethod: Mobius";

    Handle(Geom_BSplineSurface)
      occtBSurface = Handle(Geom_BSplineSurface)::DownCast(occtSurface);
    //
    if ( occtBSurface.IsNull() )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "The surface in question is not a B-surface.");
      return TCL_OK;
    }

    // Convert to Mobius surface.
    mobius::cascade_BSplineSurface converter(occtBSurface);
    converter.DirectConvert();
    //
    const mobius::ptr<mobius::bsurf>&
      mobSurface = converter.GetMobiusSurface();

    // Evaluate.
    if ( order == 0 )
    {
      mobius::xyz eval_P;
      mobSurface->Eval(u, v, eval_P);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                        gp_Pnt( eval_P.X(), eval_P.Y(), eval_P.Z() ),
                                        Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 1 )
    {
      mobius::xyz eval_P, eval_D1u, eval_D1v;
      mobSurface->Eval_D1(u, v, eval_P, eval_D1u, eval_D1v);
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                         gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                         Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1u",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D1u.X(), eval_D1u.Y(), eval_D1u.Z()),
                                             Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1v",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D1v.X(), eval_D1v.Y(), eval_D1v.Z()),
                                             Color_Green);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Z() ).c_str();
      Message += ")";
    }
    else if ( order == 2 )
    {
      mobius::xyz eval_P, eval_D1u, eval_D1v, eval_D2u, eval_D2v, eval_D2uv;
      mobSurface->Eval_D2(u, v, eval_P, eval_D1u, eval_D1v, eval_D2u, eval_D2v, eval_D2uv);
      //
      //
      interp->GetPlotter().REDRAW_POINT("eval_P",
                                         gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                         Color_Yellow);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1u",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D1u.X(), eval_D1u.Y(), eval_D1u.Z()),
                                             Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D1v",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D1v.X(), eval_D1v.Y(), eval_D1v.Z()),
                                             Color_Green);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2u",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D2u.X(), eval_D2u.Y(), eval_D2u.Z()),
                                             Color_Red);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2v",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D2v.X(), eval_D2v.Y(), eval_D2v.Z()),
                                             Color_Green);
      interp->GetPlotter().REDRAW_VECTOR_AT("eval_D2uv",
                                             gp_Pnt(eval_P.X(), eval_P.Y(), eval_P.Z()),
                                             gp_Vec(eval_D2uv.X(), eval_D2uv.Y(), eval_D2uv.Z()),
                                             Color_Yellow);
      //
      Message += "\nP(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_P.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD1v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D1v.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2u(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2u.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2v(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2v.Z() ).c_str();
      Message += ")";
      //
      Message += "\nD2uv(X, Y, Z): (";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.X() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.Y() ).c_str();
      Message += ", ";
      Message += asiAlgo_Utils::Str::ToString<double>( eval_D2uv.Z() ).c_str();
      Message += ")";
    }
    else
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Order %1 is not supported for Mobius method."
                                                          << order);
      return TCL_OK;
    }
#endif
  }
  Message += "\n----------------------------------------";
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << Message);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckToler(const Handle(asiTcl_Interp)& interp,
                      int                          argc,
                      const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  const int numRanges = ( (argc == 2) ? atoi(argv[1]) : 10 );

  Handle(asiEngine_Model)
    M = Handle(asiEngine_Model)::DownCast( interp->GetModel() );

  // Attempt to get the highlighted sub-shapes. If no sub-shapes are selected,
  // we work on the full part.
  TopTools_IndexedMapOfShape selectedSubShapes;
  //
  if ( cmdEngine::cf->ViewerPart )
  {
    asiEngine_Part PartAPI( M,
                            cmdEngine::cf->ViewerPart->PrsMgr(),
                            interp->GetProgress(),
                            interp->GetPlotter() );
    //
    PartAPI.GetHighlightedSubShapes(selectedSubShapes);
  }

  // Choose shape to check.
  TopoDS_Shape shape2Check;
  //
  if ( selectedSubShapes.Extent() )
  {
    TopoDS_Compound comp;
    BRep_Builder().MakeCompound(comp);
    //
    for ( int k = 1; k <= selectedSubShapes.Extent(); ++k )
      BRep_Builder().Add( comp, selectedSubShapes(k) );

    shape2Check = comp;
  }
  else
  {
    shape2Check = M->GetPartNode()->GetShape();
  }

  // Prepare API to analyze tolerances.
  asiEngine_TolerantShapes TolInfo( M,
                                    cmdEngine::cf->ViewerPart->PrsMgr(),
                                    interp->GetProgress(),
                                    interp->GetPlotter() );

  // Perform tolerance analysis.
  M->OpenCommand();
  {
    TolInfo.Clean_All();
    TolInfo.Populate(shape2Check, numRanges);
  }
  M->CommitCommand();

  // Update UI.
  if ( cmdEngine::cf )
    cmdEngine::cf->ObjectBrowser->Populate();

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckLength(const Handle(asiTcl_Interp)& interp,
                       int                          argc,
                       const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  Handle(asiEngine_Model)
    M = Handle(asiEngine_Model)::DownCast( interp->GetModel() );

  // Attempt to get the highlighted sub-shapes.
  TColStd_PackedMapOfInteger selectedEdgeIds;
  //
  if ( cmdEngine::cf && cmdEngine::cf->ViewerPart )
  {
    asiEngine_Part PartAPI( M,
                            cmdEngine::cf->ViewerPart->PrsMgr(),
                            interp->GetProgress(),
                            interp->GetPlotter() );
    //
    PartAPI.GetHighlightedEdges(selectedEdgeIds);
  }

  // Get total length.
  double len = 0.0;
  for ( TColStd_MapIteratorOfPackedMapOfInteger eit(selectedEdgeIds); eit.More(); eit.Next() )
  {
    const int edgeId = eit.Key();

    // Get edge.
    const TopoDS_Shape&
      edge = M->GetPartNode()->GetAAG()->RequestMapOfEdges()(edgeId);

    // Calculate global properties.
    GProp_GProps props;
    BRepGProp::LinearProperties(edge, props);
    len += props.Mass();
  }

  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Length: %1." << len);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckArea(const Handle(asiTcl_Interp)& interp,
                     int                          argc,
                     const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  Handle(asiEngine_Model)
    M = Handle(asiEngine_Model)::DownCast( interp->GetModel() );

  // Attempt to get the highlighted sub-shapes.
  TColStd_PackedMapOfInteger selectedFaceIds;
  //
  if ( cmdEngine::cf && cmdEngine::cf->ViewerPart )
  {
    asiEngine_Part PartAPI( M,
                            cmdEngine::cf->ViewerPart->PrsMgr(),
                            interp->GetProgress(),
                            interp->GetPlotter() );
    //
    PartAPI.GetHighlightedFaces(selectedFaceIds);
  }

  // Get total area.
  double area = 0.0;
  for ( TColStd_MapIteratorOfPackedMapOfInteger fit(selectedFaceIds); fit.More(); fit.Next() )
  {
    const int faceId = fit.Key();

    // Get face.
    const TopoDS_Shape&
      face = M->GetPartNode()->GetAAG()->GetMapOfFaces()(faceId);

    // Calculate global properties.
    GProp_GProps props;
    BRepGProp::SurfaceProperties(face, props);
    area += props.Mass();
  }

  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Area: %1." << area);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckValidity(const Handle(asiTcl_Interp)& interp,
                         int                          argc,
                         const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  asiAlgo_CheckValidity checker( interp->GetProgress(), interp->GetPlotter() );

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Get Part shape.
  TopoDS_Shape partSh = part_n->GetShape();

  if ( !checker.CheckBasic(partSh) )
  {
    *interp << 0;
  }
  else
  {
    *interp << 1;
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckFinite(const Handle(asiTcl_Interp)& interp,
                       int                          argc,
                       const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Get Part shape.
  TopoDS_Shape partSh = part_n->GetShape();

  // Get all contained solids.
  TopTools_IndexedMapOfShape solids;
  TopExp::MapShapes(partSh, TopAbs_SOLID, solids);
  //
  if ( solids.IsEmpty() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "There are no solids in the part.");
    return TCL_ERROR;
  }

  // Check contained solids.
  bool isOneSolid = (solids.Extent() == 1);
  bool areAllOk   = true;
  int  solidIdx   = 0;
  //
  for ( int k = 1; k <= solids.Extent(); ++k )
  {
    const TopoDS_Solid& solid = TopoDS::Solid( solids(k) );
    solidIdx++;

    // Check next solid.
    const bool isNextOk = asiAlgo_CheckValidity::IsFinite(solid);
    //
    if ( !isNextOk && areAllOk )
      areAllOk = false;

    // Print message.
    if ( isOneSolid )
    {
      if ( isNextOk )
        interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Solid is finite.");
      else
        interp->GetProgress().SendLogMessage(LogWarn(Normal) << "Solid is infinite.");
    }
    else
    {
      if ( isNextOk )
        interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Solid %1 is finite."
                                                             << solidIdx);
      else
        interp->GetProgress().SendLogMessage(LogWarn(Normal) << "Solid %1 is infinite."
                                                             << solidIdx);
    }
  }

  // Add checker result to the interpretor.
  *interp << (areAllOk ? 1 : 0);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckContours(const Handle(asiTcl_Interp)& interp,
                         int                          argc,
                         const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Read tolerance.
  double globTolerance = ( (argc == 2) ? atof(argv[1]) : 0.0 );

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Get Part shape.
  TopoDS_Shape partSh = part_n->GetShape();

  // Check each face individually.
  bool isOk = true;
  //
  for ( TopExp_Explorer exp(partSh, TopAbs_FACE); exp.More(); exp.Next() )
  {
    const TopoDS_Face& face = TopoDS::Face( exp.Current() );

    // Set default tolerance.
    double locTolerance;
    if ( globTolerance )
      locTolerance = globTolerance;
    else
      locTolerance = asiAlgo_CheckValidity::MaxTolerance(face)*5.0;

    // Check closeness.
    if ( !asiAlgo_CheckValidity::HasAllClosedWires(face, locTolerance) )
    {
      isOk = false;
      break;
    }
  }

  if ( isOk )
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "All faces have closed contours.");
  else
    interp->GetProgress().SendLogMessage(LogWarn(Normal) << "Some contours contain gaps.");

  *interp << (isOk ? 1 : 0);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_GetTolerance(const Handle(asiTcl_Interp)& interp,
                        int                          argc,
                        const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node.
  Handle(asiData_PartNode) part_n = cmdEngine::model->GetPartNode();

  // Return max tolerance to the interpreter.
  *interp << asiAlgo_CheckValidity::MaxTolerance( part_n->GetShape() );

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_GetCurveStrain(const Handle(asiTcl_Interp)& interp,
                          int                          argc,
                          const char**                 argv)
{
  if ( argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Find Curve Node by name.
  Handle(ActAPI_INode) node = cmdEngine::model->FindNodeByName(argv[1]);
  //
  if ( node.IsNull() || !node->IsKind( STANDARD_TYPE(asiData_IVCurveNode) ) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Node '%1' is not a curve."
                                                        << argv[1]);
    return TCL_OK;
  }
  //
  Handle(asiData_IVCurveNode)
    curveNode = Handle(asiData_IVCurveNode)::DownCast(node);

  // Get curve.
  double f, l;
  Handle(Geom_Curve) occtCurve = curveNode->GetCurve(f, l);
  //
  if ( occtCurve.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "The curve in question is NULL.");
    return TCL_OK;
  }

  // Calculate strain energy.
  double energy = 0;
  if ( !asiAlgo_Utils::CalculateStrainEnergy(occtCurve, energy) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot calculate strain energy.");
    return TCL_OK;
  }

  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Approximate strain energy: %1."
                                                       << energy);

  // Add to interpreter.
  *interp << energy;

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_GetSurfaceBending(const Handle(asiTcl_Interp)& interp,
                             int                          argc,
                             const char**                 argv)
{
  if ( argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Find Surface Node by name.
  Handle(ActAPI_INode) node = cmdEngine::model->FindNodeByName(argv[1]);
  //
  if ( node.IsNull() || !node->IsKind( STANDARD_TYPE(asiData_IVSurfaceNode) ) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Node '%1' is not a surface."
                                                        << argv[1]);
    return TCL_OK;
  }
  //
  Handle(asiData_IVSurfaceNode)
    surfNode = Handle(asiData_IVSurfaceNode)::DownCast(node);

  // Get surface.
  Handle(Geom_Surface) occtSurf = surfNode->GetSurface();
  //
  if ( occtSurf.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "The surface in question is NULL.");
    return TCL_OK;
  }

  // Calculate bending energy.
  double energy = 0;
  if ( !asiAlgo_Utils::CalculateBendingEnergy(occtSurf, energy) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Cannot calculate bending energy.");
    return TCL_OK;
  }

  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Approximate bending energy: %1."
                                                       << energy);

  // Add to interpreter.
  *interp << energy;

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckEdgeVexity(const Handle(asiTcl_Interp)& interp,
                           int                          argc,
                           const char**                 argv)
{
  if ( argc != 2 && argc != 3 && argc != 4 && argc != 6 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get Part Node and shape.
  Handle(asiData_PartNode) part_n    = cmdEngine::model->GetPartNode();
  TopoDS_Shape             partShape = part_n->GetShape();

  // Check whether naming service is active.
  const bool hasNaming = !part_n->GetNaming().IsNull();

  // Check whether "smooth" type is allowed.
  const bool allowSmooth = interp->HasKeyword(argc, argv, "allow-smooth");

  // Get tolerance for recognition of smooth edges.
  double smoothToler = 0.0;
  TCollection_AsciiString smoothTolerStr;
  //
  if ( interp->GetKeyValue(argc, argv, "smooth-toler", smoothTolerStr) )
    smoothToler = smoothTolerStr.RealValue();

  // Edge to check.
  TopoDS_Edge edge;

  // Check if naming service is active. If so, the user may ask to access
  // a sub-shape in question by its unique name.
  if ( hasNaming && argc == 3 )
  {
    if ( !interp->IsKeyword(argv[1], "name") )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Keyword '-name' is expected.");
      return TCL_ERROR;
    }
    else
    {
      TCollection_AsciiString name(argv[2]);
      //
      TopoDS_Shape subshape = part_n->GetNaming()->GetShape(name);
      //
      if ( subshape.IsNull() || subshape.ShapeType() != TopAbs_EDGE )
      {
        interp->GetProgress().SendLogMessage(LogErr(Normal) << "The passed sub-shape is null "
                                                               "or not of a proper type.");
        return TCL_ERROR;
      }
      //
      edge = TopoDS::Edge(subshape);
    }
  }
  else // Naming is not used to access the argument.
  {
    const int ssidx = atoi(argv[1]);
    //
    if ( ssidx < 1 )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Sub-shape index should be 1-based.");
      return TCL_ERROR;
    }

    // Get map of sub-shapes with respect to those the passed index is relevant.
    TopTools_IndexedMapOfShape subShapesOfType;
    part_n->GetAAG()->RequestMapOf(TopAbs_EDGE, subShapesOfType);

    // Get sub-shape in question.
    edge = TopoDS::Edge( subShapesOfType(ssidx) );
  }

  interp->GetPlotter().REDRAW_SHAPE("edge", edge, Color_White, 1.0, true);

  // Get faces owning the edge in question.
  TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
  TopExp::MapShapesAndAncestors(part_n->GetShape(), TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);
  //
  if ( !edgeFaceMap.Contains(edge) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "No faces for edge.");
    return TCL_ERROR;
  }
  //
  const TopTools_ListOfShape& faces = edgeFaceMap.FindFromKey(edge);
  //
  if ( faces.Extent() != 2 )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Edge should belong to exactly two faces.");
    return TCL_ERROR;
  }
  //
  const TopoDS_Face& F = TopoDS::Face( faces.First() );
  const TopoDS_Face& G = TopoDS::Face( faces.Last() );

  /* ===============
   *  Perform check
   * =============== */

  TIMER_NEW
  TIMER_GO

  // Prepare algorithm.
  asiAlgo_CheckDihedralAngle dihAngle( interp->GetProgress(),
                                       interp->GetPlotter() );
  //
  dihAngle.SetCommonEdge(edge);

  // Check.
  const asiAlgo_FeatureAngleType
    angleType = dihAngle.AngleBetweenFaces(F, G, allowSmooth, smoothToler);

  TIMER_FINISH
  TIMER_COUT_RESULT_NOTIFIER(interp->GetProgress(), "Check dihedral angle")

  // Add as a result to interpreter.
  *interp << angleType;

  interp->GetProgress().SendLogMessage( LogInfo(Normal) << "Edge vexity: %1."
                                                        << asiAlgo_Utils::FeaturAngleToString(angleType) );

  // Update UI.
  if ( cmdEngine::cf && cmdEngine::cf->ViewerPart )
  {
    double colorR, colorG, colorB;
    asiVisu_Utils::ColorForFeatureAngle(angleType, colorR, colorG, colorB);
    //
    interp->GetPlotter().REDRAW_SHAPE("vexity",
                                       edge,
                                       ActAPI_Color(colorR,
                                                    colorG,
                                                    colorB,
                                                    Quantity_TOC_RGB),
                                       1.0,
                                       true);
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_DrawCP(const Handle(asiTcl_Interp)& interp,
                  int                          argc,
                  const char**                 argv)
{
  if ( argc != 4 && argc != 5 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Find Surface Node by name.
  Handle(ActAPI_INode) node = cmdEngine::model->FindNodeByName(argv[1]);
  //
  if ( node.IsNull() || !node->IsKind( STANDARD_TYPE(asiData_IVSurfaceNode) ) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Node '%1' is not a surface."
                                                        << argv[1]);
    return TCL_OK;
  }
  //
  Handle(asiData_IVSurfaceNode)
    surfNode = Handle(asiData_IVSurfaceNode)::DownCast(node);

  // Get parametric surface.
  Handle(Geom_BSplineSurface)
    occtSurf = Handle(Geom_BSplineSurface)::DownCast( surfNode->GetSurface() );
  //
  if ( occtSurf.IsNull() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "The surface is NULL or not a spline surface.");
    return TCL_ERROR;
  }

  // Get indices of a control point in question.
  const int i = atoi(argv[2]);
  const int j = atoi(argv[3]);

  const TColgp_Array2OfPnt& poles = occtSurf->Poles();
  //
  const int numRows = poles.ColLength();
  const int numCols = poles.RowLength();
  //
  if ( i < 0 || i >= numRows )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "First index if out of range [0, %1]." << (numRows - 1));
    return TCL_ERROR;
  }
  //
  if ( j < 0 || j >= numCols )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Second index if out of range [0, %1]." << (numCols - 1));
    return TCL_ERROR;
  }

  // Check whether Mobius operator is requested.
  gp_Pnt P;
  //
  const bool isMobius = interp->HasKeyword(argc, argv, "mobius");
  //
  if ( isMobius )
  {
#if defined USE_MOBIUS
    // Convert to Mobius surface.
    mobius::cascade_BSplineSurface converter(occtSurf);
    converter.DirectConvert();
    //
    const mobius::ptr<mobius::bsurf>& mobSurf = converter.GetMobiusSurface();

    const mobius::xyz& mobP = mobSurf->GetPoles()[i][j];
    //
    P.SetX( mobP.X() );
    P.SetY( mobP.Y() );
    P.SetZ( mobP.Z() );
#else
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Mobius is not available.");
    return TCL_ERROR;
#endif
  }
  else
  {
    P = poles(i + 1, j + 1);
  }

  // Draw control point.
  TCollection_AsciiString PName("P_");
  PName += i;
  PName += "_";
  PName += j;
  //
  interp->GetPlotter().REDRAW_POINT(PName, P, isMobius ? Color_Blue : Color_Red);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckTriangulation(const Handle(asiTcl_Interp)& interp,
                              int                          argc,
                              const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get triangulation.
  Handle(asiData_TriangulationNode)
    trisNode = cmdEngine::model->GetTriangulationNode();
  //
  if ( trisNode.IsNull() || !trisNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Triangulation Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  Handle(Poly_Triangulation) tris = trisNode->GetTriangulation();

  // Perform check.
  asiAlgo_MeshCheckInter meshChecker( tris,
                                      interp->GetProgress(),
                                      interp->GetPlotter() );
  //
  asiAlgo_MeshCheckInter::t_status res = meshChecker.Perform();

  if ( res == asiAlgo_MeshCheckInter::Status_Ok )
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Triangulation is OK.");
  else if ( res == asiAlgo_MeshCheckInter::Status_HasIntersections )
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Triangulation is self-intersecting.");
  else if ( res == asiAlgo_MeshCheckInter::Status_Failed )
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Triangulation checker failed for unknown reason.");

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckTess(const Handle(asiTcl_Interp)& interp,
                     int                          argc,
                     const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get triangulation.
  Handle(asiData_TessNode)
    tessNode = cmdEngine::model->GetTessellationNode();
  //
  if ( tessNode.IsNull() || !tessNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Tessellation Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  Handle(Poly_Triangulation) tris;
  asiAlgo_MeshConvert::FromPersistent(tessNode->GetMesh(), tris);

  // Perform check.
  asiAlgo_MeshCheckInter meshChecker( tris,
                                      interp->GetProgress(),
                                      interp->GetPlotter() );
  //
  asiAlgo_MeshCheckInter::t_status res = meshChecker.Perform();

  if ( res == asiAlgo_MeshCheckInter::Status_Ok )
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Triangulation is OK.");
  else if ( res == asiAlgo_MeshCheckInter::Status_HasIntersections )
    interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Triangulation is self-intersecting.");
  else if ( res == asiAlgo_MeshCheckInter::Status_Failed )
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Triangulation checker failed for unknown reason.");

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckOpenEdges(const Handle(asiTcl_Interp)& interp,
                          int                          argc,
                          const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get part.
  Handle(asiData_PartNode)
    partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape partShape = partNode->GetShape();

  // Build map of edges to extract open ("naked") ones.
  TopTools_IndexedDataMapOfShapeListOfShape edgesFaces;
  TopExp::MapShapesAndAncestors(partShape, TopAbs_EDGE, TopAbs_FACE, edgesFaces);

  // Find open edges.
  Handle(TopTools_HSequenceOfShape) openEdgesSeq = new TopTools_HSequenceOfShape;
  //
  for ( int k = 1; k <= edgesFaces.Extent(); ++k )
  {
    const TopTools_ListOfShape& faces = edgesFaces(k);
    //
    if ( faces.Extent() == 1 )
    {
      const TopoDS_Edge& E = TopoDS::Edge( edgesFaces.FindKey(k) );
      //
      if ( BRep_Tool::Degenerated(E) )
        continue;

      openEdgesSeq->Append(E);
    }
  }

  // Compose border wires from the naked edges.
  Handle(TopTools_HSequenceOfShape) borderWires;
  ShapeAnalysis_FreeBounds::ConnectEdgesToWires(openEdgesSeq, 1e-3, 0, borderWires);

  // Dump free wires.
  if ( !borderWires.IsNull() )
  {
    for ( TopTools_HSequenceOfShape::Iterator wit(*borderWires); wit.More(); wit.Next() )
    {
      const TopoDS_Wire& wire = TopoDS::Wire( wit.Value() );

      // Calculate global properties.
      GProp_GProps props;
      BRepGProp::LinearProperties(wire, props);
      const double wireLen = props.Mass();

      interp->GetProgress().SendLogMessage(LogInfo(Normal) << "Length of the next border wire: %1." << wireLen);
      interp->GetPlotter().DRAW_SHAPE(wire, Color_Red, 1.0, true, "free_wire");
    }
  }

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckVerticesOri(const Handle(asiTcl_Interp)& interp,
                            int                          argc,
                            const char**                 argv)
{
  if ( argc != 1 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get part.
  Handle(asiData_PartNode)
    partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape partShape = partNode->GetShape();

  // Build map of edges to check them one-by-one.
  TopTools_IndexedMapOfShape edgesMap;
  TopExp::MapShapes(partShape, TopAbs_EDGE, edgesMap);

  // Check edges.
  for ( int k = 1; k <= edgesMap.Extent(); ++k )
  {
    const TopoDS_Edge& edge = TopoDS::Edge( edgesMap(k) );

    if ( !asiAlgo_CheckValidity::HasDistinctVertexOrientations(edge) )
    {
      interp->GetProgress().SendLogMessage(LogErr(Normal) << "Edge %1 has non-distinguishable vertices."
                                                          << k);
      //
      *interp << false;

      return TCL_OK;
    }
  }

  // Everything is fine.
  interp->GetProgress().SendLogMessage(LogInfo(Normal) << "All edges have distinguishable vertices.");
  //
  *interp << true;

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_CheckInternalLocations(const Handle(asiTcl_Interp)& interp,
                                  int                          argc,
                                  const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  const bool skipFirstLevel = interp->HasKeyword(argc, argv, "skip-first-level");

  // Get part.
  Handle(asiData_PartNode)
    partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape partShape = partNode->GetShape();

  // Check and result.
  *interp << asiAlgo_Utils::HasInternalLocations(partShape, skipFirstLevel);

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_RecognizeBlends(const Handle(asiTcl_Interp)& interp,
                           int                          argc,
                           const char**                 argv)
{
  if ( argc != 1 && argc != 2 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  // Get part.
  Handle(asiData_PartNode)
    partNode = cmdEngine::model->GetPartNode();
  //
  if ( partNode.IsNull() || !partNode->IsWellFormed() )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Part Node is null or ill-defined.");
    return TCL_OK;
  }
  //
  TopoDS_Shape        partShape = partNode->GetShape();
  Handle(asiAlgo_AAG) partAAG   = partNode->GetAAG();

  // Get max radius.
  const double maxRadius = ( (argc == 1) ? 1e100 : atof(argv[1]) );

  // Perform recognition.
  asiAlgo_RecognizeBlends recognizer( partShape,
                                      partAAG,
                                      interp->GetProgress(),
                                      NULL );
  //
  if ( !recognizer.Perform(maxRadius) )
  {
    interp->GetProgress().SendLogMessage(LogErr(Normal) << "Recognition failed.");
    return TCL_ERROR;
  }

  // Highlight the detected faces.
  const TColStd_PackedMapOfInteger& resIndices = recognizer.GetResultIndices();
  //
  if ( !cmdEngine::cf.IsNull() && cmdEngine::cf->ViewerPart )
    asiEngine_Part( cmdEngine::model,
                    cmdEngine::cf->ViewerPart->PrsMgr() ).HighlightFaces(resIndices);

  // Dump to result.
  *interp << resIndices;

  return TCL_OK;
}

//-----------------------------------------------------------------------------

int ENGINE_DrawPlot(const Handle(asiTcl_Interp)& interp,
                    int                          argc,
                    const char**                 argv)
{
  if ( argc < 5 )
  {
    return interp->ErrorOnWrongArgs(argv[0]);
  }

  const bool isLogScale = interp->HasKeyword(argc, argv, "log");

  // Get table of values to plot.
  std::vector<double> xvec, fxvec;
  //
  for ( int k = 1; k < argc; k += 2 )
  {
    if ( interp->IsKeyword(argv[k], "log") ) continue;

    const double x  = atof(argv[k]);
    const double fx = atof(argv[k + 1]);

    xvec.push_back(x);
    fxvec.push_back(fx);
  }

  // Open curvature plot.
  asiUI_Plot2d*
    cPlot = new asiUI_Plot2d( interp->GetProgress(),
                              interp->GetPlotter() );
  //
  cPlot->SetLogScale(isLogScale);
  //
  cPlot->Render(xvec, fxvec, "x", "fx", "Plot 2d");

  return TCL_OK;
}

//-----------------------------------------------------------------------------

void cmdEngine::Commands_Inspection(const Handle(asiTcl_Interp)&      interp,
                                    const Handle(Standard_Transient)& data)
{
  cmdEngine_NotUsed(data);
  //
  static const char* group = "cmdEngine";

  //-------------------------------------------------------------------------//
  interp->AddCommand("explode",
    //
    "explode [-vertex|-edge|-wire|-face|-shell|-solid|-noloc]\n"
    "\t Explodes active part to sub-shapes of interest. If no sub-shape\n"
    "\t qualifier is passed, this command explodes the part to its direct\n"
    "\t children (e.g. edges for wire, wires for face, etc.). If '-noloc'\n"
    "\t flag is passed, locations are not accumulated on explode.",
    //
    __FILE__, group, ENGINE_Explode);

  //-------------------------------------------------------------------------//
  interp->AddCommand("print-summary",
    //
    "print-summary\n"
    "\t Prints summary (number of sub-shapes) for the active part.",
    //
    __FILE__, group, ENGINE_PrintSummary);

  //-------------------------------------------------------------------------//
  interp->AddCommand("get-summary",
    //
    "get-summary <compounds> <compsolids> <solids> <shells> <faces> <wires> <edges> <vertices>\n"
    "\t Returns summary (number of sub-shapes) to the specified output variables.",
    //
    __FILE__, group, ENGINE_GetSummary);

  //-------------------------------------------------------------------------//
  interp->AddCommand("face-addr",
    //
    "face-addr faceIndex\n"
    "\t Prints physical address of the given face.",
    //
    __FILE__, group, ENGINE_FaceAddr);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-dist",
    //
    "check-dist varName\n"
    "\t Computes distance between the part and the given topological object.",
    //
    __FILE__, group, ENGINE_CheckDist);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-curvature",
    //
    "check-curvature [numPts [scaleFactor [curvAmpl]]] [-noplot] [-noalong]\n"
    "\t Checks curvature of the selected edge. As a result, curvature combs\n"
    "\t are visualized in 3D. You can control its scale factor with\n"
    "\t <scaleFactor> argument and also its density with <numPts> argument.\n"
    "\t To bring out the salient features of the comb, <curvAmpl> amplification\n"
    "\t factor can be used. If -noplot key is passed, the curvature plot is not\n"
    "\t constructed. If -noalong key is passed, the along-curvature value for\n"
    "\t the selected edges is not computed.",
    //
    __FILE__, group, ENGINE_CheckCurvature);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-continuity",
    //
    "check-continuity\n"
    "\t Checks continuity of the selected face.",
    //
    __FILE__, group, ENGINE_CheckContinuity);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-euler",
    //
    "check-euler [genus]\n"
    "\t Opens dialog to check Euler-Poincare property of the Part geometry.\n"
    "\t If <genus> parameter is not specified, this command will open a prompt\n"
    "\t dialog to ask the user to type genus.",
    //
    __FILE__, group, ENGINE_CheckEuler);

  //-------------------------------------------------------------------------//
  interp->AddCommand("dump-topograph-dot",
    //
    "dump-topograph-dot filename\n"
    "\t Dumps topology graph to DOT file (can be opened by Graphviz).",
    //
    __FILE__, group, ENGINE_DumpTopographDot);

  //-------------------------------------------------------------------------//
  interp->AddCommand("eval-curve",
    //
    "eval-curve curveName u order [-mobius]\n"
    "\t Evaluates curve <curveName> for the given parameter value <u>.\n"
    "\t If <-mobius> keyword is used, evaluation is performed using Mobius\n"
    "\t functions. The argument <order> specifies the order of derivatives\n"
    "\t to evaluate (0 for value evaluation).",
    //
    __FILE__, group, ENGINE_EvalCurve);

  //-------------------------------------------------------------------------//
  interp->AddCommand("eval-surf",
    //
    "eval-surf surfName u v order [-mobius]\n"
    "\t Evaluates surface <surfName> for the given parameter pair <u, v>.\n"
    "\t If <-mobius> keyword is used, evaluation is performed using Mobius\n"
    "\t functions. The argument <order> specifies the order of derivatives\n"
    "\t to evaluate (0 for value evaluation).",
    //
    __FILE__, group, ENGINE_EvalSurf);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-toler",
    //
    "check-toler [numRanges]\n"
    "\t Checks local tolerances in the part shape and distributes sub-shapes\n."
    "\t by tolerance ranges. You can control the fineness of tolerance ranges\n"
    "\t using <numRanges> optional argument.",
    //
    __FILE__, group, ENGINE_CheckToler);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-length",
    //
    "check-length\n"
    "\t Checks length of the selected edges.",
    //
    __FILE__, group, ENGINE_CheckLength);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-area",
    //
    "check-area\n"
    "\t Checks area of the selected faces.",
    //
    __FILE__, group, ENGINE_CheckArea);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-validity",
    //
    "check-validity\n"
    "\t Checks validity of the part shape.",
    //
    __FILE__, group, ENGINE_CheckValidity);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-finite",
    //
    "check-finite\n"
    "\t Checks finiteness of the part shape.",
    //
    __FILE__, group, ENGINE_CheckFinite);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-contours",
    //
    "check-contours [tolerance]\n"
    "\t Checks whether all contours of the part shape are closed or not.",
    //
    __FILE__, group, ENGINE_CheckContours);

  //-------------------------------------------------------------------------//
  interp->AddCommand("get-tolerance",
    //
    "get-tolerance\n"
    "\t Returns max geometric tolerance of the part shape.",
    //
    __FILE__, group, ENGINE_GetTolerance);

  //-------------------------------------------------------------------------//
  interp->AddCommand("get-curve-strain",
    //
    "get-curve-strain curveName\n"
    "\t Returns strain energy of the passed curve.",
    //
    __FILE__, group, ENGINE_GetCurveStrain);

  //-------------------------------------------------------------------------//
  interp->AddCommand("get-surface-bending",
    //
    "get-surface-bending surfName\n"
    "\t Returns bending energy of the passed surface.",
    //
    __FILE__, group, ENGINE_GetSurfaceBending);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-edge-vexity",
    //
    "check-edge-vexity <edgeIndex|-name 'edgeName'>\n"
    "\t Rebuilds edge with the given ID or name.",
    //
    __FILE__, group, ENGINE_CheckEdgeVexity);

  //-------------------------------------------------------------------------//
  interp->AddCommand("draw-cp",
    //
    "draw-cp surfName i_index j_index [-mobius]\n"
    "\t Draws control point of a free-form surface.",
    //
    __FILE__, group, ENGINE_DrawCP);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-triangulation",
    //
    "check-triangulation\n"
    "\t Checks triangulation for self-intersections.",
    //
    __FILE__, group, ENGINE_CheckTriangulation);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-tess",
    //
    "check-tess\n"
    "\t Checks tessellation for self-intersections.",
    //
    __FILE__, group, ENGINE_CheckTess);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-open-edges",
    //
    "check-open-edges\n"
    "\t Checks part for open edges.",
    //
    __FILE__, group, ENGINE_CheckOpenEdges);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-vertices-ori",
    //
    "check-vertices-ori\n"
    "\t Checks if the vertices of the part are distinguishable by their orientation\n"
    "\t flags. A valid edge should have one FORWARD vertex and one REVERSED\n"
    "\t vertex to denote its topological extremities. As a result of some\n"
    "\t improperly implemented modeling operators, this rule can be broken.",
    //
    __FILE__, group, ENGINE_CheckVerticesOri);

  //-------------------------------------------------------------------------//
  interp->AddCommand("check-internal-locations",
    //
    "check-internal-locations [-skip-first-level]\n"
    "\t Checks if the part contains any internal non-identity locations.",
    //
    __FILE__, group, ENGINE_CheckInternalLocations);

  //-------------------------------------------------------------------------//
  interp->AddCommand("recognize-blends",
    //
    "recognize-blends radius\n"
    "\t Recognizes all blend faces in AAG representing the part.",
    //
    __FILE__, group, ENGINE_RecognizeBlends);

  //-------------------------------------------------------------------------//
  interp->AddCommand("draw-plot",
    //
    "draw-plot x1 fx1 x2 fx2 [x3 fx3 [...]] [-log]\n"
    "\t Draws two-dimensional plot of the given values.",
    //
    __FILE__, group, ENGINE_DrawPlot);
}
