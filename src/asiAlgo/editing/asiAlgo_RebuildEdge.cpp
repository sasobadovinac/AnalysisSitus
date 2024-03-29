//-----------------------------------------------------------------------------
// Created on: 14 May (*) 2018
//-----------------------------------------------------------------------------
// Copyright (c) 2018, Sergey Slyadnev
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

// Own include
#include <asiAlgo_RebuildEdge.h>

// asiAlgo includes
#include <asiAlgo_BRepNormalizer.h>
#include <asiAlgo_ModConstructEdge.h>

//-----------------------------------------------------------------------------

asiAlgo_RebuildEdge::asiAlgo_RebuildEdge(const TopoDS_Shape&  masterCAD,
                                         ActAPI_ProgressEntry progress,
                                         ActAPI_PlotterEntry  plotter)
//
: ActAPI_IAlgorithm (progress, plotter),
  m_input           (masterCAD)
{}

//-----------------------------------------------------------------------------

bool asiAlgo_RebuildEdge::Perform(const TopoDS_Edge& edge)
{
  /* ===================
   *  Preparation stage
   * =================== */

  // Prepare history.
  if ( m_history.IsNull() )
    m_history = new asiAlgo_History;

  /* ==========================
   *  Apply geometric operator
   * ========================== */

  // Prepare Modification.
  Handle(asiAlgo_ModConstructEdge)
    Mod = new asiAlgo_ModConstructEdge(m_input,
                                       m_progress,
                                       m_plotter);

  // Initialize Modification.
  if ( !Mod->Init(edge) )
    return false;
  //
  Mod->SetFrozenVertices(m_frozenVertices.vertices);
  //
  Mod->Dump(m_plotter);

  // Initialize Modifier.
  asiAlgo_BRepNormalizer Modifier(m_progress, m_plotter);
  Modifier.Init(m_input);

  // Perform Modification.
  if ( !Modifier.Perform(Mod) )
  {
    m_progress.SendLogMessage(LogErr(Normal) << "Geometry normalization failed.");
    return false;
  }

  /* ==========
   *  Finalize
   * ========== */

  // Get result.
  if ( !Modifier.ModifiedShape(m_input, m_result) )
  {
    m_progress.SendLogMessage(LogErr(Normal) << "There is no image for the input shape.");
    return false;
  }

  // Populate history.
  {
    const asiAlgo_ModEdgeInfo& edgeInfo = Mod->GetEdgeInfo();

    /* e_s1_s2 */

    TopoDS_Shape e_s1_s2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.e_s1_s2, e_s1_s2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for e_s1_s2.");
      return false;
    }

    /* e_s1_t1 */

    TopoDS_Shape e_s1_t1_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.e_s1_t1, e_s1_t1_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for e_s1_t1.");
      return false;
    }

    /* e_s1_t2 */

    TopoDS_Shape e_s1_t2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.e_s1_t2, e_s1_t2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for e_s1_t2.");
      return false;
    }

    /* e_s2_t1 */

    TopoDS_Shape e_s2_t1_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.e_s2_t1, e_s2_t1_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for e_s2_t1.");
      return false;
    }

    /* e_s2_t2 */

    TopoDS_Shape e_s2_t2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.e_s2_t2, e_s2_t2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for e_s2_t2.");
      return false;
    }

    /* f_s1 */

    TopoDS_Shape f_s1_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.f_s1, f_s1_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for f_s1.");
      return false;
    }

    /* f_s2 */

    TopoDS_Shape f_s2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.f_s2, f_s2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for f_s2.");
      return false;
    }

    /* f_t1 */

    TopoDS_Shape f_t1_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.f_t1, f_t1_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for f_t1.");
      return false;
    }

    /* f_t2 */

    TopoDS_Shape f_t2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.f_t2, f_t2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for f_t2.");
      return false;
    }

    /* v_s1_s2_t1 */

    TopoDS_Shape v_s1_s2_t1_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.v_s1_s2_t1, v_s1_s2_t1_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for v_s1_s2_t1.");
      return false;
    }

    /* v_s1_s2_t2 */

    TopoDS_Shape v_s1_s2_t2_new;
    if ( !Modifier.ModifiedShape(edgeInfo.situation.v_s1_s2_t2, v_s1_s2_t2_new) )
    {
      m_progress.SendLogMessage(LogErr(Normal) << "There is no image for v_s1_s2_t2.");
      return false;
    }

    /* Extra edges */

    std::vector< std::pair<TopoDS_Shape, TopoDS_Shape> > eExtra_images;
    for ( int k = 1; k <= edgeInfo.situation.e_extra.Extent(); ++k )
    {
      TopoDS_Shape eExtra_old = edgeInfo.situation.e_extra(k);
      TopoDS_Shape eExtra_new;
      //
      if ( !Modifier.ModifiedShape(eExtra_old, eExtra_new) )
      {
        m_progress.SendLogMessage(LogErr(Normal) << "There is no image for extra edge %1." << k);
        return false;
      }

      // Add to the list of images.
      eExtra_images.push_back( std::pair<TopoDS_Shape, TopoDS_Shape>(eExtra_old, eExtra_new) );
    }

    /* Extra faces */

    std::vector< std::pair<TopoDS_Shape, TopoDS_Shape> > fExtra_images;
    for ( int k = 1; k <= edgeInfo.situation.f_extra.Extent(); ++k )
    {
      TopoDS_Shape fExtra_old = edgeInfo.situation.f_extra(k);
      TopoDS_Shape fExtra_new;
      //
      if ( !Modifier.ModifiedShape(fExtra_old, fExtra_new) )
      {
        m_progress.SendLogMessage(LogErr(Normal) << "There is no image for extra face %1." << k);
        return false;
      }

      // Add to the list of images.
      fExtra_images.push_back( std::pair<TopoDS_Shape, TopoDS_Shape>(fExtra_old, fExtra_new) );
    }

    // Populate history.
    m_history->AddModified(edgeInfo.situation.e_s1_s2,    e_s1_s2_new);
    m_history->AddModified(edgeInfo.situation.e_s1_t1,    e_s1_t1_new);
    m_history->AddModified(edgeInfo.situation.e_s1_t2,    e_s1_t2_new);
    m_history->AddModified(edgeInfo.situation.e_s2_t1,    e_s2_t1_new);
    m_history->AddModified(edgeInfo.situation.e_s2_t2,    e_s2_t2_new);
    m_history->AddModified(edgeInfo.situation.f_s1,       f_s1_new);
    m_history->AddModified(edgeInfo.situation.f_s2,       f_s2_new);
    m_history->AddModified(edgeInfo.situation.f_t1,       f_t1_new);
    m_history->AddModified(edgeInfo.situation.f_t2,       f_t2_new);
    m_history->AddModified(edgeInfo.situation.v_s1_s2_t1, v_s1_s2_t1_new);
    m_history->AddModified(edgeInfo.situation.v_s1_s2_t2, v_s1_s2_t2_new);
    //
    for ( size_t k = 0; k < eExtra_images.size(); ++k )
      m_history->AddModified(eExtra_images[k].first, eExtra_images[k].second);
    //
    for ( size_t k = 0; k < fExtra_images.size(); ++k )
      m_history->AddModified(fExtra_images[k].first, fExtra_images[k].second);
  }

  return true; // Success.
}

//-----------------------------------------------------------------------------

void asiAlgo_RebuildEdge::AddFrozenVertex(const TopoDS_Vertex& vertex)
{
  m_frozenVertices.Add(vertex);
}

//-----------------------------------------------------------------------------

void asiAlgo_RebuildEdge::SetFrozenVertices(const asiAlgo_FrozenVertices& vertices)
{
  m_frozenVertices = vertices;
}
