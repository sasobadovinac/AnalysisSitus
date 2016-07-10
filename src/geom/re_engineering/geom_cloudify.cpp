//-----------------------------------------------------------------------------
// Created on: 02 June 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <geom_cloudify.h>

// Geometry includes
#include <geom_classify_point_face.h>

// Qr includes
#include <QrGeom3D_PositionCloud.h>

// OCCT includes
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <gce_MakeCirc.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

static Handle(geom_point_cloud) repack(const QrPtr<pcloud>& qrCloud)
{
  // Re-pack coordinates to the conventional persistent form
  const std::vector<xyz>& coords = qrCloud->Points();
  Handle(TColStd_HArray1OfReal) packed_coords = new TColStd_HArray1OfReal(0, (int) qrCloud->NumberOfPoints()*3 - 1);
  //
  int idx = 0;
  for ( size_t p = 0; p < coords.size(); ++p )
  {
    packed_coords->ChangeValue(idx++) = coords[p].X();
    packed_coords->ChangeValue(idx++) = coords[p].Y();
    packed_coords->ChangeValue(idx++) = coords[p].Z();
  }

  // Set result
  Handle(geom_point_cloud) point_cloud = new geom_point_cloud(packed_coords);
  return point_cloud;
}

//-----------------------------------------------------------------------------

//! Constructor.
//! \param uv_step  [in] sampling step.
//! \param progress [in] Progress Notifier.
//! \param plotter  [in] Imperative Plotter.
geom_cloudify::geom_cloudify(const double         uv_step,
                             ActAPI_ProgressEntry progress,
                             ActAPI_PlotterEntry  plotter)
: ActAPI_IAlgorithm (progress, plotter),
  m_fLinStep        (uv_step)
{}

//-----------------------------------------------------------------------------

//! Builds a representative point cloud for the surfaces of the given model.
//! \param model       [in]  target CAD model.
//! \param point_cloud [out] result point cloud.
//! \return true in case of success, false -- otherwise.
bool geom_cloudify::Sample_Faces(const TopoDS_Shape&       model,
                                 Handle(geom_point_cloud)& point_cloud)
{
  QrPtr<pcloud> qrCloud = new pcloud;

  { // Progress [begin]
    TopTools_IndexedMapOfShape M;
    TopExp::MapShapes(model, TopAbs_FACE, M);
    //
    this->Progress().Init( M.Extent() );
  } // Progress [end]

  // Iterate over the faces, and sample points in their domains
  for ( TopExp_Explorer exp(model, TopAbs_FACE); exp.More(); exp.Next() )
  {
    const TopoDS_Face& face = TopoDS::Face( exp.Current() );

    // Surface adaptor
    BRepAdaptor_Surface bas(face);

    // Get parametric bounds
    double uMin, uMax, vMin, vMax;
    BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);

    // Choose an adequate sampling step in parametric space for the
    // current face. A user is not aware of any parametric spaces as a rule,
    // so we'd better treat his sampling argument as a linear metric
    // property in three dimensions. Obviously, some distortion
    // coefficient has to be applied to stick to this metric constraint
    // in 3D when we are sampling 2D

    // We do not want to use curvature since parameterization of a host
    // surface may happen to be irregular. So we are using D0-wise
    // heuristic, like we pick up two points, take the one in-between,
    // and pass a circle using three point. The radius of that circle
    // we take as a curvature radius. Having this radius R, it is easy
    // to derive an angle giving us a certain arc length:
    // L = alpha * R => alpha = L / R

    const double uStep = this->chooseParametricStep(bas, true, uMin, uMax, vMin, vMax);
    const double vStep = this->chooseParametricStep(bas, false, uMin, uMax, vMin, vMax);

    // Prepare classifier
    geom_classify_point_face classifier(face, BRep_Tool::Tolerance(face), 0.01);

    // Sample points
    double u = uMin;
    bool uStop = false;
    while ( !uStop )
    {
      if ( u > uMax )
      {
        u     = uMax;
        uStop = true;
      }

       double v = vMin;
       bool vStop = false;
       while ( !vStop )
       {
         if ( v > vMax )
         {
           v     = vMax;
           vStop = true;
         }

         // Perform point membership classification
         geom_membership pmc = classifier(gp_Pnt2d(u, v), NULL);
         //
         if ( pmc & Membership_InOn )
         {
           gp_Pnt P = bas.Value(u, v);
           qrCloud->AddPoint( xyz( P.X(), P.Y(), P.Z() ) );
         }

         v += vStep;
       }

      u += uStep;
    }

    { // Progress [begin]
      this->Progress().StepProgress(1, 1);
    } // Progress [end]
  }

  // Convert to a persistent point cloud
  point_cloud = repack(qrCloud);
  return true;
}

//-----------------------------------------------------------------------------

//! Builds a representative point cloud for the facets of the given model.
//! \param model       [in]  target CAD model to take facets from.
//! \param point_cloud [out] result point cloud.
//! \return true in case of success, false -- otherwise.
bool geom_cloudify::Sample_Facets(const TopoDS_Shape&       model,
                                  Handle(geom_point_cloud)& point_cloud)
{
  QrPtr<pcloud> qrCloud = new pcloud;

  // Constants
  const double lower = 0.0, upper = 1.0;

  // Loop over the facets
  for ( TopExp_Explorer exp(model, TopAbs_FACE); exp.More(); exp.Next() )
  {
    const TopoDS_Face& F = TopoDS::Face( exp.Current() );

    // Ask for the facets belonging to the given face
    TopLoc_Location L;
    const Handle(Poly_Triangulation)& T = BRep_Tool::Triangulation(F, L);
    //
    if ( T.IsNull() )
      continue;

    // Take data arrays
    const TColgp_Array1OfPnt&   nodes = T->Nodes();
    const Poly_Array1OfTriangle& tris = T->Triangles();

    // Loop over the array of triangles, so that we can work with each
    // individual facet independently
    for ( int t = 1; t <= tris.Length(); ++t )
    {
      const Poly_Triangle& tri = tris(t);

      // Take triangle's nodes
      int n[3];
      tri.Get(n[0], n[1], n[2]);

      // Check out the nodes
      gp_XYZ r[3] = { nodes(n[0]).XYZ(), nodes(n[1]).XYZ(), nodes(n[2]).XYZ() };

      // There are two parameters for sampling. One runs from r[1] to r[2],
      // while the second runs from r[0] to the intermediate point between
      // r[1] and r[2] (defined by the first parameter). The first parameter
      // is alpha, the second is beta
      double alpha   = lower;
      bool alphaStop = false;
      //
      while ( !alphaStop )
      {
        if ( alpha > upper )
        {
          alpha     = upper;
          alphaStop = true;
        }

        double beta   = lower;
        bool betaStop = false;
        //
        while ( !betaStop )
        {
          if ( beta > upper )
          {
            beta     = upper;
            betaStop = true;
          }

          gp_XYZ P = r[0] + beta*(r[1] + alpha*(r[2] - r[1]) - r[0]);
          //
          qrCloud->AddPoint( xyz( P.X(), P.Y(), P.Z() ) );

          beta += m_fLinStep;
        }

        alpha += m_fLinStep;
      }
    }
  }

  // Convert to a persistent point cloud
  point_cloud = repack(qrCloud);
  return true;
}

//-----------------------------------------------------------------------------

//! This method allows to choose a parametric step in adaptive way, basing on
//! the linear step in 3D and the curvature radius of the host surface.
//! \param bas  [in] host surface.
//! \param isU  [in] true for U, false for V.
//! \param uMin [in] minimal U.
//! \param uMax [in] maximal U.
//! \param vMin [in] minimal V.
//! \param vMax [in] maximal V.
//! \return calculated parametric step.
double geom_cloudify::chooseParametricStep(const BRepAdaptor_Surface& bas,
                                           const bool                 isU,
                                           const double               uMin,
                                           const double               uMax,
                                           const double               vMin,
                                           const double               vMax) const
{
  const double tMin = isU ? uMin : vMin;
  const double tMax = isU ? uMax : vMax;

  double step;
  //
  {
    const double dt = (tMax - tMin)*0.01;
    const double t0 = tMin;
    const double t1 = tMin + dt;
    const double tm = tMin + dt*0.5;
    //
    const gp_Pnt P0 = isU ? bas.Value(t0, vMin) : bas.Value(uMin, t0);
    const gp_Pnt P1 = isU ? bas.Value(t1, vMin) : bas.Value(uMin, t1);
    const gp_Pnt Pm = isU ? bas.Value(tm, vMin) : bas.Value(uMin, tm);
    //
    bool isOk;
    double R = 0.0;
    try
    {
      gce_MakeCirc mkCirc(P0, P1, Pm);
      isOk = mkCirc.IsDone() > 0;
      R = mkCirc.Value().Radius();
    }
    catch ( ... )
    {
      isOk = false;
    }
    //
    if ( isOk )
    {
      // If we're not done, it means that our points are collinear, so
      // there is no sense to apply heuristic here; just use what the user
      // has given us
      step = ( R > RealEpsilon() ) ? (m_fLinStep / R) : m_fLinStep;
    }
    else
    {
      step = m_fLinStep;
    }
  }

  return step;
}