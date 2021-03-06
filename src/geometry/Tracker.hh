//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Tracker.hh
 * \brief  Tracker class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef TRACKER_HH_
#define TRACKER_HH_

// Detran
#include "MeshMOC.hh"
#include "QuadratureMOC.hh"
#include "Track.hh"
#include "TrackDB.hh"

// Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include <vector>

namespace detran
{

/*!
 *  \class Tracker
 *  \brief Track a mesh.
 */
class Tracker
{

public:

  typedef SP<Tracker>                   SP_tracker;
  typedef Mesh::SP_mesh                 SP_mesh;
  typedef QuadratureMOC::SP_quadrature  SP_quadrature;
  typedef Track::SP_track               SP_track;
  typedef TrackDB::SP_trackdb           SP_trackdb;

  Tracker(SP_mesh mesh, SP_quadrature quadrature);

  static SP_tracker Create(SP<Mesh>       mesh,
                           SP<Quadrature> quadrature)
  {
    SP_tracker p;
    p = new Tracker(mesh, quadrature);
    return p;
  }

  SP_trackdb trackdb() const
  {
    Require(d_trackdb);
    return d_trackdb;
  }

  SP_mesh meshmoc()
  {
    // Create the MOC mesh
    SP_mesh newmesh(new MeshMOC(d_mesh, d_trackdb));
    return newmesh;
  }

  // Normalize the track segments based on actual volumes.
  void normalize();

private:

  /// \name Private Data
  /// \{

  SP_mesh d_mesh;

  SP_quadrature d_quadrature;

  SP_trackdb d_trackdb;

  // Number azimuths per octant
  int d_number_azimuths;

  vec_dbl d_x;

  vec_dbl d_y;

  /// \}

  /// \name Implementation
  /// \{

  void generate_tracks();

  void find_starting_cell(Point enter, double tan_phi, int *IJ);

  /// \}
};

} // end namespace detran

#endif /* TRACKER_HH_ */

//---------------------------------------------------------------------------//
//              end of file Tracker.hh
//---------------------------------------------------------------------------//
