//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_fixture.hh
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef MESH_FIXTURE_HH_
#define MESH_FIXTURE_HH_

// Detran includes
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "PinCell.hh"
#include "Assembly.hh"
#include "Core.hh"

// Detran utilities
#include "DBC.hh"
#include "Definitions.hh"

// System includes

namespace detran_test
{

typedef detran::Mesh::SP_mesh SP_mesh;
typedef detran::PinCell::SP_pincell SP_pincell;
typedef detran::Assembly::SP_assembly SP_assembly;
typedef detran::Core::SP_core SP_core;
/*!
 *  \brief Create a common 1D mesh for transport tests.
 */
static SP_mesh mesh_1d_fixture()
{
  // Define the mesh coarse and fine meshes, and the
  // coarse mesh material map.
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;
  detran::vec_int fm(2, 5);
  detran::vec_int mat_map(2, 1);
  mat_map[0] = 0;

  // Create the new database.
  SP_mesh mesh;
  mesh = new detran::Mesh1D(fm, cm, mat_map);

  // Ensure a valid mesh.
  Ensure(mesh->is_valid());

  // Return the fixture.
  return mesh;
}


/*!
 *  \brief Create a common 2D mesh for transport tests.
 *
 *  This is a square mesh with 2x2 coarse meshes of
 *  width 10 cm sides, each with 10x10 fine meshes.
 *  The bottom left coarse mesh has material 0, while
 *  the rest have material 1.
 */
static SP_mesh mesh_2d_fixture(int id = 0)
{
  // Mesh pointer
  SP_mesh mesh;

  if (id == 0)
  {
    // Define the mesh coarse and fine meshes, and the
    // coarse mesh material map.
    detran::vec_dbl cm(3, 0.0);
    cm[1] = 10.0;
    cm[2] = 20.0;
    detran::vec_int fm(2, 10);
    // All maps are stored in a single dimension.  The
    // memory order follows the transpose of the physical
    // order.  Here, the entries are (i=0,j=0), (i=0,j=1),
    // etc.  Hence, all entries are 1 but the first.
    detran::vec_int mat_map(4, 1);
    mat_map[0] = 0;

    // Construct the mesh.
    mesh = new detran::Mesh2D(fm, fm, cm, cm, mat_map);

    // Ensure a valid mesh.
    Ensure(mesh->is_valid());

    // Return the fixture.
    return mesh;
  }
  else if (id == 1)
  {
    detran::vec_dbl cm(2, 0.0);
    cm[1] = 1.0;
    detran::vec_int fm(1, 3);
    detran::vec_int mat_map(1, 1);
    // Construct the mesh.
    mesh = new detran::Mesh2D(fm, fm, cm, cm, mat_map);
    // Ensure a valid mesh.
    Ensure(mesh->is_valid());

  }

  // Return the fixture.
  return mesh;
}

/*!
 *  \brief Create a common 3D mesh for transport tests.
 */
static SP_mesh mesh_3d_fixture()
{
  // Define the mesh coarse and fine meshes, and the
  // coarse mesh material map.
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 5.0;
  cm[2] = 10.0;
  detran::vec_int fm(2, 5);
  detran::vec_int mat_map(8, 1);
  mat_map[0] = 0;

  // Create the new database.
  SP_mesh mesh;
  mesh = new detran::Mesh3D(fm, fm, fm, cm, cm, cm, mat_map);

  // Ensure a valid mesh.
  Ensure(mesh->is_valid());

  // Return the fixture.
  return mesh;
}

// Defines a two region pin cell for use with material_fixture_2g
static SP_pincell pincell_fixture()
{
  // PinCell(double pitch, vec_dbl radii,
  //         vec_int mat_map, bool fuel_flag = true);
  double pitch = 1.26;
  detran::vec_dbl radii(1, 0.54);
  detran::vec_int mat_map(2, 0);
  mat_map[0] = 1; // pin (fuel I)
  mat_map[1] = 0; // mod
  bool fuel_flag = true; // this is a fuel pin (i.e. not a moderator box)
  SP_pincell pin(new detran::PinCell(pitch, radii, mat_map, fuel_flag));
  pin->meshify(7, true);
  return pin;
}

// Defines a 2x2 assembly for use with material_fixture_2g
static SP_assembly assembly_fixture()
{
  // Define the 2x2 assembly
  SP_assembly assembly(new detran::Assembly(2));

  // Get the pin fixture
  SP_pincell pin = pincell_fixture();

  // Add the pin to the assembly
  assembly->add_pincell(pin);

  // Create a pincell map (all the same pin)
  detran::vec_int pincell_map(4, 0);

  // Add the map
  assembly->finalize(pincell_map);

  return assembly;
}

// Defines a 2x2 core for use with material_fixture_2g
static SP_core core_fixture()
{
  // Define the 3x3 assembly
  SP_core core(new detran::Core(2));

  // Get the assembly fixture
  SP_assembly assembly = assembly_fixture();

  // Add the pin to the assembly
  core->add_assembly(assembly);

  // Create a pincell map (all the same pin)
  detran::vec_int assembly_map(4, 0);

  // Add the map
  core->finalize(assembly_map);

  return core;
}

} // end namespace detran_test

#endif /* MESH_FIXTURE_HH_ */

//---------------------------------------------------------------------------//
//              end of mesh_fixture.hh
//---------------------------------------------------------------------------//
