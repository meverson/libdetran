//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SiloOutput.cc
 * \brief  SiloOutput 
 * \author Jeremy Roberts
 * \date   Jul 27, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_SILO

// Postprocess
#include "SiloOutput.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"

// System
#include <string>
#include <stdio.h>

namespace detran_ioutils
{

SiloOutput::SiloOutput(SP_mesh mesh)
  : d_mesh(mesh)
  , d_initialized(false)
{
  // Preconditions
  Require(d_mesh);
}

SiloOutput::~SiloOutput()
{
  if(d_initialized) finalize();
}

bool SiloOutput::initialize(std::string filename)
{
  // Simply return if already initialized.
  if (d_initialized) return false;

  // Create the file
  d_silofile = DBCreate((char *)filename.c_str(), DB_CLOBBER,
                         DB_LOCAL, NULL, DB_PDB );
  Insist(d_silofile, "Error creating Silo file.");
  d_initialized = true;

  // Problem dimension.  Note, while 1-D Silo makes little sense,
  // we include it as 2-D for simplicity.
  d_dimension = d_mesh->dimension();
  if (d_dimension == 1) d_dimension = 2;

  // Names of the coordinates
  char *coordnames[3];

  // Array of coordinate arrays
  double *coords[3];

  // Name the coordinates
  char c1[5] = "x";
  coordnames[0] = c1;
  char c2[5] = "y";
  coordnames[1] = c2;
  char c3[5] = "z";
  coordnames[2] = c3;

  // Set coordinate dimensions
  d_dims[0] = d_mesh->number_cells_x() + 1;
  d_dims[1] = d_mesh->number_cells_y() + 1;
  d_dims[2] = d_mesh->number_cells_z() + 1;

  // Create the coordinate vectors (mesh edges)
  detran::vec_dbl x(d_mesh->number_cells_x() + 1, 0.0);
  detran::vec_dbl y(d_mesh->number_cells_y() + 1, 0.0);
  detran::vec_dbl z(d_mesh->number_cells_z() + 1, 0.0);
  for (int i = 0; i < d_mesh->number_cells_x(); i++)
    x[i + 1] = x[i] + d_mesh->dx(i);
  for (int i = 0; i < d_mesh->number_cells_y(); i++)
    y[i + 1] = y[i] + d_mesh->dy(i);
  for (int i = 0; i < d_mesh->number_cells_z(); i++)
    z[i + 1] = z[i] + d_mesh->dz(i);
  coords[0] = &x[0];
  coords[1] = &y[0];
  coords[2] = &z[0];

  // Add the mesh
  DBPutQuadmesh(d_silofile, "mesh", NULL, coords, d_dims,
                d_dimension, DB_DOUBLE, DB_COLLINEAR, NULL);

  // Decrement the dimensions for zonal quantities.
  d_dims[0]--;
  d_dims[1]--;
  d_dims[2]--;

}

void SiloOutput::finalize()
{
  if (d_initialized)
  {
    DBClose(d_silofile);
    d_initialized = false;
  }
}

bool SiloOutput::write_mesh_map(const std::string &key)
{
  using std::cout;
  using std::endl;

  if (!d_initialized)
  {
    cout << "Silo not initialized; not writing mesh map" << endl;
    return false;
  }

  // Get the mesh map
  if (!d_mesh->mesh_map_exists(key))
  {
    cout << "Mesh map doesn't exist; not writing mesh map" << endl;
    return false;
  }
  detran::vec_int map = d_mesh->mesh_map(key);

  // Write to silo
  DBPutQuadvar1(d_silofile, key.c_str(), "mesh", &map[0],
                d_dims, d_dimension, NULL, 0, DB_INT,
                DB_ZONECENT, NULL);

  return true;
}

bool SiloOutput::write_scalar_flux(SP_state state)
{
  // Preconditions
  Require(state);

  if (!d_initialized)
  {
    std::cout << "Silo not initialized; not writing fluxes" << std::endl;
    return false;
  }

  for (int g = 0; g < state->number_groups(); g++)
  {

    // Write the group
    char buffer[14];
    sprintf(buffer, "group%i", g);

    // Get a pointer to the group flux
    double *phi = &state->phi(g)[0];

    // Write to silo
    DBPutQuadvar1(d_silofile, buffer, "mesh", phi,
                  d_dims, d_dimension, NULL, 0, DB_DOUBLE,
                  DB_ZONECENT, NULL);
  }

  return true;
}

bool SiloOutput::write_angular_flux(SP_state state, SP_quadrature quad)
{
  // Preconditions
  Require(state);
  Require(quad);

  if (!d_initialized)
  {
    std::cout << "Silo not initialized; not writing fluxes" << std::endl;
    return false;
  }
  if (!state->store_angular_flux())
  {
    std::cout << "Angular fluxes not stored; not writing angular fluxes" << std::endl;
    return false;
  }

  for (int g = 0; g < state->number_groups(); g++)
  {
    for (int o = 0; o < quad->number_octants(); o++)
    {
      for (int a = 0; a < quad->number_angles_octant(); a++)
      {
        // Example: g000_o0_a000
        char buffer[14];
        sprintf(buffer, "g%i_o%i_a%i", g, o, a);

        // Get a pointer to the group flux
        double *psi = &state->psi(g, o, a)[0];

        // Write to silo
        DBPutQuadvar1(d_silofile, buffer, "mesh", psi,
                      d_dims, d_dimension, NULL, 0, DB_DOUBLE,
                      DB_ZONECENT, NULL);
      }
    }
  }

  return true;
}

} // end namespace detran_ioutils


#endif // DETRAN_ENABLE_SILO

//---------------------------------------------------------------------------//
//              end of file SiloOutput.cc
//---------------------------------------------------------------------------//
