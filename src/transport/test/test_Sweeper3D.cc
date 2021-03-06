//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Sweeper3D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of test_Sweeper3D
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Sweeper3D_basic)

// Detran headers
#include "TestDriver.hh"
#include "Sweeper3D.hh"
#include "ConstantSource.hh"
#include "LevelSymmetric.hh"
#include "Mesh3D.hh"
//#include "SiloOutput.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"

using namespace detran;
//using namespace detran_ioutils;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Sweeper3D_basic(int argc, char *argv[])
{
  typedef Sweeper3D<Equation_DD_3D> Sweeper_T;

  vec_int fm(1, 3);
  vec_dbl cm(2, 0.0);
  vec_int fmz(1, 1);
  cm[1] = 1.0;
  vec_int mt(1, 0);

  // Test fixtures
  Sweeper_T::SP_mesh mesh       = Mesh3D::Create(fm, fm, fm, cm, cm, cm, mt);
  Sweeper_T::SP_material mat    = material_fixture_1g();
  Sweeper_T::SP_quadrature quad = LevelSymmetric::Create(2, 3);

  // Input
  Sweeper_T::SP_input input(new InputDB());
  input->put<int>("number_groups", 1);

  // State
  Sweeper_T::SP_state state(new State(input, mesh, quad));

  // Boundary
  Sweeper_T::SP_boundary
    bound(new Sweeper_T::Boundary_T(input, mesh, quad));

  // Moment to Discrete
  MomentToDiscrete<_3D>::SP_MtoD m2d(new MomentToDiscrete<_3D>(0));
  m2d->build(quad);

  // External
  ConstantSource::SP_source q_e(new ConstantSource(mesh, quad, 1, 1.0));

  // Sweep source
  Sweeper_T::SP_sweepsource
    source(new SweepSource<_3D>(state, mesh, quad, mat, m2d));
  source->set_moment_source(q_e);
  source->build_fixed(0);

  // Sweeper
  Sweeper_T sweeper(input, mesh, mat, quad, state, bound, source);

  // Sweep
  State::moments_type phi(mesh->number_cells(0), 0.0);
  sweeper.setup_group(0);
  sweeper.sweep(phi);
  state->set_moments(0, phi);
  for (int i = 0; i < phi.size(); i++)
  {
    cout << " i = " << i << " phi = " << phi[i] << endl;
  }

//  // Create the SiloOutput.
//  SiloOutput out(mesh);
//  cout << "...created out." << endl;
//
//  // Open the file, and write the mesh.
//  TEST(out.initialize("test.silo"));
//  cout << "...initialized out." << endl;
//
//  // Write out the scalar flux.
//  TEST(out.write_scalar_flux(state));
//  cout << "...wrote scalar flux." << endl;
//
//  out.finalize();
  return 0;
}
