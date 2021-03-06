//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_CoarseMesh.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Equation_DD_2D
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_CoarseMesh)

// Detran headers
#include "TestDriver.hh"
#include "CoarseMesh.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

// Setup
#include "coarsemesh_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_CoarseMesh(int argc, char *argv[])
{

  //-------------------------------------------------------------------------//
  // SHARED DATA
  //-------------------------------------------------------------------------//

  // Coarse mesh dx's in all directions
  double delta[] = {3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  //-------------------------------------------------------------------------//
  // 1D TESTS
  //-------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh1 = coarsemesh_1d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh1 = coarsemesh1->get_coarse_mesh();

  TEST(cmesh1->number_cells()   == 7);
  TEST(cmesh1->number_cells_x() == 7);
  TEST(cmesh1->number_cells_y() == 1);
  TEST(cmesh1->number_cells_z() == 1);
  for (int i = 0; i < cmesh1->number_cells(); i++)
  {
    TEST(soft_equiv(cmesh1->dx(i), delta[i]));
  }
  TEST(soft_equiv(cmesh1->dy(0), 1.0));
  TEST(soft_equiv(cmesh1->dz(0), 1.0));

  //-------------------------------------------//
  // 2D TESTS
  //-------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh2 = coarsemesh_2d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh2 = coarsemesh2->get_coarse_mesh();

  TEST(cmesh2->number_cells()   == 49);
  TEST(cmesh2->number_cells_x() == 7);
  TEST(cmesh2->number_cells_y() == 7);
  TEST(cmesh2->number_cells_z() == 1);
  for (int i = 0; i < cmesh2->number_cells_x(); i++)
  {
    TEST(soft_equiv(cmesh2->dx(i), delta[i]));
    TEST(soft_equiv(cmesh2->dy(i), delta[i]));
  }
  TEST(soft_equiv(cmesh2->dz(0), 1.0));

  //-------------------------------------------//
  // 3D TESTS
  //-------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh3 = coarsemesh_3d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh3 = coarsemesh3->get_coarse_mesh();

  TEST(cmesh3->number_cells()   == 343);
  TEST(cmesh3->number_cells_x() == 7);
  TEST(cmesh3->number_cells_y() == 7);
  TEST(cmesh3->number_cells_z() == 7);
  for (int i = 0; i < cmesh3->number_cells_x(); i++)
  {
    TEST(soft_equiv(cmesh3->dx(i), delta[i]));
    TEST(soft_equiv(cmesh3->dy(i), delta[i]));
    TEST(soft_equiv(cmesh3->dz(i), delta[i]));
  }

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_CoarseMesh.cc
//---------------------------------------------------------------------------//
