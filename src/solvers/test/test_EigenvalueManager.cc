//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_EigenvalueManager.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  test_EigenvalueManager class definition.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                           \
        FUNC(test_EigenvalueManager_1D)     \
        FUNC(test_EigenvalueManager_2D)     \
        FUNC(test_EigenvalueManager_3D)

// Detran headers
#include "TestDriver.hh"
#include "EigenvalueManager.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "callow/utils/Initialization.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "external_source/test/external_source_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

SP_material test_EigenvalueManager_material()
{
  // Material (kinf = 1.0)
  SP_material mat(new Material(1, 1));
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_s(0, 0, 0, 0.5);
  mat->set_sigma_f(0, 0, 0.5);
  mat->set_chi(0, 0, 1.0);
  mat->set_diff_coef(0, 0, 0.33);
  mat->compute_sigma_a();
  mat->finalize();
  return mat;
}

SP_mesh test_EigenvalueManager_mesh(int d)
{
  vec_dbl cm(2, 0.0); cm[1] = 5.0;
  vec_int fm(1, 10);
  vec_int mat_map(1, 0);
  SP_mesh mesh;
  if (d == 1) mesh = new Mesh1D(fm, cm, mat_map);
  if (d == 2) mesh = new Mesh2D(fm, fm, cm, cm, mat_map);
  if (d == 3) mesh = new Mesh3D(fm, fm, fm, cm, cm, cm, mat_map);
  Assert(mesh);
  return mesh;
}

InputDB::SP_input test_EigenvalueManager_input()
{
  InputDB::SP_input inp(new InputDB());
  inp->put<int>("number_groups", 1);
  inp->put<string>("problem_type", "eigenvalue");
  inp->put<string>("equation", "diffusion");
  inp->put<string>("bc_west", "reflect");
  inp->put<string>("bc_east", "reflect");
  inp->put<string>("bc_south", "reflect");
  inp->put<string>("bc_north", "reflect");
  inp->put<double>("inner_tolerance", 1e-17);
  inp->put<int>("inner_print_level", 0);
  inp->put<int>("outer_print_level", 0);
  inp->put<int>("eigen_print_level", 2);
  inp->put<int>("eigen_print_interval", 5);
  inp->put<int>("eigen_max_iters", 100);
  inp->put<int>("quad_number_polar_octant", 2);
  return inp;
}

template <class D>
int test_EigenvalueManager_T()
{
  typedef EigenvalueManager<D>       Manager_T;

  // Input, materials, and mesh
  typename Manager_T::SP_input input = test_EigenvalueManager_input();
  SP_material mat = test_EigenvalueManager_material();
  SP_mesh mesh = test_EigenvalueManager_mesh(D::dimension);

  // Manager
  Manager_T manager(input, mat, mesh);

  // Solve.
  bool flag = manager.solve();
  TEST(flag);

  return 0;
}

int test_EigenvalueManager_1D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_1D>();
}

int test_EigenvalueManager_2D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_2D>();
}

int test_EigenvalueManager_3D(int argc, char *argv[])
{
  return test_EigenvalueManager_T<_3D>();
}

//---------------------------------------------------------------------------//
//              end of test_EigenvalueManager.cc
//---------------------------------------------------------------------------//