//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_IO_HDF5.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of IO_HDF5 class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_IO_HDF5_input)      \
        FUNC(test_IO_HDF5_material)   \
        FUNC(test_IO_HDF5_mesh)

// Detran headers
#include "TestDriver.hh"
#include "IO_HDF5.hh"
#include "Mesh3D.hh"

// Setup
/* ... */

using namespace detran;
using namespace detran_ioutils;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_IO_HDF5_input(int argc, char *argv[])
{

  // Create the input.
  InputDB::SP_input input(new InputDB());

  // Put some integers
  input->put<int>("dimension", 2);
  input->put<int>("number_groups", 7);


  // Put some doubles
  input->put<double>("double0_123", 0.123);
  input->put<double>("double0_234", 0.234);

  // Put some strings
  input->put<string>("bc_left",  "vacuum");
  input->put<string>("bc_right", "vacuum");

  // Put some vectors

  vec_int v_int0(2, 1);
  v_int0[1] = 2;
  vec_int v_int1(3, 3);
  v_int1[1] = 6;
  v_int1[2] = 9;

  vec_dbl v_dbl0(2, 0.0);
  v_dbl0[1] = 1.0;
  vec_dbl v_dbl1(4, 0.0);
  v_dbl1[1] = 2.0;
  v_dbl1[2] = 4.0;
  v_dbl1[3] = 6.0;

  input->put<vec_int>("vec_int_12",   v_int0);
  input->put<vec_int>("vec_int_369",  v_int1);
  input->put<vec_dbl>("vec_dbl_01",   v_dbl0);
  input->put<vec_dbl>("vec_dbl_0246", v_dbl1);

  // Create an IO_HDF5
  IO_HDF5 io("test.h5");

  // Write to the HDF5 file.  The input has the filename
  // to write out, or a default is used.
  io.write(input);
  io.close();

  // Create a new input database and read in from the
  InputDB::SP_input input2;
  input2 = io.read_input();
  return 0;
  io.close();

  //-------------------------------------------------------------------------//
  // TEST
  //-------------------------------------------------------------------------//

  // Integers
  TEST(input2->check("dimension"));
  TEST(input2->get<int>("dimension")     == 2);
  TEST(input2->check("number_groups"));
  TEST(input2->get<int>("number_groups") == 7);

  // Doubles
  TEST(input2->check("double0_123"));
  TEST(soft_equiv(input2->get<double>("double0_123"), 0.123));
  TEST(input2->check("double0_234"));
  TEST(soft_equiv(input2->get<double>("double0_234"), 0.234));

  // Strings
  TEST(input2->check("bc_left"));
  TEST(input2->get<string>("bc_left") == "vacuum");
  TEST(input2->check("bc_right"));
  TEST(input2->get<string>("bc_right") == "vacuum");

  // Integer vectors
  TEST(input2->check("vec_int_12"));
  vec_int vitest = input2->get<vec_int>("vec_int_12");
  TEST(vitest.size() == v_int0.size());
  for (int i = 0; i < vitest.size(); i++)
    TEST(vitest[i] == v_int0[i]);
  //
  TEST(input2->check("vec_int_369"));
  vitest = input2->get<vec_int>("vec_int_369");
  TEST(vitest.size() == v_int1.size());
  for (int i = 0; i < vitest.size(); i++)
    TEST(vitest[i] == v_int1[i]);

  // Double vectors
  TEST(input2->check("vec_dbl_01"));
  vec_dbl vdtest = input2->get<vec_dbl>("vec_dbl_01");
  TEST(vdtest.size() == v_dbl0.size());
  for (int i = 0; i < vdtest.size(); i++)
    TEST(soft_equiv(vdtest[i], v_dbl0[i]));
  //
  TEST(input2->check("vec_dbl_0246"));
  vdtest = input2->get<vec_dbl>("vec_dbl_0246");
  TEST(vdtest.size() == v_dbl1.size());
  for (int i = 0; i < vdtest.size(); i++)
    TEST(soft_equiv(vdtest[i], v_dbl1[i]));

  return 0;
}

int test_IO_HDF5_material(int argc, char *argv[])
{

  // Create test material
  Material::SP_material mat(new Material(2, 2));
  // material 0
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_t(0, 1, 2.0);
  mat->set_sigma_s(0, 0, 0, 1.1);
  mat->set_sigma_s(0, 1, 0, 2.1);
  mat->set_sigma_s(0, 0, 1, 3.1);
  mat->set_sigma_s(0, 1, 1, 4.1);
  // material 1
  mat->set_sigma_t(1, 0, 1.2);
  mat->set_sigma_t(1, 1, 2.2);
  mat->set_sigma_s(1, 0, 0, 1.3);
  mat->set_sigma_s(1, 1, 0, 2.3);
  mat->set_sigma_s(1, 0, 1, 3.3);
  mat->set_sigma_s(1, 1, 1, 4.3);
  mat->set_sigma_f(1, 0, 5.2);
  mat->set_sigma_f(1, 1, 6.2);
  mat->set_chi(1, 0, 1.0);
  mat->finalize();

  // Create an IO_HDF5
  IO_HDF5 io("test.h5");

  // Write to the HDF5 file.  The input has the filename
  // to write out, or a default is used.
  io.write(mat);
  io.close();

  // Create a new material object.  This must be empty.
  Material::SP_material mat2;
  mat2 = io.read_material();
  io.close();

  // Tests
  TEST(mat2);
  TEST(mat2->number_groups()    == 2);
  TEST(mat2->number_materials() == 2);
  //
  TEST(soft_equiv(mat2->sigma_t(0, 0),    1.0));
  TEST(soft_equiv(mat2->sigma_t(0, 1),    2.0));
  TEST(soft_equiv(mat2->sigma_s(0, 0, 0), 1.1));
  TEST(soft_equiv(mat2->sigma_s(0, 1, 0), 2.1));
  TEST(soft_equiv(mat2->sigma_s(0, 0, 1), 3.1));
  TEST(soft_equiv(mat2->sigma_s(0, 1, 1), 4.1));
  //
  TEST(soft_equiv(mat2->sigma_t(1, 0),    1.2));
  TEST(soft_equiv(mat2->sigma_t(1, 1),    2.2));
  TEST(soft_equiv(mat2->sigma_s(1, 0, 0), 1.3));
  TEST(soft_equiv(mat2->sigma_s(1, 1, 0), 2.3));
  TEST(soft_equiv(mat2->sigma_s(1, 0, 1), 3.3));
  TEST(soft_equiv(mat2->sigma_s(1, 1, 1), 4.3));
  TEST(soft_equiv(mat2->sigma_f(1, 0),    5.2));
  TEST(soft_equiv(mat2->sigma_f(1, 1),    6.2));
  TEST(soft_equiv(mat2->nu(1, 0),         1.0));
  TEST(soft_equiv(mat2->nu(1, 1),         1.0));
  TEST(soft_equiv(mat2->chi(1, 0),        1.0));
  TEST(soft_equiv(mat2->chi(1, 1),        0.0));

  return 0;
}

int test_IO_HDF5_mesh(int argc, char *argv[])
{

  vec_dbl cm(3, 0.0);
  cm[1] = 1.0;
  cm[2] = 2.0;
  vec_int fm(2, 2);
  vec_int mt(8, 0);
  mt[7] = 1; // Right-Top-North corner is different.

  Mesh::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create an IO_HDF5
  IO_HDF5 io("test.h5");

  // Write to the HDF5 file.  The input has the filename
  // to write out, or a default is used.
  io.write(mesh);
  io.close();
return 0;
  // Create a new mesh.
  Mesh::SP_mesh mesh2;
  mesh2 = io.read_mesh();
  io.close();
  TEST(mesh2);

  // Tests.
  TEST(mesh2->dimension()      == 3);
  TEST(mesh2->number_cells()   == 64);
  TEST(mesh2->number_cells_x() == 4);
  TEST(soft_equiv(mesh2->dx(0),   0.5));
  TEST(soft_equiv(mesh2->dy(1),   0.5));
  TEST(soft_equiv(mesh2->dz(0),   0.5));
  TEST(mesh->mesh_map_exists("MATERIAL"));
  vec_int map = mesh->mesh_map("MATERIAL");
  int ref[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1};
  for (int i = 0; i < 64; i++)
  {
    TEST(map[i] == ref[i]);
  }

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_IO_HDF5.cc
//---------------------------------------------------------------------------//
