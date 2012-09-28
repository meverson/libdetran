//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   callow.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow library.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"
%include "callow_config.hh"

//---------------------------------------------------------------------------//
// initialization
//---------------------------------------------------------------------------//

void callow_initialize(int argc, char *argv[]);
void callow_finalize();

//---------------------------------------------------------------------------//
// setup for numerical arrays
//---------------------------------------------------------------------------//

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
  import_array();
%}
%numpy_typemaps(double, NPY_DOUBLE, int)
%numpy_typemaps(in,     NPY_INT, int)
%apply (int* IN_ARRAY1, double* IN_ARRAY2, int DIM1) 
       {(int *j, double *v, int n)}
%apply (int* IN_ARRAY1, int* IN_ARRAY2, double* IN_ARRAY3, int DIM1) 
       {(int *i, int *j, double* v, int n)}

//---------------------------------------------------------------------------//
// vector
//---------------------------------------------------------------------------//

%include "vector/Vector.i"

//---------------------------------------------------------------------------//
// matrix
//---------------------------------------------------------------------------//

%include "matrix/Matrix.i"

//%include "matrix/MatrixBase.hh"
//%template(MatrixBaseDouble)     callow::MatrixBase<double>;
//%template(MatrixBaseDoubleSP)   detran_utilities::SP<callow::MatrixBase<double> >;
//%include "matrix/Matrix.hh"
//%template(MatrixDouble)     callow::Matrix<double>;
//%template(MatrixDoubleSP)   detran_utilities::SP<callow::Matrix<double> >;

//---------------------------------------------------------------------------//
// linear solver
//---------------------------------------------------------------------------//
