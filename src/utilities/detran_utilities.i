//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_utilities.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran utilities.
 */
//---------------------------------------------------------------------------//

%module(package="detran") utilities
%{
#include <stddef.h>
#include "utilities/utilities_export.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/SP.hh"
%}

// STL
%include std_container.i
%include std_map.i
%include std_string.i
%include std_vector.i

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define UTILITIES_EXPORT
#define UTILITIES_TEMPLATE_EXPORT(...)
#define UTILITIES_INSTANTIATE_EXPORT(...)
}

// Dummy SP interface to stop some annoying warnings from SWIG
namespace detran_utilities
{

template<class T>
class SP 
{
public:
  SP();
  inline explicit SP(T *p_in);
  template<class X>
  inline explicit SP(X *px_in);
  inline SP(const SP<T> &sp_in);
  template<class X>
  inline SP(const SP<X> &spx_in);
  ~SP();
  T* operator->() const;
  T& operator*() const;
  T* bp() const;
  operator bool() const;
  bool operator==(const T *p_in) const;
  bool operator!=(const T *p_in) const;
  bool operator==(const SP<T> &sp_in) const;
  bool operator!=(const SP<T> &sp_in) const;
};

} // end namespace detran_utilities

// Vectors
%include "Definitions.hh"
namespace std
{
  %template(vec_int)      vector<int>;
  %template(vec2_int)     vector<vector<int> >;
  %template(vec3_int)     vector<vector<vector<int> > >;
  %template(vec_dbl)      vector<double>;
  %template(vec2_dbl)     vector<vector<double> >;
  %template(vec3_dbl)     vector<vector<vector<double> > >;
  %template(vec_size_t)   vector<unsigned int>;
  %template(vec2_size_t)  vector<vector<unsigned int> >;
  %template(vec3_size_t)  vector<vector<vector<unsigned int> > >;
}

%include "InputDB.hh"
%template(InputDBSP)    detran_utilities::SP<detran_utilities::InputDB>;
%template(get_int)      detran_utilities::InputDB::get<int>;
%template(get_dbl)      detran_utilities::InputDB::get<double>;
%template(get_vec_int)  detran_utilities::InputDB::get<detran_utilities::vec_int>;
%template(get_vec_dbl)  detran_utilities::InputDB::get<detran_utilities::vec_dbl>;
%template(get_str)      detran_utilities::InputDB::get<std::string>;
%template(get_spdb)     detran_utilities::InputDB::get<detran_utilities::InputDB::SP_input>;
%template(put_int)      detran_utilities::InputDB::put<int>;
%template(put_dbl)      detran_utilities::InputDB::put<double>;
%template(put_vec_int)  detran_utilities::InputDB::put<detran_utilities::vec_int>;
%template(put_vec_dbl)  detran_utilities::InputDB::put<detran_utilities::vec_dbl>;
%template(put_str)      detran_utilities::InputDB::put<std::string>;
%template(put_spdb)     detran_utilities::InputDB::put<detran_utilities::InputDB::SP_input>;

%include "MathUtilities.hh"

// Anyhere in C/C++ that we need (argc, argv)
%typemap(in) (int argc, char *argv[]) 
{
  // Check if is a list
  if (PyList_Check($input)) 
  {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) 
    {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $2[i] = PyString_AsString(PyList_GetItem($input,i));
      else 
      {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($2);
        return NULL;
      }
    }
    $2[i] = 0;
  } 
  else 
  {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(freearg) (int argc, char *argv[]) {
  free((char *) $2);
}

// see http://embedded.eecs.berkeley.edu/Alumni/pinhong/scriptEDA/pyTypemapFAQ.html#32

//---------------------------------------------------------------------------//
//              end of detran_utilities.i
//---------------------------------------------------------------------------//
