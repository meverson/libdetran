//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_materials.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran materials.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"

%include "Material.hh"

%template(MaterialSP) detran::SP<detran::Material>;

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//