//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_materials.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran materials.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"

%include "Material.hh"
%include "MaterialDGM.hh"

%template(MaterialSP)    detran::SP<detran::Material>;
%template(MaterialDGMSP) detran::SP<detran::MaterialDGM>;

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//