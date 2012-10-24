//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SourceIteration.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  SourceIteration class definition.
 */
//---------------------------------------------------------------------------//

#include "WGSolverSI.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
WGSolverSI<D>::WGSolverSI(SP_state                  state,
                          SP_material               material,
                          SP_quadrature             quadrature,
                          SP_boundary               boundary,
                          const vec_externalsource &q_e,
                          SP_fissionsource          q_f,
                          bool                      multiply)
  : Base(state, material, boundary, quadrature, q_e, q_f, multiply)
{
  d_sweeper->set_update_boundary(true);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of SourceIteration.cc
//---------------------------------------------------------------------------//