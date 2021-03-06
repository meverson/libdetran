//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   State.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  State class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef STATE_HH_
#define STATE_HH_

// Other libtran headers
#include "Definitions.hh"
#include "SP.hh"
#include "InputDB.hh"
#include "Quadrature.hh"
#include "Mesh.hh"

// System
#include <vector>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class State
 * \brief Problem state.
 *
 * The state of the problem at any point in the solution process can be
 * defined by the unknown flux moments (scalar flux and higher moments),
 * as these quantities are sufficient to describe reaction rates, which is
 * typically what we need (e.g. doses or fission rates).  For eigenvalue
 * problems, keff is also important.
 *
 * When needed, the underlying data array for a group can be accessed as
 * demonstrated by the following (I <b>think!</b>).
 * \code
     // Get moments by reference.
     moments_type moments = state.moments();
     double *localgrouparray;
     localgrouparray = &moments[g][0];
 * \endcode
 * This would be useful e.g. for filling a PETSc Vec object to use in one
 * of their Krylov schemes (rather than copying into another array).
 *
 * \todo Test whether referencing a Moments_Field object at [0] actually
 *       yields the array underneath.
 *
 * Relevant input entries:
 * - number_groups (int)
 * - store_angular_flux (int)
 */
//---------------------------------------------------------------------------//
class State : public Object
{

public:

  typedef SP<State>                                     SP_state;
  typedef InputDB::SP_input                             SP_input;
  typedef Mesh::SP_mesh                                 SP_mesh;
  typedef Quadrature::SP_quadrature                     SP_quadrature;
  typedef vec_dbl                                       moments_type;
  typedef std::vector<moments_type>                     vec_moments_type;
  typedef std::vector<moments_type>                     group_moments_type;
  typedef vec_dbl                                       angular_flux_type;
  typedef std::vector<std::vector<angular_flux_type> >  vec_angular_flux_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    quadrature  Angular quadrature.
   */
  State(SP_input        input,
        SP_mesh         mesh,
        SP_quadrature   quadrature);

  /*!
   *  \brief Constructor.
   *
   *  This variant has no quadrature and can be used in
   *  diffusion problems.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   */
  State(SP_input        input,
        SP_mesh         mesh);

  /// SP constructor.
  static SP_state Create(SP<detran::InputDB>       input,
                         SP<detran::Mesh>          mesh,
                         SP<detran::Quadrature>    quadrature)
  {
    SP_state p;
    p = new State(input, mesh, quadrature);
    return p;
  }

  /// SP constructor with no quadrature.
  static SP_state Create(SP<detran::InputDB>       input,
                         SP<detran::Mesh>          mesh)
  {
    SP_state p;
    p = new State(input, mesh);
    return p;
  }

  /*
   *  \brief Const accessor to a group moments field.
   *
   *  \param    g   Group of field requested.
   *  \return       Constant reference to group moment vector.
   */
  const moments_type& phi(int g) const
  {
    Require(g >= 0);
    Require(g < d_number_groups);
    return d_moments[g];
  }

  /*
   *  \brief Mutable accessor to a group moments field.
   *
   *  This is to be used for copying, i.e.
   *  \code
   *    State::moments_type new_phi;
   *    // compute new_phi and update
   *    moments->phi(g) = new_phi;
   *  \endcode
   *
   *  \param    g   Group of field requested.
   *  \return       Mutable reference to group moment vector.
   */
  moments_type& phi(int g)
  {
    // Cast away return type
    return const_cast<moments_type&>
    (
      // Add const to *this's type and call const version
      static_cast<const State*>(this)->phi(g)
    );

  }

  /*
   *  \brief Const accessor to all group moments.
   */
  const group_moments_type& all_phi() const
  {
    return d_moments;
  }

  /*
   *  \brief Mutable accessor to all group moments.
   */
  group_moments_type& all_phi()
  {
    // Cast away return type
    return const_cast<group_moments_type&>
    (
      // Add const to *this's type and call const version
      static_cast<const State*>(this)->all_phi()
    );

  }

  /*!
   *   \brief Set a group moment vector.
   *   \param   f   User-defined moment vector.
   */
  void set_moments(int g, std::vector<double> f)
  {
    d_moments[g] = f;
  }

  /*
   *  \brief Const accessor to a group angular flux.
   *
   *  \param    g   Group of field requested.
   *  \param    o   Octant
   *  \param    a   Angle within octant
   *  \return       Constant reference to group angular flux vector.
   */
  const angular_flux_type& psi(int g, int o, int a) const
  {
    Require(d_store_angular_flux);
    Require(d_angular_flux.size() > 0);
    Require(o >= 0);
    Require(o < d_quadrature->number_octants());
    Require(a >= 0);
    Require(a < d_quadrature->number_angles_octant());
    Require(g >= 0);
    Require(g < d_number_groups);
    int angle = d_quadrature->index(o, a);
    return d_angular_flux[g][angle];
  }

  /*
   *  \brief Mutable accessor to a group angular flux.
   *
   *  \param    g   Group of field requested.
   *  \param    o   Octant
   *  \param    a   Angle within octant
   *  \return       Mutable reference to group angular flux vector.
   */
  angular_flux_type& psi(int g, int o, int a)
  {
    // Cast away return type
    return const_cast<angular_flux_type&>
    (
      // Add const to *this's type and call const version
      static_cast<const State*>(this)->psi(g, o, a)
    );
  }

  double eigenvalue() const
  {
    return d_eigenvalue;
  }

  void set_eigenvalue(double v)
  {
    d_eigenvalue = v;
  }

  SP_input get_input()
  {
    return d_input;
  }

  SP_mesh get_mesh()
  {
    return d_mesh;
  }

  int moments_size() const
  {
    return d_moments[0].size();
  }

  SP_quadrature get_quadrature()
  {
    return d_quadrature;
  }

  int number_groups() const
  {
    return d_number_groups;
  }

  bool store_angular_flux() const
  {
    return d_store_angular_flux;
  }

  bool is_valid() const
  {
    return true;
  }

private:

  /// Input database
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// Boundary fluxes
  //SP_boundary d_boundary;

  /// Number of energy groups
  int d_number_groups;

  /// Cell-center scalar flux moments, [energy, (space-moment)]
  vec_moments_type d_moments;

  /// Cell-center angular flux, [energy, angle, (space)]
  vec_angular_flux_type d_angular_flux;

  /// k-eigenvalue
  double d_eigenvalue;

  /// Store the angular flux?
  bool d_store_angular_flux;

};

} // end namespace detran

#endif /* STATE_HH_ */

//---------------------------------------------------------------------------//
//              end of State.hh
//---------------------------------------------------------------------------//
