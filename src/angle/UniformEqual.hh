//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UniformEqual.hh
 * \author robertsj
 * \date   May 22, 2012
 * \brief  UniformEqual class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef UNIFORMEQUAL_HH_
#define UNIFORMEQUAL_HH_

// Detran
#include "Quadrature.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class Uniform_Equal_3D
 * \brief 2D/3D Uniform, Equal Weight (UEn) quadrature class.
 *
 * As mentioned in \ref LevelSymmetric, a fundamental problem inherent
 * to LQn quadrature is presence of negative weights at high order.  These
 * weights produce unphysical solutions (and may inhibit convergence).  For
 * problems where an increased quadrature order (i.e. more angles) is
 * required to study convergence or simply to get better answers, we
 * require an arbitrarily high order, positive weight quadrature.
 *
 * Here, we implement the "uniform, equal weight" (UEn) quadrature of Carew
 * and Zamonsky.  The basic idea is to choose uniform azimuthal divisions
 * and uniform polar cosines.  The result is product quadrature that
 * can be extended on-the-fly to arbitrary numbers of angles.
 *
 * Here, the quadrature order defines the number of polar angles.
 * We take the number of azimuthal angles to be twice this
 * number.  Hence, the total number of angles is
 * twice the number of polar angles squared.
 *
 *
 * \refs
 * - Carew, J. and Zamonsky. G., <em>Nuclear Science and Engineering</em>
 *     <b>131</b>, 199-207 (1999).
 *
 */
//===========================================================================//
class UniformEqual : public Quadrature
{

public:

  typedef SP<UniformEqual>  SP_quadrature;
  typedef SP<Quadrature>    SP_base;

  /*!
   *  \brief Constructor.
   *
   *  \param    order       Quadrature order.
   *  \param    dim         Problem dimension
   */
  UniformEqual(int order, int dim);

  /*!
   *  \brief SP Constructor.
   */
  static SP<Quadrature> Create(int order, int dim)
  {
    SP_quadrature p;
    p = new UniformEqual(order, dim);
    return p;
  }

  void display() const;

};

} // end namespace detran


#endif /* UNIFORMEQUAL_HH_ */
