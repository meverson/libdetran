//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ChebyshevDPN.cc
 *  @brief  ChebyshevDPN
 *  @author Jeremy Roberts
 *  @date   Oct 11, 2012
 */
//---------------------------------------------------------------------------//

#include "ChebyshevDPN.hh"
#include "GenerateGaussLegendre.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
ChebyshevDPN::ChebyshevDPN(const size_t dim,
                           const size_t na,
                           const size_t np)
  : ProductQuadrature(dim, na, np, "ChebyshevDPN")
{

  //-------------------------------------------------------------------------//
  // AZIMUTH QUADRATURE
  //-------------------------------------------------------------------------//

  // Chebyshev points are just equally spaced phi's in the first quadrant
  using detran_utilities::pi;
  for (int i = 0; i < na; ++i)
  {
    d_phi[i] = 0.25 * (2*i + 1)*pi / na;
    d_cos_phi[i] = std::cos(d_phi[i]);
    d_sin_phi[i] = std::sin(d_phi[i]);
    d_azimuth_weight[i] = 0.5 * pi / na;
  }

  //-------------------------------------------------------------------------//
  // POLAR QUADRATURE
  //-------------------------------------------------------------------------//

  // temporary arrays
  vec_dbl x(np, 0.0);
  vec_dbl w(np, 0.0);

  // generate parameters
  generate_gl_parameters(np, x, w);

  // fill array
  for (int i = 0; i < np; ++i)
  {
    d_cos_theta[i]    = 0.5*x[i] + 0.5; // shift and scale
    d_sin_theta[i]    = std::sqrt(1.0 - d_cos_theta[i]*d_cos_theta[i]);
    d_polar_weight[i] = 0.5*w[i];
  }

  //-------------------------------------------------------------------------//
  // PRODUCT QUADRATURE
  //-------------------------------------------------------------------------//

  double scale = 1.0;
  if (d_dimension == 2) scale = 2.0;
  double weight_tot = 0.0;
  int n = 0;
  for (int a = 0; a < na; ++a)
  {
    for (int p = 0; p < np; ++p, ++n)
    {
      d_mu[n]     = d_cos_theta[p] * d_cos_phi[a];
      d_eta[n]    = d_cos_theta[p] * d_sin_phi[a];
      d_xi[n]     = d_sin_theta[p];
      d_weight[n] = scale * d_polar_weight[p] * d_azimuth_weight[a];
      weight_tot += d_weight[n];
    } // end polar loop
  } // end azimuth loop

}

//---------------------------------------------------------------------------//
ChebyshevDPN::SP_quadrature
ChebyshevDPN::Create(const size_t dim,
                     const size_t na,
                     const size_t np)
{
  SP_quadrature p(new ChebyshevDPN(dim, na, np));
  return p;
}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of file ChebyshevDPN.cc
//---------------------------------------------------------------------------//