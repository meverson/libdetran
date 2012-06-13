/*
 * Material.i.hh
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#ifndef MATERIAL_I_HH_
#define MATERIAL_I_HH_

namespace detran
{

double Material::sigma_t(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_sigma_t[m][g];
}

double Material::sigma_a(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_sigma_a[m][g];
}

double Material::nu_sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_nu_sigma_f[m][g];
}

double Material::sigma_f(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_sigma_f[m][g];
}

double Material::nu(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_nu[m][g];
}

double Material::chi(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_chi[m][g];
}

double Material::sigma_s(int m, int g, int gp) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  Require(gp >= 0);
  Require(gp < b_number_groups);
  return b_sigma_s[m][g][gp];
}

double Material::diff_coef(int m, int g) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  return b_diff_coef[m][g];
}

///

double Material::delta(int m, int g, int a) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  Require(a >= 0);
  Require(a < b_number_angles);
  return d_delta[m][g][a];
}

} // end namespace detran

#endif /* MATERIAL_I_HH_ */
