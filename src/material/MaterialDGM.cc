/*
 * MaterialDGM.cc
 *
 *  Created on: Jun 11, 2012
 *      Author: robertsj
 */

// Detran
#include "MaterialDGM.hh"

// Utilities
#include "Warning.hh"

// System headers
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran
{

MaterialDGM::MaterialDGM(int number_groups,
                         int number_materials,
                         int number_angles,
                         bool downscatter)
 : Base(number_groups, number_materials, downscatter)
 , d_number_angles(number_angles)
 , d_delta(b_number_materials,
           vec2_dbl(b_number_groups,
                    vec_dbl(d_number_angles, 0.0)))
{
  // Preconditions
  Require(d_number_angles >= 0);

  // Postconditions
  Ensure(d_delta.size()       == b_number_materials);
  Ensure(d_delta[0].size()    == b_number_groups);
  Ensure(d_delta[0][0].size() == d_number_angles);
}

void MaterialDGM::set_delta(int m, int g, int a, double v)
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  Require(a >= 0);
  Require(a < d_number_angles);
  d_delta[m][g][a] = v;
}

} // end namespace detran

