/*
 * MaterialDGM.i.hh
 *
 *  Created on: Jun 11, 2012
 *      Author: robertsj
 */

#ifndef MATERIALDGM_I_HH_
#define MATERIALDGM_I_HH_

#include "MaterialDGM.hh"

namespace detran
{

inline double MaterialDGM::delta(int m, int g, int a) const
{
  Require(m >= 0);
  Require(m < b_number_materials);
  Require(g >= 0);
  Require(g < b_number_groups);
  Require(a >= 0);
  Require(a < d_number_angles);
  return d_delta[m][g][a];
}

} // end namespace detran

#endif /* MATERIALDGM_I_HH_ */
