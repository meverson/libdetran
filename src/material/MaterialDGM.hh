/*
 * MaterialDGM.hh
 *
 *  Created on: Jun 11, 2012
 *      Author: robertsj
 */

#ifndef MATERIAL_DGM_HH_
#define MATERIAL_DGM_HH_

#include "Material.hh"

namespace detran
{

class MaterialDGM : public Material
{

public:

  typedef Material Base;
  typedef SP<MaterialDGM> SP_material;

  MaterialDGM(int g, int m, int a, bool downscatter);

  static SP<Material> Create(int number_groups,
                                int number_materials,
                                int number_angles,
                                bool downscatter)
  {
    SP_material p;
    p = new MaterialDGM(number_groups, number_materials, number_angles, downscatter);
    return p;
  }

  double delta(int m, int g, int a) const;

  void set_delta(int m, int g, int a, double v);

private:

  using Base::b_number_groups;
  using Base::b_number_materials;

  // Number of angles
  int d_number_angles;

  // Delta [material, group, angle]
  vec3_dbl d_delta;

};

} // end namespace detran

#include "MaterialDGM.i.hh"

#endif /* MATERIAL_DGM_HH_ */
