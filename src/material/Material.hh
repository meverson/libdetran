//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Material.hh
 * \author Jeremy Roberts
 * \brief  Material class definition.
 */
//---------------------------------------------------------------------------//


#ifndef MATERIAL_HH_
#define MATERIAL_HH_

// Detran
#include "Definitions.hh"
#include "DBC.hh"
#include "SP.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Material
 * \brief Simple cross section container.
 */
//---------------------------------------------------------------------------//
class Material : public Object
{

public:

  /// \name Useful typedefs when using Material
  /// \{
  typedef SP<Material> SP_material;
  /// \}

  /*!
   *  \brief Constructor.
   *
   *  \param    number_groups       Number of energy groups.
   *  \param    number_materials    Number of materials.
   *  \param    downscatter         Switch on to use only downscatter.
   */
  Material(int  number_groups,
           int  number_materials,
           bool downscatter = false);

  virtual ~Material(){}

  /*!
   *  \brief SP Constructor.
   *
   *  \param    number_groups       Number of energy groups.
   *  \param    number_materials    Number of materials.
   *  \param    downscatter         Switch on to use only downscatter.
   *  \return                       Smart pointer to Material object.
   */
  static SP_material Create(int number_groups,
                            int number_materials,
                            bool downscatter)
  {
    SP_material p;
    p = new Material(number_groups, number_materials, downscatter);
    return p;
  }

  //--------------------------------------------------------------------------//
  // Setters
  //--------------------------------------------------------------------------//

  void set_sigma_t(int m, int g, double v);

  void set_sigma_a(int m, int g, double v);

  void set_nu_sigma_f(int m, int g, double v);

  void set_sigma_f(int m, int g, double v);

  void set_nu(int m, int g, double v);

  void set_chi(int m, int g, double v);

  void set_sigma_s(int m, int g, int gp, double v);

  void set_diff_coef(int m, int g, double v);

  void set_delta(int m, int g, int a, double v);

  // Vectorized setters

  void set_sigma_t(int m, vec_dbl &v);

  void set_sigma_a(int m, vec_dbl &v);

  void set_nu_sigma_f(int m, vec_dbl &v);

  void set_sigma_f(int m, vec_dbl &v);

  void set_nu(int m, vec_dbl &v);

  void set_chi(int m, vec_dbl &v);

  void set_sigma_s(int m, int g, vec_dbl &v);

  void set_diff_coef(int m, vec_dbl &v);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  inline double sigma_t(int m, int g) const;

  inline double sigma_a(int m, int g) const;

  inline double nu_sigma_f(int m, int g) const;

  inline double sigma_f(int m, int g) const;

  inline double nu(int m, int g) const;

  inline double chi(int m, int g) const;

  inline double sigma_s(int m, int g, int gp) const;

  inline double diff_coef(int m, int g) const;

  inline double delta(int m, int g, int a) const;

  int number_groups()
  {
    return b_number_groups;
  }

  int number_materials()
  {
    return b_number_materials;
  }

  /*!
   *  \brief Lower scatter group bound.
   *
   *  This is the *lowest* index (highest energy) \f$ g' \f$
   *  that leads to downscatter for a given outgoing group \f$ g \f$.
   */
  int lower(int g);

  /*!
   *  \brief Upper scatter group bound.
   *
   *  This is the *highest* index (lowest energy) \f$ g' \f$
   *  that upscatters into the outgoing group \f$ g \f$.
   */
  int upper(int g);

  /// Do we do only downscatter?
  bool downscatter() const
  {
    return b_downscatter;
  }

  /// Index below which upscatter doesn't occur for any material.
  int upscatter_cutoff() const
  {
    return b_upscatter_cutoff;
  }

  void set_dgm(int number_angles)
  {
    b_number_angles = number_angles;
    d_delta.resize(b_number_materials,
                   vec2_dbl(b_number_groups,
                            vec_dbl(b_number_angles, 0.0)));
    b_dgm = true;
  }

  bool has_dgm() const
  {
    return b_dgm;
  }

  /// Computes scattering bounds and absorption cross section.
  void finalize();

  /// Pretty print the material database.
  void display();

  /// Incomplete implementation of DBC function.
  bool is_valid() const
  {
    Ensure(b_number_groups >= 0);
    Ensure(b_finalized);
    return true;
  };

protected:

  /// Number of groups
  int b_number_groups;

  /// Number of materials
  int b_number_materials;

  /// Downscatter switch (when true, upscatter ignored)
  bool b_downscatter;

  /// Total cross section [material, group]
  vec2_dbl b_sigma_t;

  /// Absorption cross section [material, group]
  vec2_dbl b_sigma_a;

  /// nu * Fission [material, group]
  vec2_dbl b_nu_sigma_f;

  /// Fission [material, group]
  vec2_dbl b_sigma_f;

  /// nu [material, group]
  vec2_dbl b_nu;

  /// Fission spectrum [material, group]
  vec2_dbl b_chi;

  /// Scatter [material, group<-, group']
  vec3_dbl b_sigma_s;

  /// Diffusion coefficient [material, group]
  vec2_dbl b_diff_coef;

  /// Scatter bounds applied to all materials [group, 2]
  vec2_int b_scatter_bounds;

  /*!
   *  Upscatter cutoff.  Only groups equal to or above this cutoff are
   *  subject to upscatter iterations.
   */
  int b_upscatter_cutoff;

  /// Are we ready to be used?
  bool b_finalized;

  /// DGM

  bool b_dgm;

  int b_number_angles;

  // Delta [material, group, angle]
  vec3_dbl d_delta;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Material.i.hh"

#endif /* MATERIAL_HH_ */
