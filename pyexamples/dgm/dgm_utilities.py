# dgm_utilities
#
# A collection of post-process utilities for generating DGM
# data, etc.
import numpy as np

def calculate_sigma(phi, mat, mesh, coarse_structure) :
  """ Calculate the DGM total and nu*fission cross-section.

  Follows Eqs. 17 and 28 in the NSE paper.
  """
  ncg = len(coarse_structure)
  sigma_t    = np.zeros((ncg, mesh.number_cells()))
  nu_sigma_f = np.zeros((ncg, mesh.number_cells()))

  mm = mesh.mesh_map("MATERIAL")

  # Lower fine group index
  lower_g = 0

  # Loop over coarse groups
  for cg in range(0, ncg) :

    # Upper fine group index
    upper_g = lower_g + coarse_structure[cg]

    # Loop over cells
    for cell in xrange(0, mesh.number_cells()) :
      num1 = 0.0
      num2 = 0.0
      den  = 0.0
      
      # Loop over fine groups
      for fg in range(lower_g, upper_g) :
        num1 += mat.sigma_t(mm[cell], fg) * phi[fg][cell]
        num2 += mat.nu_sigma_f(mm[cell], fg) * phi[fg][cell]
        den  += phi[fg][cell]
      sigma_t[cg][cell]    = num1 / den
      nu_sigma_f[cg][cell] = num2 / den

    lower_g = upper_g
    
  return sigma_t, nu_sigma_f

def calculate_chi(dlp, mat, mesh, coarse_structure) :
  """ Calculate the DGM fission spectrum

  Follows Eqs. 14 in the NSE paper.
  """
  ncg = len(coarse_structure)
  # list of coarse group chi moments
  chi = []
  mm = mesh.mesh_map("MATERIAL")
  lower_K = 0
  # loop over coarse groups
  for cg, coarse in enumerate(coarse_structure):
    upper_K = lower_K + coarse
    chi_cg = np.zeros((coarse, mesh.number_cells()))
    # loop over fine groups
    for i in xrange(lower_K, upper_K):
      # loop over space
      for cell in xrange(0, mesh.number_cells()) :
        val = 0.0
        # loop over fine groups
        for K in xrange(lower_K,upper_K):
          val += dlp[(K-lower_K,i-lower_K, cg)]*mat.chi(mm[cell], K)
        chi_cg[i - lower_K][cell] = val
    chi.append(chi_cg)
    lower_K = upper_K
  return chi

def calculate_sigma_s(phi,dlp,mat,mesh,coarse_structure):
  """
        [g, i, gp, cell]
        sigma_si = [g_prime1, g_prime2, g_prime3]
                    g_prime1 = [sigma_si_cell1, sigma_si_cell2, ... ]
  g - coarse group index (from 0 to len(coarse_structure))
  i - moment order index within chosen coarse group G (zero is leading order) (from 0 to coarse_structure[g])
  from eqn 29 of NSE paper, using only scalar flux for isotropic scattering.
  L and K are global finegroup group numbers here, whereas in eqn 29 they're relative within the coarse group
  """
  sigma_s = []
  mm = mesh.mesh_map("MATERIAL")
  # number coarse groups
  ncg = len(coarse_structure)

  for g in range(0, ncg) :

    # number group moments
    ngm = coarse_structure[g]
    
    sigma_sg = np.zeros( (ngm, ncg, mesh.number_cells()) )

    for i in range(0, ngm) :
      lower_K = 0
      if not g == 0:
        for h in range(g):
          lower_K += coarse_structure[h]
      upper_K = lower_K+coarse_structure[g]
      lower_L = 0
      for g_prime, coarse in enumerate(coarse_structure):
        upper_L = lower_L+coarse
        cellVals = []
        for cell in xrange(0, mesh.number_cells()):
          num = 0.0
          denom = 0.0
          # K <-- L
          for L in xrange(lower_L, upper_L):
            for K in xrange(lower_K, upper_K):
              if K >= L :
                num += dlp[(K-lower_K,i,g)] * mat.sigma_s(mm[cell],K,L)*phi[L][cell]
            denom += phi[L][cell]
          #print i, g_prime, cell, np.size(sigma_sg, 0), np.size(sigma_sg, 1), np.size(sigma_sg, 2)
          sigma_sg[i][g_prime][cell] = num/denom
        lower_L = upper_L
    sigma_s.append(sigma_sg)

  return sigma_s

def calculate_delta(phi,psi,dlp,mat,mesh,quad,sigma_t,coarse_structure):
  """
        [g, i, gp, cell]
        sigma_si = [g_prime1, g_prime2, g_prime3]
                    g_prime1 = [sigma_si_cell1, sigma_si_cell2, ... ]
  g - coarse group index (from 0 to len(coarse_structure))
  i - moment order index within chosen coarse group G (zero is leading order) (from 0 to coarse_structure[g])
  from eqn 29 of NSE paper, using only scalar flux for isotropic scattering.
  L and K are global finegroup group numbers here, whereas in eqn 29 they're relative within the coarse group
  """
  delta = []
  # material map
  mm = mesh.mesh_map("MATERIAL")
  # number coarse groups
  ncg = len(coarse_structure)
  # number angles
  na = quad.number_angles()
  # number cells
  nc = mesh.number_cells()

  # loop over coarse groups
  for g in range(0, ncg) :

    # number group moments
    ngm = coarse_structure[g]
    
    delta_g = np.zeros( (ngm, na, nc) )

    # loop over group moments
    for i in range(0, ngm) :

      lower_K = 0
      if not g == 0:
          for h in range(g):
              lower_K += coarse_structure[h]
      upper_K = lower_K+coarse_structure[g]

      for a in xrange(0, na):
        for cell in xrange(0, nc) :
          num = 0
          denom = 0
          for K in xrange(lower_K, upper_K):
            sigma_g = mat.sigma_t(mm[cell], K) - sigma_t[g][cell]
            num    += dlp[(K-lower_K,i,g)] * sigma_g*psi[K][a][cell]
            denom  += phi[K][cell]
          delta_g[i][a][cell] = num/denom

    delta.append(delta_g)
  return delta

  return delta
