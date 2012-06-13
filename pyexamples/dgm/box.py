# pyexamples/box.py
#
# A simple 2-d square region, 1 or 7 group reactor.  The goal here is
# to provide a very simple benchmark for the Denovo-DGM effort.

import numpy as np
import time
from detran import *
from dgm_utilities import *
import cPickle as pickle
import box_input

#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#

inp = InputDB.Create()
box_input.fill_input(inp)

#-----------------------------------------------------------------------------#
# MATERIAL
#-----------------------------------------------------------------------------#

if inp.get_int("number_groups") == 1 :

  mat = Material.Create(1, 1, False)
  mat.set_sigma_t(0, 0,    1.0)
  mat.set_sigma_s(0, 0, 0, 0.5)
  mat.set_nu_sigma_f(0, 0, 0.5)
  mat.set_chi(0, 0,        1.0)

else :

  # Material
  mat = Material.Create(7, 2, True)
  m = 0
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.77949E-01)
  mat.set_sigma_t(m, 1, 3.29805E-01)
  mat.set_sigma_t(m, 2, 4.80388E-01)
  mat.set_sigma_t(m, 3, 5.54367E-01)
  mat.set_sigma_t(m, 4, 3.11801E-01)
  mat.set_sigma_t(m, 5, 3.95168E-01)
  mat.set_sigma_t(m, 6, 5.64406E-01)
  # Fission times nu
  mat.set_sigma_f(m, 0, 7.21206E-03)
  mat.set_sigma_f(m, 1, 8.19301E-04)
  mat.set_sigma_f(m, 2, 6.45320E-03)
  mat.set_sigma_f(m, 3, 1.85648E-02)
  mat.set_sigma_f(m, 4, 1.78084E-02)
  mat.set_sigma_f(m, 5, 8.30348E-02)
  mat.set_sigma_f(m, 6, 2.16004E-01)
  mat.set_nu(m, 0, 2.78145E+00)
  mat.set_nu(m, 1, 2.47443E+00)
  mat.set_nu(m, 2, 2.43383E+00)
  mat.set_nu(m, 3, 2.43380E+00)
  mat.set_nu(m, 4, 2.43380E+00)
  mat.set_nu(m, 5, 2.43380E+00)
  mat.set_nu(m, 6, 2.43380E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 0 <- g'
  mat.set_sigma_s(m, 0, 0, 1.27537E-01)
  # 1 <- g'
  mat.set_sigma_s(m, 1, 0, 4.23780E-02)
  mat.set_sigma_s(m, 1, 1, 3.24456E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 2, 0, 9.43740E-06)
  mat.set_sigma_s(m, 2, 1, 1.63140E-03)
  mat.set_sigma_s(m, 2, 2, 4.50940E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 3, 0, 5.51630E-09)
  mat.set_sigma_s(m, 3, 1, 3.14270E-09)
  mat.set_sigma_s(m, 3, 2, 2.67920E-03)
  mat.set_sigma_s(m, 3, 3, 4.52565E-01)
  mat.set_sigma_s(m, 3, 4, 1.25250E-04*0.0)
  # 4 <- g'
  mat.set_sigma_s(m, 4, 3, 5.56640E-03)
  mat.set_sigma_s(m, 4, 4, 2.71401E-01)
  mat.set_sigma_s(m, 4, 5, 1.29680E-03*0.0)
  # 5 <- g'
  mat.set_sigma_s(m, 5, 4, 1.02550E-02)
  mat.set_sigma_s(m, 5, 5, 2.65802E-01)
  mat.set_sigma_s(m, 5, 6, 8.54580E-03*0.0)
  # 6 <- g'
  mat.set_sigma_s(m, 6, 4, 1.00210E-08)
  mat.set_sigma_s(m, 6, 5, 1.68090E-02)
  mat.set_sigma_s(m, 6, 6, 2.73080E-01)
  # --------------------------------------------
  # Material 1: moderator
  # --------------------------------------------
  m = 1
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.59206E-01)
  mat.set_sigma_t(m, 1, 4.12970E-01)
  mat.set_sigma_t(m, 2, 5.90310E-01)
  mat.set_sigma_t(m, 3, 5.84350E-01)
  mat.set_sigma_t(m, 4, 7.18000E-01)
  mat.set_sigma_t(m, 5, 1.25445E+00)
  mat.set_sigma_t(m, 6, 2.65038E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 4.44777E-02)
  # 2 <- g
  mat.set_sigma_s(m, 1, 0, 1.13400E-01)
  mat.set_sigma_s(m, 1, 1, 2.82334E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 7.23470E-04)
  mat.set_sigma_s(m, 2, 1, 1.29940E-01)
  mat.set_sigma_s(m, 2, 2, 3.45256E-01)
  # 4 <- g'
  aa = 0.0
  mat.set_sigma_s(m, 3, 0, 3.74990E-06)
  mat.set_sigma_s(m, 3, 1, 6.23400E-04)
  mat.set_sigma_s(m, 3, 2, 2.24570E-01)
  mat.set_sigma_s(m, 3, 3, 9.10284E-02)
  mat.set_sigma_s(m, 3, 4, 7.14370E-05*aa)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 0, 5.31840E-08)
  mat.set_sigma_s(m, 4, 1, 4.80020E-05)
  mat.set_sigma_s(m, 4, 2, 1.69990E-02)
  mat.set_sigma_s(m, 4, 3, 4.15510E-01)
  mat.set_sigma_s(m, 4, 4, 1.39138E-01)
  mat.set_sigma_s(m, 4, 5, 2.21570E-03*aa)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 1, 7.44860E-06)
  mat.set_sigma_s(m, 5, 2, 2.64430E-03)
  mat.set_sigma_s(m, 5, 3, 6.37320E-02)
  mat.set_sigma_s(m, 5, 4, 5.11820E-01)
  mat.set_sigma_s(m, 5, 5, 6.99913E-01)
  mat.set_sigma_s(m, 5, 6, 1.32440E-01*aa)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 1, 1.04550E-06)
  mat.set_sigma_s(m, 6, 2, 5.03440E-04)
  mat.set_sigma_s(m, 6, 3, 1.21390E-02)
  mat.set_sigma_s(m, 6, 4, 6.12290E-02)
  mat.set_sigma_s(m, 6, 5, 5.37320E-01)
  mat.set_sigma_s(m, 6, 6, 2.48070E+00)
  
mat.finalize()

#-----------------------------------------------------------------------------#
# MESH
#-----------------------------------------------------------------------------#

cm, fm, mt = box_input.get_geom()
mesh = Mesh2D.Create(fm, fm, cm, cm, mt)

#-----------------------------------------------------------------------------#
# QUADRATURE
#-----------------------------------------------------------------------------#

quad_order = 8
quad = LevelSymmetric.Create(quad_order, 2)

#-----------------------------------------------------------------------------#
# MISC. SETUP
#-----------------------------------------------------------------------------#

# State
state = State.Create(inp, mesh, quad)
# Empty external source.
q_e = ExternalSourceSP()
# Fission source
q_f = FissionSource.Create(state, mesh, mat)
q_f.initialize()
# Boundary
bound = Boundary2D.Create(inp, mesh, quad)

#-----------------------------------------------------------------------------#
# SOLVE
#-----------------------------------------------------------------------------#

solver = PowerIteration2D.Create(inp, state, mesh, mat, quad, bound, q_e, q_f)

#start = time.time()
solver.solve()
#elapsed = (time.time() - start)

#print elapsed, " seconds"

#phi_0 = np.asarray(state.phi(0))

#psi_0= np.asarray(state.psi(0, 0, 0))
#psi_1= np.asarray(state.psi(0, 0, 1))
#psi_2= np.asarray(state.psi(0, 0, 2))
#psi_3= np.asarray(state.psi(0, 0, 3))

#for i in range(0, 16):
#  print '%10.8f %10.8f %10.8f %10.8f %10.8f' % (phi_0[i]/phi_0[0], psi_0[i]/psi_0[0], psi_1[i]/psi_0[0], psi_2[i]/psi_0[0], psi_3[i]/psi_0[0])

# DGM tests

from dlp import *
coarse_structure = [3,4]
dlp, rhos = get_dlop(coarse_structure)


phi = np.zeros((7, mesh.number_cells()))
psi = np.zeros((7, quad.number_angles(), mesh.number_cells()))
for g in range(0, 7) :
  phi[g] = np.asarray(state.phi(g))
  angle = 0
  for o in range(0, 4) :
    for a in range(0, quad.number_angles_octant()) :
      psi[g][angle] = np.asarray(state.psi(g, o, a))
      angle = angle + 1
  #print mat.sigma_t(0, g), mat.nu_sigma_f(0, g), phi[g][0]

# Get DGM moments
sigma_t, nu_sigma_f = calculate_sigma(phi, mat, mesh, coarse_structure)
chi = calculate_chi(dlp, mat, mesh, coarse_structure)
sigma_s = calculate_sigma_s(phi, dlp, mat, mesh, coarse_structure)
delta = calculate_delta(phi, psi, dlp, mat, mesh, quad, sigma_t, coarse_structure)

dgm_data = {}
dgm_data['quad_order']       = quad_order
dgm_data['coarse_structure'] = coarse_structure
dgm_data['sigma_t']     = sigma_t
dgm_data['nu_sigma_f']  = nu_sigma_f
dgm_data['chi']         = chi
dgm_data['sigma_s']     = sigma_s
dgm_data['delta']       = delta
pickle.dump(dgm_data, open('dgm_data_ref.p', 'wb'))

#print ""
#print "sigma_t"
#for i in range(0, 16) :
#  print format(sigma_t[0][i], '10.8f'), format(sigma_t[1][i], '10.8f')
#print "sigma_f"
#for i in range(0, 16) :
#  print format(sigma_t[0][i], '10.8f'), format(sigma_t[1][i], '10.8f')
#print "delta"

#for cg in range(0, len(coarse_structure)) :
#  print "cg = ", cg
#  for i in range(0, 1) :
#    #print " i = ", i
#    for a in range(0, quad.number_angles()) :
#      #print " a = ", a
#      for c in range(0, 1) :
#        print format(delta[cg][i][a][c], '10.8f')

