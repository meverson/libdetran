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

data = pickle.load(open('dgm_data_ref.p', 'rb'))
quad_order = data['quad_order']
coarse_structure = data['coarse_structure']
number_groups = len(coarse_structure)

#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#

inp = InputDB.Create()
box_input.fill_input(inp)
inp.put_int("number_groups",      number_groups)

#-----------------------------------------------------------------------------#
# QUADRATURE
#-----------------------------------------------------------------------------#

quad = LevelSymmetric.Create(quad_order, 2)

#-----------------------------------------------------------------------------#
# MESH
#-----------------------------------------------------------------------------#

cm, fm, mt = box_input.get_geom()
L = max(cm)

cm = np.linspace(0, L, sum(fm)+1)
mt = range(0, sum(fm)*sum(fm))
print cm
mesh = Mesh2D.Create(cm, cm, mt)

#-----------------------------------------------------------------------------#
# MATERIAL
#-----------------------------------------------------------------------------#

mat = Material.Create(number_groups, mesh.number_cells(), False)
mat.set_dgm(quad.number_angles())


for g in range(0, number_groups) :
  for m in range(0, mesh.number_cells()) :
    mat.set_sigma_t( m, g, data['sigma_t'][g][m] )
    #print data['sigma_t'][g][m]
    mat.set_sigma_f( m, g, data['nu_sigma_f'][g][m] )
    
    for gp in range(0, g + 1) :
      mat.set_sigma_s( m, g, gp, data['sigma_s'][g][0][gp][m] )
    for a in range(0, quad.number_angles()) :
      mat.set_delta( m, g, a, data['delta'][g][0][a][m] )

    mat.set_chi( m, g, data['chi'][g][0][m] )

mat.finalize()
#exit()

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
solver.solve()
#print mat.sigma_t(0, 0), mat.sigma_t(1, 0)
#print np.asarray(mesh.mesh_map("MATERIAL"))

#v0 = np.asarray(state.phi(0))
#print v0
#v1 = np.asarray(state.phi(1))
#print v1
