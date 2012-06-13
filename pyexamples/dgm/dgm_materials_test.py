#  denovo/src/dgm/examples/simple_box_core/dgm_materials_test.py

# This script tests construction of a dgm materials object

import os, sys, math, string
from sc_2d import *
import XSMoments
import dlp

initialize(sys.argv)

dgm = False

#---------------------------------------------------------------------------##
# INPUT
#---------------------------------------------------------------------------##
db = DB("dgm_materials_test") 
db.insert("num_groups",           2)
db.insert("num_blocks_i",         1)
db.insert("num_blocks_j",         1)
db.insert("dgm",                  dgm)
# boundary conditions 
db.insert("boundary",             "reflect")
db.add_db("boundary_db",          "bnd_conditions")
db.insert("boundary_db",          "reflect", [1, 1, 1, 1], 1)
# solve parameters
db.insert("downscatter",          0, 1)
db.insert("Pn_order",             0)
db.insert("Pn_order_flux",        0)
# Solver information
db.insert("problem_type",         "EIGENVALUE")
db.insert("eigen_solver",         "power_iteration")
db.insert("within_group_solver",  "SI")
db.insert("mg_solver",            "gauss_seidel")
db.insert("tolerance",            1.0e-9)
db.insert("max_itr",              10000)
# Eigenvalue information
db.add_db("eigenvalue_db",        "eigenvalue")
db.insert("eigenvalue_db",        "k_tolerance",      1.0e-12)
db.insert("eigenvalue_db",        "max_tolerance",    0.001)
db.insert("eigenvalue_db",        "diagnostic_level", 2)
db.insert("eigenvalue_db",        "keff", 1.0)
# Upscatter database
db.add_db("upscatter_db",         "upscatter")
db.insert("upscatter_db",         "tolerance", 1.0e-3)
# Angular options
db.add_db("quadrature_db",        "quadrature")
db.insert("quadrature_db",        "Sn_order",     2)
db.insert("quadrature_db",        "quad_type",    "levelsym")
# Energy decomposition
#db.insert("num_sets",             1)
#db.insert("partition_upscatter",  1, 1)
## Dimension of Problem
db.insert("dimension",            2)

#-----------------------------------------------------------------------------#
# MESH
#-----------------------------------------------------------------------------#

length  = 10.0
nx      = 2
dx      = length / nx
x = []
for n in range(0, nx+1):
	x.append(dx * n)
db.insert("x_edges", x)
db.insert("y_edges", x)

#-----------------------------------------------------------------------------#
# MANAGER
#-----------------------------------------------------------------------------#

manager = Manager()                             # 
mat     = Mat(dgm)                              # TRUE for DGM
source  = Zero_Source()                         # zero source still needed?
angles  = Angles()                              #
manager.partition_2d(db, mat, angles)

#-----------------------------------------------------------------------------#
# MATERIAL
#-----------------------------------------------------------------------------#

# kinf = 1.2782865117e+00
# Total cross section 
SigmaT = [0.2263, 1.0119] 
# Scattering cross section -- DOWNSCATTER ONLY
        #   g<--0  g<--1        
SigmaS = [ [[0.2006]              ], \
           [[0.0161],[0.9355]], \
         ]
# Fission cross section
nuSigmaF = [0.0067, 0.1241]
# Fission spectrum
chi = [1.0, 0.0]
# One material
mat.set_num(1)
R1 = SigmaT[0] - SigmaS[0][0][0]
R2 = SigmaT[1] - SigmaS[1][1][0]
F1 = nuSigmaF[0]
F2 = nuSigmaF[1]
S12 = SigmaS[1][0][0]

print "kinf", (F1 + F2*S12/R2)/R1


if dgm :
  for g in range(0, 2) :
    #    void create_xs_dgm(int matid,
    #                       int g,
    #                       int num_energy_moments,
    #                       int num_angles,
    #                       int num_scattering_moments);
    #    void assign_total_dgm(int       matid,
    #                          int       g,
    #                          double    sigma_g);
    #    void assign_delta_dgm(int       matid,
    #                          int       g,
    #                          int       order,
    #                          int       angle,
    #                          double    delta);
    #    void assign_scatter_dgm(int    matid,
    #                            int    g,
    #                            int    order,
    #                            int    gp,
    #                            int    lm,
    #                            double sigma);
    #    void assign_nu_fission_dgm(int     matid,
    #                               int     g,
    #                               double  sigma);
    #    void assign_chi_dgm(int     matid,
    #                        int     g,
    #                        int     order,
    #                        double  chi);
    #    void complete_xs_dgm(int matid, int g);
    mat.create_xs_dgm(0, g, 1, 4, 1)
    mat.assign_total_dgm(0, g, SigmaT[g])
    for a in range(0, 4) :
      mat.assign_delta_dgm(0, g, 0, a, 0.0)
    for gp in range(0, g+1) :
      mat.assign_scatter_dgm(0, g, 0, gp, 0, SigmaS[g][gp][0])
    mat.assign_nu_fission_dgm(0, g, nuSigmaF[g])
    mat.assign_chi_dgm(0, g, 0, chi[g])
    mat.complete_xs_dgm(0, g)
else :
  for g in range(0, 2):
    mat.assign_xs(0, g, SigmaT[g], SigmaS[g])
  mat.assign_fission(0, nuSigmaF, chi)            # assign fission data



num_cells = (len(x)-1)*(len(x)-1)               # number of cells
ids = Vec_Int(num_cells, 0)                     # all cells get 
mat.assign_global_ids(ids)                      #   material 0

#-----------------------------------------------------------------------------#
# Setup, Verify, and Solve.
#-----------------------------------------------------------------------------#

manager.partition_energy(mat, angles) 
manager.setup(source)
manager.verify()
manager.solve(angles)
manager.close()

